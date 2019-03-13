use bio::io::fastq;
use bio::io::fastq::{FastqRead, Record};
use csv;
use ordered_float::NotNaN;
use rocksdb::DB;
use serde_json;
use std::error::Error;
use std::io;
use std::io::Write;
use std::mem;
use std::process::{Command, Stdio};
use std::str;
use tempfile::tempdir;
use uuid::Uuid;

use super::calc_consensus::{calc_consensus, calc_paired_consensus};

const HAMMING_THRESHOLD: f64 = 10.0;

/// Interpret a cluster returned by starcode
fn parse_cluster(record: csv::StringRecord) -> Result<Vec<usize>, Box<dyn Error>> {
    let seqids = &record[2];
    Ok(csv::ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(false)
        .from_reader(seqids.as_bytes())
        .deserialize()
        .next()
        .unwrap()?)
}

/// Calculates the median hamming distance for all records by deriving the overlap from insert size
fn median_hamming_distance(
    insert_size: &usize,
    f_recs: &Vec<fastq::Record>,
    r_recs: &Vec<fastq::Record>,
) -> Option<f64> {
    let mut distances: Vec<f64> = Vec::new();

    for (f_rec, r_rec) in f_recs.into_iter().zip(r_recs).map(|(a, b)| (a, b)) {
        // check if reads overlap within insert size
        if (insert_size < &f_rec.seq().len()) | (insert_size < &r_rec.seq().len()) {
            return None;
        }
        if insert_size >= &(&f_rec.seq().len() + &r_rec.seq().len()) {
            return None;
        }
        let overlap = (f_rec.seq().len() + r_rec.seq().len()) - insert_size;
        let suffix_start_idx: usize = f_rec.seq().len() - overlap;
        let mut distance = bio::alignment::distance::hamming(
            &f_rec.seq()[suffix_start_idx..],
            &bio::alphabets::dna::revcomp(r_rec.seq())[..overlap],
        ) as f64;
        distance = distance / overlap as f64; //Use absolute or relative distance?!
        distances.push(distance);
    }
    stats::median(distances.into_iter())
}

/// Used to store a mapping of read index to read sequence
#[derive(Debug)]
struct FASTQStorage {
    db: DB,
    storage_dir: std::path::PathBuf,
}

impl FASTQStorage {
    /// Create a new FASTQStorage using a Rocksdb database
    /// that maps read indices to read seqeunces.
    pub fn new() -> Result<Self, Box<Error>> {
        // Save storage_dir to prevent it from leaving scope and
        // in turn deleting the tempdir
        let storage_dir = tempdir()?.path().join("db");
        Ok(FASTQStorage {
            db: DB::open_default(storage_dir.clone())?,
            storage_dir,
        })
    }

    fn as_key<'a>(i: u64) -> [u8; 8] {
        unsafe { mem::transmute::<u64, [u8; 8]>(i) }
    }

    /// Enter a (read index, read sequence) pair into the database.
    pub fn put(
        &mut self,
        i: usize,
        f_rec: &fastq::Record,
        r_rec: &fastq::Record,
    ) -> Result<(), Box<dyn Error>> {
        Ok(self.db.put(
            &Self::as_key(i as u64),
            serde_json::to_string(&(f_rec, r_rec))?.as_bytes(),
        )?)
    }

    /// Retrieve the read sequence of the read with index `i`.
    pub fn get(&self, i: usize) -> Result<(fastq::Record, fastq::Record), Box<dyn Error>> {
        Ok(serde_json::from_str(
            str::from_utf8(&self.db.get(&Self::as_key(i as u64))?.unwrap()).unwrap(),
        )?)
    }
}

pub trait CallConsensusReads<'a, R: io::Read + 'a, W: io::Write + 'a> {
    /// Cluster reads from fastq readers according to their sequence
    /// and UMI, then compute a consensus sequence.
    ///
    /// Cluster the reads in the input file according to their sequence
    /// (concatenated p5 and p7 reads without UMI). Read the
    /// identified clusters, and cluster all reds in a cluster by UMI,
    /// creating groups of very likely PCR duplicates.
    /// Next, compute a consensus read for each unique read,
    /// i.e. a cluster with similar sequences and identical UMI,
    /// and write it into the output files.
    fn call_consensus_reads(&'a mut self) -> Result<(), Box<dyn Error>> {
        // cluster by sequence
        // Note: If starcode is not installed, this throws a
        // hard to interpret error:
        // (No such file or directory (os error 2))
        // The expect added below should make this more clear.
        // cluster within in this cluster by umi
        let mut umi_cluster = Command::new("starcode")
            .arg("--dist")
            .arg(format!("{}", self.umi_dist()))
            .arg("--seq-id")
            .arg("-s")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("Error in starcode call. Starcode might not be installed.");

        let mut f_rec = fastq::Record::new();
        let mut r_rec = fastq::Record::new();
        // init temp storage for reads
        let mut read_storage = FASTQStorage::new()?;
        let mut i = 0;
        let mut seqs = Vec::new();
        loop {
            self.fq1_reader().read(&mut f_rec)?;
            self.fq2_reader().read(&mut r_rec)?;

            match (f_rec.is_empty(), r_rec.is_empty()) {
                (true, true) => break,
                (false, false) => (),
                _ => panic!("Given FASTQ files have unequal lengths"),
            }
            // save umis for second (intra cluster) clustering
            let umi = if self.reverse_umi() {
                r_rec.seq()[..self.umi_len()].to_owned()
            } else {
                f_rec.seq()[..self.umi_len()].to_owned()
            };
            umi_cluster.stdin.as_mut().unwrap().write(&umi)?;
            umi_cluster.stdin.as_mut().unwrap().write(b"\n")?;

            read_storage.put(i, &f_rec, &r_rec)?;
            let seq = if self.reverse_umi() {
                [&f_rec.seq()[..], &r_rec.seq()[self.umi_len()..]].concat()
            } else {
                [&f_rec.seq()[self.umi_len()..], &r_rec.seq()[..]].concat()
            };
            seqs.push(seq);
            i += 1;
        }
        umi_cluster.stdin.as_mut().unwrap().flush()?;
        drop(umi_cluster.stdin.take());

        eprint!("Read starcode results");
        // read clusters identified by the first starcode run
        for record in csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(umi_cluster.stdout.as_mut().unwrap())
            .records()
        {
            let seqids = parse_cluster(record?)?;
            // cluster within in this cluster by umi
            let mut seq_cluster = Command::new("starcode")
                .arg("--dist")
                .arg(format!("{}", self.seq_dist()))
                .arg("--seq-id")
                .arg("-s")
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()?;
            for &seqid in &seqids {
                seq_cluster
                    .stdin
                    .as_mut()
                    .unwrap()
                    .write(&seqs[seqid - 1])?;
                seq_cluster.stdin.as_mut().unwrap().write(b"\n")?;
            }
            seq_cluster.stdin.as_mut().unwrap().flush()?;
            drop(seq_cluster.stdin.take());

            // handle each potential unique read
            for record in csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .from_reader(seq_cluster.stdout.as_mut().unwrap())
                .records()
            {
                let inner_seqids = parse_cluster(record?)?;
                // this is a proper cluster
                // calculate consensus reads and write to output FASTQs
                let mut f_recs = Vec::new();
                let mut r_recs = Vec::new();
                let mut outer_seqids = Vec::new();

                for inner_seqid in inner_seqids {
                    let seqid = seqids[inner_seqid - 1];
                    let (f_rec, r_rec) = read_storage.get(seqid - 1)?;
                    f_recs.push(f_rec);
                    r_recs.push(r_rec);
                    outer_seqids.push(seqid);
                }
                self.write_records(f_recs, r_recs, outer_seqids)?;
            }

            match seq_cluster
                .wait()
                .expect("process did not even start")
                .code()
            {
                Some(0) => (),
                Some(s) => println!("Starcode failed with error code {}", s),
                None => println!("Starcode was terminated by signal"),
            }
        }
        Ok(())
    }

    fn write_records(
        &mut self,
        f_recs: Vec<Record>,
        r_recs: Vec<Record>,
        outer_seqids: Vec<usize>,
    ) -> Result<(), Box<Error>>;
    fn fq1_reader(&mut self) -> &mut fastq::Reader<R>;
    fn fq2_reader(&mut self) -> &mut fastq::Reader<R>;
    fn umi_len(&self) -> usize;
    fn seq_dist(&self) -> usize;
    fn umi_dist(&self) -> usize;
    fn reverse_umi(&self) -> bool;
}

/// Struct for calling non-overlapping consensus reads
/// Implements CallConsensusReads-Trait
pub struct CallNonOverlappingConsensusRead<'a, R: io::Read, W: io::Write> {
    fq1_reader: &'a mut fastq::Reader<R>,
    fq2_reader: &'a mut fastq::Reader<R>,
    fq1_writer: &'a mut fastq::Writer<W>,
    fq2_writer: &'a mut fastq::Writer<W>,
    umi_len: usize,
    seq_dist: usize,
    umi_dist: usize,
    reverse_umi: bool,
}

impl<'a, R: io::Read, W: io::Write> CallNonOverlappingConsensusRead<'a, R, W> {
    pub fn new(
        fq1_reader: &'a mut fastq::Reader<R>,
        fq2_reader: &'a mut fastq::Reader<R>,
        fq1_writer: &'a mut fastq::Writer<W>,
        fq2_writer: &'a mut fastq::Writer<W>,
        umi_len: usize,
        seq_dist: usize,
        umi_dist: usize,
        reverse_umi: bool,
    ) -> Self {
        CallNonOverlappingConsensusRead {
            fq1_reader,
            fq2_reader,
            fq1_writer,
            fq2_writer,
            umi_len,
            seq_dist,
            umi_dist,
            reverse_umi,
        }
    }
}

impl<'a, R: io::Read, W: io::Write> CallConsensusReads<'a, R, W>
    for CallNonOverlappingConsensusRead<'a, R, W>
{
    fn write_records(
        &mut self,
        f_recs: Vec<Record>,
        r_recs: Vec<Record>,
        outer_seqids: Vec<usize>,
    ) -> Result<(), Box<dyn Error>> {
        if f_recs.len() > 1 {
            let uuid = &Uuid::new_v4().to_hyphenated().to_string();
            self.fq1_writer
                .write_record(&calc_consensus(&f_recs, &outer_seqids, uuid).0)?;
            self.fq2_writer
                .write_record(&calc_consensus(&r_recs, &outer_seqids, uuid).0)?;
        } else {
            self.fq1_writer.write_record(&f_recs[0])?;
            self.fq2_writer.write_record(&r_recs[0])?;
        }
        Ok(())
    }
    fn fq1_reader(&mut self) -> &mut fastq::Reader<R> {
        &mut self.fq1_reader
    }
    fn fq2_reader(&mut self) -> &mut fastq::Reader<R> {
        &mut self.fq2_reader
    }
    fn umi_len(&self) -> usize {
        self.umi_len
    }
    fn seq_dist(&self) -> usize {
        self.seq_dist
    }
    fn umi_dist(&self) -> usize {
        self.umi_dist
    }
    fn reverse_umi(&self) -> bool {
        self.reverse_umi
    }
}

///Clusters fastq reads by UMIs and calls consensus for overlapping reads
///
pub struct CallOverlappingConsensusRead<'a, R: io::Read, W: io::Write> {
    fq1_reader: &'a mut fastq::Reader<R>,
    fq2_reader: &'a mut fastq::Reader<R>,
    fq_writer: &'a mut fastq::Writer<W>,
    umi_len: usize,
    seq_dist: usize,
    umi_dist: usize,
    insert_size: usize,
    std_dev: usize,
    reverse_umi: bool,
}

impl<'a, R: io::Read, W: io::Write> CallOverlappingConsensusRead<'a, R, W> {
    pub fn new(
        fq1_reader: &'a mut fastq::Reader<R>,
        fq2_reader: &'a mut fastq::Reader<R>,
        fq_writer: &'a mut fastq::Writer<W>,
        umi_len: usize,
        seq_dist: usize,
        umi_dist: usize,
        insert_size: usize,
        std_dev: usize,
        reverse_umi: bool,
    ) -> Self {
        CallOverlappingConsensusRead {
            fq1_reader,
            fq2_reader,
            fq_writer,
            umi_len,
            seq_dist,
            umi_dist,
            insert_size,
            std_dev,
            reverse_umi,
        }
    }
}

impl<'a, R: io::Read, W: io::Write> CallConsensusReads<'a, R, W>
    for CallOverlappingConsensusRead<'a, R, W>
{
    fn write_records(
        &mut self,
        f_recs: Vec<Record>,
        r_recs: Vec<Record>,
        outer_seqids: Vec<usize>,
    ) -> Result<(), Box<Error>> {
        // Determine hamming distance for different insert sizes
        let mut median_distances: Vec<(f64, usize)> = Vec::new();
        for insert_size in
            (self.insert_size - 2 * self.std_dev)..(self.insert_size + 2 * self.std_dev)
        {
            let median_hamming_distance = median_hamming_distance(&insert_size, &f_recs, &r_recs);
            if let Some(median_hamming_distance) = median_hamming_distance {
                median_distances.push((median_hamming_distance, insert_size))
            }
        }

        //TODO Add deterministic uuid considering read ids
        let uuid = &Uuid::new_v4().to_hyphenated().to_string();
        if let Some(consensus_record) = median_distances
            .iter()
            .filter_map(|(mean_distance, insert_size)| {
                if *mean_distance < HAMMING_THRESHOLD {
                    let overlap = (f_recs[0].seq().len() + r_recs[0].seq().len()) - insert_size;
                    Some(calc_paired_consensus(
                        &f_recs,
                        &r_recs,
                        &overlap,
                        &outer_seqids,
                        &uuid,
                    ))
                } else {
                    None
                }
            })
            .max_by_key(|&(_, lh)| NotNaN::new(*lh).unwrap())
        {
            self.fq_writer.write_record(&consensus_record.0)?;
        } else {
            //Replace this by calling non overlapping consensus read in future implementation
            eprintln!("No read pairs with hamming distance < {} found in cluster! No consensus read created. If no read is created at all check insert size.", HAMMING_THRESHOLD);
        }
        Ok(())
    }
    fn fq1_reader(&mut self) -> &mut fastq::Reader<R> {
        &mut self.fq1_reader
    }
    fn fq2_reader(&mut self) -> &mut fastq::Reader<R> {
        &mut self.fq2_reader
    }
    fn umi_len(&self) -> usize {
        self.umi_len
    }
    fn seq_dist(&self) -> usize {
        self.seq_dist
    }
    fn umi_dist(&self) -> usize {
        self.umi_dist
    }
    fn reverse_umi(&self) -> bool {
        self.reverse_umi
    }
}
