use bio::io::fastq;
use bio::io::fastq::{FastqRead, Record};
use csv;
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

use super::calc_consensus::calc_consensus;

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

/// Used to store a mapping of read index to read sequence
#[derive(Debug)]
pub struct FASTQStorage {
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
        let mut seq_cluster = Command::new("starcode")
            .arg("--dist")
            .arg(format!("{}", self.seq_dist()))
            .arg("--seq-id")
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

        let mut umis = Vec::new();
        loop {
            self.fq1_reader().read(&mut f_rec)?;
            self.fq2_reader().read(&mut r_rec)?;

            match (f_rec.is_empty(), r_rec.is_empty()) {
                (true, true) => break,
                (false, false) => (),
                _ => panic!("Given FASTQ files have unequal lengths"),
            }
            // save umis for second (intra cluster) clustering
            let umi = r_rec.seq()[..self.umi_len()].to_owned();
            umis.push(umi);

            read_storage.put(i, &f_rec, &r_rec)?;

            let seq = [f_rec.seq(), &r_rec.seq()[self.umi_len()..]].concat();
            seq_cluster.stdin.as_mut().unwrap().write(&seq)?;
            seq_cluster.stdin.as_mut().unwrap().write(b"\n")?;
            i += 1;
        }

        seq_cluster.stdin.as_mut().unwrap().flush()?;
        drop(seq_cluster.stdin.take());

        eprint!("Read starcode results");
        // read clusters identified by the first starcode run
        for record in csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(seq_cluster.stdout.as_mut().unwrap())
            .records()
        {
            let seqids = parse_cluster(record?)?;
            // cluster within in this cluster by umi
            let mut umi_cluster = Command::new("starcode")
                .arg("--dist")
                .arg(format!("{}", self.umi_dist()))
                .arg("--seq-id")
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()?;
            for &seqid in &seqids {
                umi_cluster
                    .stdin
                    .as_mut()
                    .unwrap()
                    .write(&umis[seqid - 1])?;
                umi_cluster.stdin.as_mut().unwrap().write(b"\n")?;
            }
            umi_cluster.stdin.as_mut().unwrap().flush()?;
            drop(umi_cluster.stdin.take());

            // handle each potential unique read
            for record in csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .from_reader(umi_cluster.stdout.as_mut().unwrap())
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

            match umi_cluster
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
    ) -> Self {
        CallNonOverlappingConsensusRead {
            fq1_reader,
            fq2_reader,
            fq1_writer,
            fq2_writer,
            umi_len,
            seq_dist,
            umi_dist,
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
    ) -> Result<(), Box<Error>> {
        if f_recs.len() > 1 {
            let uuid = &Uuid::new_v4().to_hyphenated().to_string();
            self.fq1_writer
                .write_record(&calc_consensus(&f_recs, &outer_seqids, uuid))?;
            self.fq2_writer
                .write_record(&calc_consensus(&r_recs, &outer_seqids, uuid))?;
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
}
