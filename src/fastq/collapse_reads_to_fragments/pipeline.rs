use bio::io::fastq;
use bio::io::fastq::{FastqRead, Record};
use bio::stats::probs::LogProb;
use csv;
use derive_new::new;
use indicatif;
use ordered_float::NotNaN;
use rgsl::randist::gaussian::ugaussian_P;
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

use super::calc_consensus::{CalcNonOverlappingConsensus, CalcOverlappingConsensus};

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
    insert_size: usize,
    f_recs: &[fastq::Record],
    r_recs: &[fastq::Record],
) -> Option<f64> {
    let distances = f_recs.iter().zip(r_recs).filter_map(|(f_rec, r_rec)| {
        // check if reads overlap within insert size
        if (insert_size < f_rec.seq().len()) | (insert_size < r_rec.seq().len()) {
            return None;
        }
        if insert_size >= (&f_rec.seq().len() + &r_rec.seq().len()) {
            return None;
        }
        let overlap = (f_rec.seq().len() + r_rec.seq().len()) - insert_size;
        let suffix_start_idx: usize = f_rec.seq().len() - overlap;
        Some(bio::alignment::distance::hamming(
            &f_rec.seq()[suffix_start_idx..],
            &bio::alphabets::dna::revcomp(r_rec.seq())[..overlap],
        ))
    });
    stats::median(distances)
}

/// as shown in http://www.milefoot.com/math/stat/pdfc-normaldisc.htm
fn isize_pmf(value: f64, mean: f64, sd: f64) -> LogProb {
    LogProb((ugaussian_P((value + 0.5 - mean) / sd) - ugaussian_P((value - 0.5 - mean) / sd)).ln())
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
    pub fn new() -> Result<Self, Box<dyn Error>> {
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

pub struct OverlappingConsensus {
    record: Record,
    likelihood: LogProb,
}

pub struct NonOverlappingConsensus {
    f_record: Record,
    r_record: Record,
    likelihood: LogProb,
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
        let spinner_style = indicatif::ProgressStyle::default_spinner()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ")
            .template("{prefix:.bold.dim} {spinner} {wide_msg}");

        // cluster by umi
        // Note: If starcode is not installed, this throws a
        // hard to interpret error:
        // (No such file or directory (os error 2))
        // The expect added below should make this more clear.
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

        // prepare spinner for user feedback
        let pb = indicatif::ProgressBar::new_spinner();
        pb.set_style(spinner_style.clone());
        pb.set_prefix(&format!(
            "[1/2] Clustering input reads by UMI using starcode."
        ));

        loop {
            // update spinner
            pb.set_message(&format!("  Processed {:>10} reads", i));
            pb.inc(1);
            self.fq1_reader().read(&mut f_rec)?;
            self.fq2_reader().read(&mut r_rec)?;

            match (f_rec.is_empty(), r_rec.is_empty()) {
                (true, true) => break,
                (false, false) => (),
                (true, false) => {
                    let error_message = format!("Given FASTQ files have unequal lengths. Forward file returned record {} as empty, reverse record is not: id:'{}' seq:'{:?}'.", i, r_rec.id(), str::from_utf8(r_rec.seq()));
                    panic!(error_message);
                }
                (false, true) => {
                    let error_message = format!("Given FASTQ files have unequal lengths. Reverse file returned record {} as empty, forward record is not: id:'{}' seq:'{:?}'.", i, f_rec.id(), str::from_utf8(f_rec.seq()));
                    panic!(error_message);
                }
            }
            // extract umi for clustering
            let umi = if self.reverse_umi() {
                r_rec.seq()[..self.umi_len()].to_owned()
            } else {
                f_rec.seq()[..self.umi_len()].to_owned()
            };
            umi_cluster.stdin.as_mut().unwrap().write(&umi)?;
            umi_cluster.stdin.as_mut().unwrap().write(b"\n")?;
            // remove umi from read sequence for all further clustering steps
            if self.reverse_umi() {
                r_rec = self.strip_umi_from_record(&r_rec)
            } else {
                f_rec = self.strip_umi_from_record(&f_rec)
            }
            // store read sequences in an on-disk key value store for random access
            read_storage.put(i, &f_rec, &r_rec)?;
            i += 1;
        }
        umi_cluster.stdin.as_mut().unwrap().flush()?;
        drop(umi_cluster.stdin.take());
        pb.finish_with_message(&format!("Done. Analyzed {} reads.", i));

        // prepare user feedback
        let mut j = 0;
        let pb = indicatif::ProgressBar::new_spinner();
        pb.set_style(spinner_style.clone());
        pb.set_prefix(&format!(
            "[2/2] Merge eligable reads within each cluster.    "
        ));
        // read clusters identified by the first starcode run
        // the first run clustered by UMI, hence all reads in
        // the clusters handled here had similar UMIs
        for record in csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(umi_cluster.stdout.as_mut().unwrap())
            .records()
        {
            // update spinner
            pb.inc(1);
            pb.set_message(&format!("Processed {:>10} cluster", j));
            let seqids = parse_cluster(record?)?;
            // cluster within in this cluster by read sequence
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
                // get sequences from rocksdb (key value store)
                let (f_rec, r_rec) = read_storage.get(seqid - 1).unwrap();
                // perform clustering using the concatenated read sequences
                // without the UMIs (remove in the first clustering step)
                seq_cluster
                    .stdin
                    .as_mut()
                    .unwrap()
                    .write(&[&f_rec.seq()[..], &r_rec.seq()[..]].concat())?;
                seq_cluster.stdin.as_mut().unwrap().write(b"\n")?;
            }
            seq_cluster.stdin.as_mut().unwrap().flush()?;
            drop(seq_cluster.stdin.take());

            // handle each potential unique read, i.e. clusters with similar
            // UMI and similar sequence
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
                Some(s) => eprintln!("Starcode failed with error code {}", s),
                None => eprintln!("Starcode was terminated by signal"),
            }
            j += 1;
        }
        pb.finish_with_message(&format!("Done. Processed {} cluster.", j));
        Ok(())
    }
    fn strip_umi_from_record(&mut self, record: &Record) -> Record {
        let rec_seq = &record.seq()[self.umi_len()..];
        let rec_qual = &record.qual()[self.umi_len()..];
        Record::with_attrs(record.id(), record.desc(), rec_seq, rec_qual)
    }
    fn write_records(
        &mut self,
        f_recs: Vec<Record>,
        r_recs: Vec<Record>,
        outer_seqids: Vec<usize>,
    ) -> Result<(), Box<dyn Error>>;
    fn fq1_reader(&mut self) -> &mut fastq::Reader<R>;
    fn fq2_reader(&mut self) -> &mut fastq::Reader<R>;
    fn umi_len(&self) -> usize;
    fn seq_dist(&self) -> usize;
    fn umi_dist(&self) -> usize;
    fn reverse_umi(&self) -> bool;
}

/// Struct for calling non-overlapping consensus reads
/// Implements Trait CallConsensusReads
#[derive(new)]
pub struct CallNonOverlappingConsensusRead<'a, R: io::Read, W: io::Write> {
    fq1_reader: &'a mut fastq::Reader<R>,
    fq2_reader: &'a mut fastq::Reader<R>,
    fq1_writer: &'a mut fastq::Writer<W>,
    fq2_writer: &'a mut fastq::Writer<W>,
    umi_len: usize,
    seq_dist: usize,
    umi_dist: usize,
    reverse_umi: bool,
    verbose_read_names: bool,
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
            self.fq1_writer.write_record(
                &CalcNonOverlappingConsensus::new(
                    &f_recs,
                    &outer_seqids,
                    uuid,
                    self.verbose_read_names,
                )
                .calc_consensus()
                .0,
            )?;
            self.fq2_writer.write_record(
                &CalcNonOverlappingConsensus::new(
                    &r_recs,
                    &outer_seqids,
                    uuid,
                    self.verbose_read_names,
                )
                .calc_consensus()
                .0,
            )?;
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
#[derive(new)]
pub struct CallOverlappingConsensusRead<'a, R: io::Read, W: io::Write> {
    fq1_reader: &'a mut fastq::Reader<R>,
    fq2_reader: &'a mut fastq::Reader<R>,
    fq1_writer: &'a mut fastq::Writer<W>,
    fq2_writer: &'a mut fastq::Writer<W>,
    fq3_writer: &'a mut fastq::Writer<W>,
    umi_len: usize,
    seq_dist: usize,
    umi_dist: usize,
    insert_size: usize,
    std_dev: usize,
    reverse_umi: bool,
    verbose_read_names: bool,
}

impl<'a, R: io::Read, W: io::Write> CallOverlappingConsensusRead<'a, R, W> {
    fn isize_highest_probability(&mut self, f_seq_len: usize, r_seq_len: usize) -> f64 {
        if f_seq_len + f_seq_len < self.insert_size {
            return self.insert_size as f64;
        } else if f_seq_len + r_seq_len > self.insert_size + 2 * self.std_dev {
            return (self.insert_size + 2 * self.std_dev) as f64;
        } else {
            return (f_seq_len + r_seq_len) as f64;
        }
    }

    fn maximum_likelihood_overlapping_consensus(
        &mut self,
        f_recs: &Vec<Record>,
        r_recs: &Vec<Record>,
        outer_seqids: &Vec<usize>,
        uuid: &str,
    ) -> OverlappingConsensus {
        //Returns consensus record by filtering overlaps with lowest hamming distance.
        //For these overlaps(insert sizes) the consensus reads and their likelihoods are calculated.
        //The read with maximum likelihood will be returned.
        let insert_sizes = ((self.insert_size - 2 * self.std_dev)
            ..(self.insert_size + 2 * self.std_dev))
            .filter_map(|insert_size| {
                median_hamming_distance(insert_size, &f_recs, &r_recs)
                    .filter(|&median_distance| median_distance < HAMMING_THRESHOLD)
                    .map(|_| insert_size)
            });
        insert_sizes
            .map(|insert_size| {
                let overlap = (f_recs[0].seq().len() + r_recs[0].seq().len()) - insert_size;
                let (consensus_record, lh_isize) = CalcOverlappingConsensus::new(
                    &f_recs,
                    &r_recs,
                    overlap,
                    &outer_seqids,
                    &uuid,
                    self.verbose_read_names,
                )
                .calc_consensus();
                let likelihood = lh_isize
                    + isize_pmf(
                        insert_size as f64,
                        self.insert_size as f64,
                        self.std_dev as f64,
                    );
                OverlappingConsensus {
                    record: consensus_record,
                    likelihood,
                }
            })
            .max_by_key(|consensus| NotNaN::new(*consensus.likelihood).unwrap())
            .unwrap()
    }

    fn maximum_likelihood_nonoverlapping_consensus(
        &mut self,
        f_recs: &Vec<Record>,
        r_recs: &Vec<Record>,
        outer_seqids: &Vec<usize>,
        uuid: &str,
    ) -> NonOverlappingConsensus {
        //Calculate non-overlapping consensus records and shared lh
        let (f_consensus_rec, f_lh) =
            CalcNonOverlappingConsensus::new(&f_recs, &outer_seqids, uuid, self.verbose_read_names)
                .calc_consensus();
        let (r_consensus_rec, r_lh) =
            CalcNonOverlappingConsensus::new(&r_recs, &outer_seqids, uuid, self.verbose_read_names)
                .calc_consensus();
        let overall_lh_isize = f_lh + r_lh;
        //Determine insert size with highest probability for non-overlapping records based on expected insert size
        let likeliest_isize =
            self.isize_highest_probability(f_recs[0].seq().len(), r_recs[0].seq().len());
        let overall_lh = overall_lh_isize
            + isize_pmf(
                likeliest_isize,
                self.insert_size as f64,
                self.std_dev as f64,
            );
        NonOverlappingConsensus {
            f_record: f_consensus_rec,
            r_record: r_consensus_rec,
            likelihood: overall_lh,
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
    ) -> Result<(), Box<dyn Error>> {
        //TODO Add deterministic uuid considering read ids
        let uuid = &Uuid::new_v4().to_hyphenated().to_string();
        let ol_consensus =
            self.maximum_likelihood_overlapping_consensus(&f_recs, &r_recs, &outer_seqids, uuid);
        let non_ol_consensus =
            self.maximum_likelihood_nonoverlapping_consensus(&f_recs, &r_recs, &outer_seqids, uuid);
        match ol_consensus.likelihood > non_ol_consensus.likelihood {
            true => self.fq3_writer.write_record(&ol_consensus.record)?,
            false => {
                self.fq1_writer.write_record(&non_ol_consensus.f_record)?;
                self.fq2_writer.write_record(&non_ol_consensus.r_record)?;
            }
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
