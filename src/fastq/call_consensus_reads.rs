//! Tool remove PCR duplicates from UMI-tagged reads.
//!
//! This tool takes two FASTQ files (forward and reverse)
//! and returns two FASTQ files in which all PCR duplicates
//! have been merged into a consensus read.
//! Duplicates are identified by a Unique Molecular Identifier (UMI).
//!
//! ## Requirements:
//!
//! * starcode
//!
//! 
//! ## Usage:
//!
//! ```bash
//! $ rbt call-consensus-reads \
//!   <Path to FASTQ file with forward reads> \
//!   <Path to FASTQ file with reverse reads> \
//!   <Path for output forward FASTQ file> \
//!   <Path for output reverse FASTQ file> \
//!   -l <Length of UMI sequence> \
//!   -D <Maximum distance between sequences in a cluster> \  # See step 1 below
//!   -d <Maximum distance between UMIs in a cluster>   # See step 2 below
//! ```
//!
//! ## Assumptions:
//!
//!  - Reads are of equal length
//!  - UMI is the prefix of the reverse reads
//!
//! ## Workflow:
//!
//! The three main steps are:
//!
//! 1. Cluster all Reads by their sequence  
//! 
//!    1. Remove UMI sequence from reverse read (and save it for later use)
//!    2. Concatenate forward and reverse sequence
//!    3. Cluster by concatenated sequence using starcode
//!
//!        ```
//!        Forward Read: [================]
//!        Reverse Read: [(UMI)-----------]
//!        Sequence for clustering: [================-----------]
//!        ```
//!   After this, each cluster contains reads with similar sequences,
//!   but potentially different UMIs. Hence, each cluster contains candidates
//!   for merging.
//!
//! 2. For each cluster from step one: Cluster all reads within the cluster
//!    by their UMI (again, using starcode). Each of clusters generated in
//!    This step contain reads with similar read sequences as well as
//!    similar UMIs.
//!    Hence, each (new) cluster contains reads with similar sequence and
//!    UMI which can be merged into a consensus read as PCR duplicates.
//!
//! 3. For each cluster from step two: Compute a consensus sequence.
//!
//!    At each position in the read, all bases and quality values are used
//!    to compute the base with Maximum a-posteriori probability (MAP).
//!
//!      1. For one position, compute the likelihood for the four alleles
//!         A, C, G, and T, incorporating the number of bases as well as
//!         their quality values.
//!      2. Choose the allele with the largest likelihood for the consensus read.
//!      3. Compute the quality value of the consensus read from the maximum posterior
//!         probability used to select the allele.
//! 
//! 4. Write consensus reads to output file.
//!
//!
//!
// Since this is a binary crate, documentation needs to be compiled with this 'ancient incantation':
// https://github.com/rust-lang/cargo/issues/1865#issuecomment-394179125
use std::cmp;
use std::error::Error;
use std::fs;
use std::io;
use std::io::{BufReader, Write};
use std::mem;
use std::process::{Command, Stdio};
use std::str;
use tempfile::tempdir;

use bio::io::fastq;
use bio::io::fastq::FastqRead;
use bio::stats::probs::{LogProb, PHREDProb};
use csv;
use flate2::bufread::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use itertools::Itertools;
use ordered_float::NotNaN;
use rocksdb::DB;
use serde_json;
use uuid::Uuid;

const ALLELES: &'static [u8] = b"ACGT";

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
    storage_dir: String,
}

impl FASTQStorage {
    /// Create a new FASTQStorage using a Rocksdb database
    /// that maps read indices to read seqeunces.
    pub fn new() -> Result<Self, Box<dyn Error>> {
        let storage_dir = tempdir()?;
        Ok(FASTQStorage {
            db: DB::open_default(storage_dir.join("db"))?,
            storage_dir: tmp_path.clone(),
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

const PROB_CONFUSION: LogProb = LogProb(-1.0986122886681098); // (1 / 3).ln()

/// Compute a consensus sequence for a collection of FASTQ reads.
///
/// For each position, compute the likelihood of each allele and
/// choose the most likely one. Write the most likely allele i.e. base
/// as sequence into the consensus sequence. The quality value is the
/// likelihood for this allele, encoded in PHRED+33.
pub fn calc_consensus(recs: &[fastq::Record], seqids: &[usize], uuid: &str) -> fastq::Record {
    let seq_len = recs[0].seq().len();
    let mut consensus_seq = Vec::with_capacity(seq_len);
    let mut consensus_qual = Vec::with_capacity(seq_len);

    // assert that all reads have the same length here
    let identical_lengths = || {
        let reference_length = recs[0].seq().len();
        recs.iter()
            .map(|rec| rec.seq().len())
            .all(|len| len == reference_length)
    };
    assert_eq!(
        identical_lengths(),
        true,
        "Read length of FASTQ records {:?} differ. Cannot compute consensus sequence.",
        seqids
    );
    // Potential workflow for different read lengths
    // compute consensus of all reads with max len
    // compute offset of all shorter reads
    // pad shorter reads
    // drop first consensus, compute consensus of full length reads and padded reads
    // ignore padded bases for consensus computation

    for i in 0..seq_len {
        /// Compute the likelihood for the given allele and read position.
        /// The allele (A, C, G, or T) is an explicit parameter,
        /// the position i is captured by the closure.
        ///
        /// Likelihoods are managed in log space.
        /// A matching base is scored with (1 - PHRED score), a mismatch
        /// with PHRED score + confusion constant.
        let likelihood = |allele: &u8| {
            let mut lh = LogProb::ln_one(); // posterior: log(P(theta)) = 1
            for rec in recs {
                let q = LogProb::from(PHREDProb::from((rec.qual()[i] - 33) as f64));
                lh += if *allele == rec.seq()[i].to_ascii_uppercase() {
                    q.ln_one_minus_exp()
                } else {
                    q + PROB_CONFUSION
                };
            }
            lh
        };

        // Maximum a-posteriori estimate for the consensus base.
        // Find the allele (theta \in ACGT) with the highest likelihood
        // given the bases at this position, weighted with their quality values
        let likelihoods = ALLELES.iter().map(&likelihood).collect_vec();
        let max_posterior = likelihoods
            .iter()
            .enumerate()
            .max_by_key(|&(_, &lh)| NotNaN::new(*lh).unwrap()) // argmax of MAP
            .unwrap()
            .0;

        // 
        let marginal = LogProb::ln_sum_exp(&likelihoods);
        // new base: MAP
        consensus_seq.push(ALLELES[max_posterior]);
        // new qual: (1 - MAP)
        let qual = (likelihoods[max_posterior] - marginal).ln_one_minus_exp();

        // Assume the maximal quality, if the likelihood is infinite
        let truncated_quality: f64;
        if (*PHREDProb::from(qual)).is_infinite() {
            truncated_quality = 41.0;
        } else {
            truncated_quality = *PHREDProb::from(qual);
        }
        // Truncate quality values to PHRED+33 range
        consensus_qual.push(cmp::min(74, (truncated_quality + 33.0) as u64) as u8);
    }
    let name = format!(
        "{}_consensus-read-from:{}",
        uuid,
        seqids.iter().map(|i| format!("{}", i)).join(",")
    );
    fastq::Record::with_attrs(&name, None, &consensus_seq, &consensus_qual)
}


/// Build readers for the given input and output FASTQ files and pass them to
/// `call_consensus_reads`.
///
/// The type of the readers (writers) depends on the file ending.
/// If the input file names end with '.gz' a gzipped reader (writer) is used.
pub fn call_consensus_reads_from_paths(
    fq1: &str,
    fq2: &str,
    fq1_out: &str,
    fq2_out: &str,
    umi_len: usize,
    seq_dist: usize,
    umi_dist: usize,
) -> Result<(), Box<dyn Error>> {
    eprintln!("Reading input files:\n    {}\n    {}", fq1, fq2);
    eprintln!("Writing output to:\n    {}\n    {}", fq1_out, fq2_out);
    
    match (fq1.ends_with(".gz"), fq2.ends_with(".gz"), fq1_out.ends_with(".gz"), fq2_out.ends_with(".gz")) {
        (false, false, false, false) => call_consensus_reads(
            &mut fastq::Reader::from_file(fq1)?,
            &mut fastq::Reader::from_file(fq2)?,
            &mut fastq::Writer::to_file(fq1_out)?,
            &mut fastq::Writer::to_file(fq2_out)?,
            umi_len,
            seq_dist,
            umi_dist,
        ),
        (true, true, false, false) => call_consensus_reads(
            &mut fastq::Reader::new(fs::File::open(fq1).map(BufReader::new).map(GzDecoder::new).expect("Couldn't read fq1")),
            &mut fastq::Reader::new(fs::File::open(fq2).map(BufReader::new).map(GzDecoder::new).expect("Couldn't read fq2")),
            &mut fastq::Writer::to_file(fq1_out)?,
            &mut fastq::Writer::to_file(fq2_out)?,
            umi_len,
            seq_dist,
            umi_dist,
        ),
        (false, false, true, true) => call_consensus_reads(
            &mut fastq::Reader::from_file(fq1)?,
            &mut fastq::Reader::from_file(fq2)?,
            &mut fastq::Writer::new(GzEncoder::new(fs::File::create(fq1_out).expect("Couldn't open fq1_out"), Compression::default())),
            &mut fastq::Writer::new(GzEncoder::new(fs::File::create(fq2_out).expect("Couldn't open fq2_out"), Compression::default())),
            umi_len,
            seq_dist,
            umi_dist,
        ),
        (true, true, true, true) => call_consensus_reads(
            &mut fastq::Reader::new(fs::File::open(fq1).map(BufReader::new).map(GzDecoder::new).expect("Couldn't read fq1")),
            &mut fastq::Reader::new(fs::File::open(fq2).map(BufReader::new).map(GzDecoder::new).expect("Couldn't read fq2")),
            &mut fastq::Writer::new(GzEncoder::new(fs::File::create(fq1_out).expect("Couldn't open fq1_out"), Compression::default())),
            &mut fastq::Writer::new(GzEncoder::new(fs::File::create(fq2_out).expect("Couldn't open fq2_out"), Compression::default())),
            umi_len,
            seq_dist,
            umi_dist,
        ),
        _ => panic!("Invalid combination of files. Each pair of files (input and output) need to be both gzipped or both not zipped.")
    }
}

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
pub fn call_consensus_reads<R: io::Read, W: io::Write>(
    fq1_reader: &mut fastq::Reader<R>,
    fq2_reader: &mut fastq::Reader<R>,
    fq1_writer: &mut fastq::Writer<W>,
    fq2_writer: &mut fastq::Writer<W>,
    umi_len: usize,
    seq_dist: usize,
    umi_dist: usize,
) -> Result<(), Box<dyn Error>> {
    // cluster by sequence
    // Note: If starcode is not installed, this throws a
    // hard to interpret error:
    // (No such file or directory (os error 2))
    let mut seq_cluster = Command::new("starcode")
        .arg("--dist")
        .arg(format!("{}", seq_dist))
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
        fq1_reader.read(&mut f_rec)?;
        fq2_reader.read(&mut r_rec)?;
        match (f_rec.is_empty(), r_rec.is_empty()) {
            (true, true) => break,
            (false, false) => (),
            _ => panic!("Given FASTQ files have unequal lengths"),
        }
        // save umis for second (intra cluster) clustering
        let umi = r_rec.seq()[..umi_len].to_owned();
        umis.push(umi);

        read_storage.put(i, &f_rec, &r_rec)?;

        let seq = [f_rec.seq(), &r_rec.seq()[umi_len..]].concat();
        seq_cluster.stdin.as_mut().unwrap().write(&seq)?;
        seq_cluster.stdin.as_mut().unwrap().write(b"\n")?;
        i += 1;
    }
    seq_cluster.stdin.as_mut().unwrap().flush()?;
    drop(seq_cluster.stdin.take());

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
            .arg(format!("{}", umi_dist))
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

            if f_recs.len() > 1 {
                let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                fq1_writer.write_record(&calc_consensus(&f_recs, &outer_seqids, uuid))?;
                fq2_writer.write_record(&calc_consensus(&r_recs, &outer_seqids, uuid))?;
            } else {
                fq1_writer.write_record(&f_recs[0])?;
                fq2_writer.write_record(&r_recs[0])?;
            }
        }
    }

    Ok(())
}
