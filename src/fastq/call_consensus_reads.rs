use std::cmp;
use std::error::Error;
use std::io;
use std::io::{BufReader, Write};
use std::mem;
use std::process::{Command, Stdio};
use std::str;
use std::fs;
use tempfile::tempdir;

use bio::io::fastq;
use bio::stats::probs::{LogProb, PHREDProb};
use itertools::Itertools;
use ordered_float::NotNaN;
use rocksdb::DB;
use serde_json;
use uuid::Uuid;
use csv;
use flate2::bufread::GzDecoder;


const ALLELES: &'static [u8] = b"ACGT";


/// Collects UMIs from a reader on a p7 FASTQ file and returns them in a vector.
///
/// This takes the first umi_len characters from each sequence in the file and
/// hence assumes that the UMI are the first characters in the line.
/// If the read other sequences before the UMI, these need to be trimmed before
/// this.
///
/// The UMI sequences are cloned into the vector.
///
/// # Errors
/// Passes on errors from bio::io::fastq::Reader i.e. fails if the file cannot
/// be opened.
///
/// # Examples
///
/// ```
/// let p7_fq = fastq::Reader::from_file(fq2);
/// let umi_len = 13;
/// let umis = umis(&mut p7_fq, umi_len)?;
/// ```
///
fn umis<R: io::Read>(
    fq_reader: &mut fastq::Reader<R>,
    umi_len: usize,
) -> Result<Vec<Vec<u8>>, Box<Error>> {
    let mut record = fastq::Record::new();
    let mut umis = Vec::new();
    loop {
        fq_reader.read(&mut record)?;
        if record.is_empty() {
            break;
        }

        let umi = record.seq()[..umi_len].to_owned();
        umis.push(umi);
    }

    Ok(umis)
}


fn parse_cluster(record: csv::StringRecord) -> Result<Vec<usize>, Box<Error>> {
    let seqids = &record[2];
    Ok(csv::ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(false)
        .from_reader(seqids.as_bytes())
        .deserialize()
        .next()
        .unwrap()?)
}

#[derive(Debug)]
pub struct FASTQStorage {
    db: DB,
}

impl FASTQStorage {
    pub fn new() -> Result<Self, Box<Error>> {
        let storage_dir = tempdir()?;
        Ok(FASTQStorage {
            db: DB::open_default(storage_dir.path().join("db"))?,
        })
    }

    fn as_key<'a>(i: u64) -> [u8; 8] {
        unsafe { mem::transmute::<u64, [u8; 8]>(i) }
    }

    pub fn put(
        &mut self,
        i: usize,
        f_rec: &fastq::Record,
        r_rec: &fastq::Record,
    ) -> Result<(), Box<Error>> {
        Ok(self.db.put(
            &Self::as_key(i as u64),
            serde_json::to_string(&(f_rec, r_rec))?.as_bytes(),
        )?)
    }

    pub fn get(&self, i: usize) -> Result<(fastq::Record, fastq::Record), Box<Error>> {
        Ok(serde_json::from_str(
            str::from_utf8(&self.db.get(&Self::as_key(i as u64))?.unwrap()).unwrap(),
        )?)
    }
}

const PROB_CONFUSION: LogProb = LogProb(-1.0986122886681098); // (1 / 3).ln()

pub fn calc_consensus(recs: &[fastq::Record], seqids: &[usize]) -> fastq::Record {
    let seq_len = recs[0].seq().len();
    let mut consensus_seq = Vec::with_capacity(seq_len);
    let mut consensus_qual = Vec::with_capacity(seq_len);

    for r in recs {
        eprintln!("{:?}", std::str::from_utf8(r.seq()));
    }
    for r in recs {
        eprintln!("{:?}", std::str::from_utf8(r.qual()));
    }
    
    for i in 0..seq_len {
        let likelihood = |allele: &u8| {
            let mut lh = LogProb::ln_one();
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
        
        let likelihoods = ALLELES.iter().map(&likelihood).collect_vec();
        eprintln!("{:?}", likelihoods);

        let max_posterior = likelihoods
            .iter()
            .enumerate()
            .max_by_key(|&(_, &lh)| NotNaN::new(*lh).unwrap())
            .unwrap()
            .0;

        let marginal = LogProb::ln_sum_exp(&likelihoods);
        // new base: MAP
        consensus_seq.push(ALLELES[max_posterior]);
        // new qual: (1 - MAP)
        let qual = (likelihoods[max_posterior] - marginal).ln_one_minus_exp();

        let truncated_quality: f64;
        if (*PHREDProb::from(qual)).is_infinite() {
            truncated_quality = 41.0;
        } else {
            truncated_quality = *PHREDProb::from(qual);
        }
        eprintln!("LL(MAP) - marginal = {:?} - {:?} = {:?}", likelihoods[max_posterior], marginal, likelihoods[max_posterior] - marginal);
        eprintln!("=> Quality: {:?} => {:?}", qual, PHREDProb::from(qual));
        eprintln!("as u64: {:?}", (*PHREDProb::from(qual) + 33.0) as u64);
        eprintln!("as u8: {:?}\n", ((*PHREDProb::from(qual) + 33.0) as u64) as u8);
        eprintln!("truncated_quality {:?}", truncated_quality);
        consensus_qual.push(cmp::min(74, (truncated_quality + 33.0) as u64) as u8);

        // This is the old code:
        // consensus_qual.push(cmp::min(255, (*PHREDProb::from(qual) + 33.0) as u64) as u8);
    }
    eprintln!("{:?}", consensus_seq);
    eprintln!("{:?}", consensus_qual);

    fastq::Record::with_attrs(
        &Uuid::new_v4().to_hyphenated().to_string(),
        Some(&format!(
            "consensus-read-from: {}",
            seqids.iter().map(|i| format!("{}", i)).join(",")
        )),
        &consensus_seq,
        &consensus_qual,
    )
}

pub fn call_consensus_reads(
    fq1: &str,
    fq2: &str,
    fq1_out: &str,
    fq2_out: &str,
    umi_len: usize,
    seq_dist: usize,
    umi_dist: usize,
) -> Result<(), Box<Error>> {
    // TODO Opening a gz file should be optional.
    // Right now, this fails for non-gzipped files
    // Below, opening the p5 file has the same issue.
    let load_fq2 = || fastq::Reader::new(fs::File::open(fq2)
                                         .map(BufReader::new)
                                         .map(GzDecoder::new).unwrap());
    let umis = umis(&mut load_fq2(), umi_len)?;
    // cluster by sequence
    // if starcode is not installed, this throws a hard to interpret error:
    // (No such file or directory (os error 2))
    let mut seq_cluster = Command::new("starcode")
        .arg("--dist")
        .arg(format!("{}", seq_dist))
        .arg("--seq-id")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;
    let mut fq1_reader = fastq::Reader::new(fs::File::open(fq1)
                                            .map(BufReader::new)
                                            .map(GzDecoder::new).unwrap());
    let mut fq2_reader = load_fq2();
    let mut f_rec = fastq::Record::new();
    let mut r_rec = fastq::Record::new();

    // init temp storage for reads
    let mut read_storage = FASTQStorage::new()?;
    let mut i = 0;
    loop {
        fq1_reader.read(&mut f_rec)?;
        fq2_reader.read(&mut r_rec)?;
        match (f_rec.is_empty(), r_rec.is_empty()) {
            (true, true) => break,
            (false, false) => (),
            _ => panic!("Given FASTQ files have unequal lengths"),
        }

        read_storage.put(i, &f_rec, &r_rec)?;

        let seq = [f_rec.seq(), &r_rec.seq()[umi_len..]].concat();
        seq_cluster.stdin.as_mut().unwrap().write(&seq)?;
        seq_cluster.stdin.as_mut().unwrap().write(b"\n")?;
        i += 1;
    }
    eprintln!("Read Storage {:?}", read_storage);
    seq_cluster.stdin.as_mut().unwrap().flush()?;
    drop(seq_cluster.stdin.take());

    let mut fq1_writer = fastq::Writer::to_file(fq1_out)?;
    let mut fq2_writer = fastq::Writer::to_file(fq2_out)?;
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
        eprintln!("this is reached");
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
            eprintln!("{:?}", inner_seqids);
            for inner_seqid in inner_seqids {
                eprintln!("{:?} {:?}", seqids, inner_seqid);
                let seqid = seqids[inner_seqid - 1];
                let (f_rec, r_rec) = read_storage.get(seqid - 1)?;
                f_recs.push(f_rec);
                r_recs.push(r_rec);
                outer_seqids.push(seqid);
            }
            eprintln!("Even this is reached: {}", f_recs.len());
            // This is where the ERROR happends:
            if f_recs.len() > 1 {
                fq1_writer.write_record(&calc_consensus(&f_recs, &outer_seqids))?;
                fq2_writer.write_record(&calc_consensus(&r_recs, &outer_seqids))?;
            } else {
                fq1_writer.write_record(&f_recs[0])?;
                fq2_writer.write_record(&r_recs[0])?;
            }
        }
    }

    Ok(())
}
