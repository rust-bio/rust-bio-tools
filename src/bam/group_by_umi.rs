
use std::error::Error;
use std::collections::HashMap;
use std::str;

use rust_htslib::bam;
use rust_htslib::bam::Read;

use rand;
use rand::Rng;
use cogset::{Dbscan, BruteScan, Point};


/// Four-bit encode a given DNA sequence. Non-ACGT characters are randomly replaced with A, C, G or T.
///
/// Note:
///     Right now, IUPAC ambiguity codes are not interpreted.
fn encode_seq<R: Rng>(seq: &[u8], rng: &mut R) -> u32 {
    let mut enc = 0;
    for c in seq {
        enc = enc << 4;
        enc |= match c {
            b'A' | b'a' => 0b1,
            b'C' | b'c' => 0b10,
            b'G' | b'g' => 0b100,
            b'T' | b't' => 0b1000,
            _ => 1 << rng.gen_range(0, 3)
        };
    }

    enc
}


custom_derive! {
    #[derive(NewtypeDeref, Copy, Clone, Debug, PartialEq, Eq, NewtypeBitXor)]
    pub struct Fingerprint(u64);
}


impl Fingerprint {
    /// Compute hamming distance between two 4bit-encoded k-mers.
    ///
    /// Note:
    ///     We enforce at most 1 bit per base
    ///     -> hamming_dist of 0b_0001 (A) and 0b_0010 (C) yields count_ones(0011) = 2.
    pub fn hamming_dist(&self, other: &Fingerprint) -> u64 {
        (*self ^ *other).count_ones() as u64 / 2
    }
}


/// Implement the cogset Point trait for u32 integers, i.e. 4bit-encoded k-mers.
impl Point for Fingerprint {
    fn dist(&self, other: &Fingerprint) -> f64 {
        self.hamming_dist(other) as f64
    }
}


pub fn group_by_umi(in_bam: &str, max_hamming_dist: u64) -> Result<(), Box<Error>> {
    let get_reader = || bam::Reader::from_path(in_bam);
    let mut reader = get_reader()?;
    let header = bam::header::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_stdout(&header)?;
    let mut rng = rand::thread_rng();

    let mut rec = bam::Record::new();
    let mut i = 0;
    let mut umis_and_infix = Vec::new();

    loop {
        if let Err(e) = reader.read(&mut rec) {
            if e.is_eof() {
                break;
            } else {
                return Err(Box::new(e));
            }
        }

        // read UMI sequence from BAM file
        let umi = rec.aux(b"RX").expect(&format!("Expected tag RX in record {}", i)).string();
        let umi = encode_seq(&umi, &mut rng);

        // extract slice from the midpoint of the read for clustering
        // for this to work well, the spacers must have been removed beforehand
        let seq = rec.seq().as_bytes();
        assert!(seq.len() >= 8, "Unsupported read length (<8 bases) in record {}.", i);

        let mid_point = seq.len() / 2;
        let read_infix = encode_seq(&seq[mid_point-4..mid_point+4], &mut rng);

        // combine 32 bits from the umi with 32 bits () from the read infix
        //    UMI               Read seq
        // ##########=====================*=====================
        // ##########          mid-4 /....*..../ mid+4
        // ##########/....*..../  <- umi + infix
        //
        umis_and_infix.push(Fingerprint((umi as u64) << 32 | read_infix as u64));

        i += 1;
    }

    let scanner = BruteScan::new(&umis_and_infix);
    // the only reasonable value for minPTS is 2, because even 1 PCR duplicate shall be recognized
    let mut dbscan = Dbscan::new(scanner, max_hamming_dist as f64, 2);

    let mut clusters = HashMap::new();
    let mut cluster_sizes = HashMap::new();
    for (cluster_id, cluster) in dbscan.by_ref().enumerate() {
        let mut count = cluster_sizes.entry(cluster.len()).or_insert(0);
        *count += 1;
        for read_id in cluster {
            clusters.insert(read_id, cluster_id);
        }
    }

    info!("Number of clusters: {}", cluster_sizes.values().sum::<i32>());
    info!("Cluster sizes:");
    for (size, count) in cluster_sizes {
        info!("{}: {}", size, count);
    }

    for read_id in dbscan.noise_points() {
        // current size is equivalent to next id
        let cluster_id = clusters.len();
        clusters.insert(*read_id, cluster_id);
    }

    let mut reader = get_reader()?;
    let mut read_id = 0;
    loop {
        if let Err(e) = reader.read(&mut rec) {
            if e.is_eof() {
                break;
            } else {
                return Err(Box::new(e));
            }
        }
        let cluster_id = clusters.get(&read_id).unwrap();

        let mi = bam::record::Aux::Integer(*cluster_id as i32);
        rec.push_aux(b"MI", &mi);
        writer.write(&rec)?;

        read_id += 1;
    }

    Ok(())
}
