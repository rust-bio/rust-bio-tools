
use std::error::Error;
use std::collections::HashMap;

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
    ///     Right now this is overestimating the distance. We enforce at most 1 bit per base
    ///     -> hamming_dist of 0b_0001 (A) and 0b_0010 (C) yields count_ones(0011) = 2.
    ///     This should not inhibit the distance computation in any way, since everything is
    ///     scaled by the factor 2, but might leed to unintuitive results.
    pub fn hamming_dist(&self, other: &Fingerprint) -> u64 {
        (*self ^ *other).count_ones() as u64
    }
}


/// Implement the cogset Point trait for u32 integers, i.e. 4bit-encoded k-mers.
impl Point for Fingerprint {
    fn dist(&self, other: &Fingerprint) -> f64 {
        self.hamming_dist(other) as f64
    }
}


/// Read a dbr pattern encoded as IUPAC ambiguity codes, compute the number
/// of possible configurations and the number of variable positions.
fn parse_dbr_pattern(dbr_pattern: &[u8]) -> (usize, usize) {
    let mut n = 1;
    let mut variable = 0;
    for c in dbr_pattern {
        n *= match c {
            b'N' => { variable += 1; 4 },
            b'M' | b'R' | b'W' | b'S' | b'Y' | b'K' => { variable += 1; 2},
            b'V' | b'H' | b'D' | b'B' => { variable += 1; 3 },
            b'A' | b'C' | b'G' | b'T' => 1,
            _ => panic!("unexpected character in sequence")
        }
    }
    (n, variable)
}


pub fn group_by_dbr(in_bam: &str, dbr_pattern: &[u8], max_hamming_dist: u64, min_cluster_size: usize) -> Result<(), Box<Error>> {
    let (dbr_n, dbr_len) = parse_dbr_pattern(dbr_pattern);

    let get_reader = || bam::Reader::from_path(in_bam);
    let mut reader = get_reader()?;
    let header = bam::header::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_stdout(&header)?;
    let mut rng = rand::thread_rng();

    let mut rec = bam::Record::new();
    let mut i = 0;
    let mut dbrs_and_infix = Vec::with_capacity(dbr_n);

    loop {
        if let Err(e) = reader.read(&mut rec) {
            if e.is_eof() {
                break;
            } else {
                return Err(Box::new(e));
            }
        }

        // read DBR sequence from BAM file
        let dbr = rec.aux(b"RX").expect(&format!("Expected tag RX in record {}", i)).string();
        let dbr = encode_seq(&dbr[..dbr_len], &mut rng);

        // extract slice from the midpoint of the read for clustering
        // for this to work well, the spacers must have been removed beforehand
        let seq = rec.seq().as_bytes();
        let mid_point = seq.len() / 2;
        let read_infix = encode_seq(&seq[mid_point-4..mid_point+4], &mut rng);

        // combine 32 bits from the dbr with 32 bits () from the read infix
        //    DBR               Read seq
        // ##########=====================*=====================
        // ##########          mid-4 /....*..../ mid+4
        // ##########/....*..../  <- dbr + infix
        //
        dbrs_and_infix.push(Fingerprint((dbr as u64) << 32 | read_infix as u64));


        i += 1;
    }

    let scanner = BruteScan::new(&dbrs_and_infix);
    let mut dbscan = Dbscan::new(scanner, max_hamming_dist as f64, min_cluster_size);

    let mut clusters = HashMap::new();
    for (cluster_id, cluster) in dbscan.by_ref().enumerate() {
        for read_id in cluster {
            clusters.insert(read_id, cluster_id);
        }
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
