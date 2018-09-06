
use std::io;
use std::error::Error;
use bio::io::fastq;
use rust_htslib::bam;

use rand;
use rand::Rng;
use cogset::{Dbscan, BruteScan, Point};


const BIT: u32 = 1u32;


fn encode_seq(seq: &[u8], rng: &mut Rng) -> u32 {
    let mut enc = 0;
    for c in seq {
        enc = enc << 4;
        enc |= match c {
            b'A' | b'a' => 0b1,
            b'C' | b'c' => 0b10,
            b'G' | b'g' => 0b100,
            b'T' | b't' => 0b1000,
            _ => 1 << rng.get_range(0, 3)
        };
    }

    enc
}

fn hamming_dist(a: u32, b: u32) -> u32 {
    (a ^ b).count_ones()
}


impl Point for u32 {
    fn dist(other: u32) -> f64 {
        hamming_dist(self, other) as f64
    }
}


fn parse_dbr_pattern(dbr_pattern: &[u8]) -> (u32, u32) {
    let mut n = 1;
    let mut variable = 0;
    for c in dbr_pattern {
        n *= match c {
            b'N' => variable += 1; 4,
            b'M' | b'R' | b'W' | b'S' | b'Y' | b'K' => variable += 1; 2,
            b'V' | b'H' | b'D' | b'B' => variable += 1; 3,
            b'A' | b'C' | b'G' | b'T' => 1
        }
    }
    (n, variable)
}


pub fn group_by_dbr(in_bam: &str, dbr_pattern: &[u8]) -> Result<(), Box<Error>> {
    let (dbr_n, dbr_len) = parse_dbr_pattern(dbr_pattern);

    let mut reader = bam::Reader::from_path(in_bam)?;
    let mut writer = bam::Writer::from_stdout(reader.header());
    let mut rng = rand::thread_rng();

    let mut rec = bam::Record::new();
    let mut i = 0;
    let mut dbrs_and_infix = Vec::with_capacity(dbr_n);

    loop {
        if let Err(e) = reader.read(&mut rec) {
            if e.is_eof() {
                break;
            } else {
                return Err(e);
            }
        }

        // read DBR sequence from BAM file
        let dbr = rec.aux(b"RX").expect(&format!("Expected tag RX in record {}", i));
        let dbr = encode_seq(dbr[..dbr_len], &mut rng);

        // extract slice from the midpoint of the read for clustering
        // for this to work well, the spacers must have been removed beforehand
        let mid_point = (rec.seq_len / 2);
        let read_infix = encode_seq(rec.seq[mid_point-4..mid_point+4], &mut rng);

        // combine 32 bits from the dbr with 32 bits from the read infix
        dbrs_and_infix.push((dbr as u64) << 32 | read_infix as u64);

        
        i += 1;
    }

    let scanner = BruteScan::new(&dbrs_and_infix);
    // TODO what to set for epsilon?
    let mut dbscan = Dbscan::new(scanner, 1.0, 1);
    // TODO go on
}
