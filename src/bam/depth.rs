//! Compute the depth of coverage in a BAM file for a list of reference sequences and positions.
//!
//! ## Input:
//! A BAM file and a positions file.
//! The positions file contains the name of one reference sequence and one position per line (tab separated).
//! Example:
//! ```
//! 16	1
//! 17	1
//! 17	2
//! 17	38
//! 17	39
//! ```
//!
//! Positions are read from stdin, the BAM file is the first argument.
//!
//! ## Output:
//! Depth are written to stdout as tab-separated lines, similar to the positions input.
//! Example:
//! ```
//! 16	1	0
//! 17	1	5
//! 17	2	5
//! 17	38	14
//! 17	39	13
//! ```
//!
//! ## Usage:
//!
//! ```bash
//! $ rbt bam-depth tests/test.bam < tests/pos.txt > tests/depth.txt
//! ```
//! Where `pos.txt` is a positions file, as described above.
//!
//!
use log::info;
use std::cmp;
use std::error::Error;
use std::io;

use csv;
use serde::Deserialize;

use rust_htslib::bam;
use rust_htslib::bam::Read;

#[derive(Deserialize, Debug)]
struct PosRecord {
    chrom: String,
    pos: u32,
}

pub fn depth(
    bam_path: &str,
    max_read_length: u32,
    include_flags: u16,
    exclude_flags: u16,
    min_mapq: u8,
) -> Result<(), Box<dyn Error>> {
    let mut bam_reader = bam::IndexedReader::from_path(&bam_path)?;
    let bam_header = bam_reader.header().clone();
    let mut pos_reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(io::stdin());
    let mut csv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(io::BufWriter::new(io::stdout()));

    for (i, record) in pos_reader.deserialize().enumerate() {
        let record: PosRecord = record?;

        // jump to correct position
        let tid = bam_header.tid(record.chrom.as_bytes()).unwrap();
        let start = cmp::max(record.pos as i32 - max_read_length as i32 - 1, 0) as u32;
        bam_reader.fetch(tid, start, start + max_read_length * 2)?;

        // iterate over pileups
        let mut covered = false;
        for pileup in bam_reader.pileup() {
            let pileup = pileup?;
            covered = pileup.pos() == record.pos - 1;

            if covered {
                let depth = pileup
                    .alignments()
                    .filter(|alignment| {
                        let record = alignment.record();
                        let flags = record.flags();
                        (!flags) & include_flags == 0
                            && flags & exclude_flags == 0
                            && record.mapq() >= min_mapq
                    })
                    .count();

                csv_writer.serialize((&record.chrom, record.pos, depth))?;
                break;
            } else if pileup.pos() > record.pos {
                break;
            }
        }
        if !covered {
            csv_writer.serialize((&record.chrom, record.pos, 0))?;
        }

        if (i + 1) % 100 == 0 {
            info!("{} records written.", i + 1);
        }
    }
    Ok(())
}
