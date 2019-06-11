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
use std::io;

use csv;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use serde::Deserialize;
use snafu::ResultExt;

use crate::errors;

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
) -> errors::Result<()> {
    let mut bam_reader = bam::IndexedReader::from_path(&bam_path)
        .context(errors::BamIndexedReaderError { filepath: bam_path })?;
    let bam_header = bam_reader.header().clone();
    let mut pos_reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(io::stdin());
    let mut csv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(io::BufWriter::new(io::stdout()));

    for (i, record) in pos_reader.deserialize().enumerate() {
        let record: PosRecord = record.context(errors::CsvReadError)?;

        // jump to correct position
        let tid = bam_header.tid(record.chrom.as_bytes()).unwrap();
        let start = cmp::max(record.pos as i32 - max_read_length as i32 - 1, 0) as u32;
        bam_reader
            .fetch(tid, start, start + max_read_length * 2)
            .context(errors::BamFetchError { tid })?;

        // iterate over pileups
        let mut covered = false;
        for pileup in bam_reader.pileup() {
            let pileup = pileup.context(errors::BamPileupError { path: bam_path })?;
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

                r#try!(csv_writer
                    .serialize((&record.chrom, record.pos, depth))
                    .context(errors::CsvWriteError));
                break;
            } else if pileup.pos() > record.pos {
                break;
            }
        }
        if !covered {
            r#try!(csv_writer
                .serialize((&record.chrom, record.pos, 0))
                .context(errors::CsvWriteError));
        }

        if (i + 1) % 100 == 0 {
            info!("{} records written.", i + 1);
        }
    }
    Ok(())
}
