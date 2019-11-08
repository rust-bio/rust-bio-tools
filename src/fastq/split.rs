//! Split reads from stdin up into the given files.
//!
//! ## Usage:
//!
//! Distribute reads from `test.fastq` into the files `A.fastq` and `B.fastq`.
//! ```bash
//! $ rbt fastq-split A.fastq B.fastq < test.fastq
//! ```
//!
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use log::info;
use std::error::Error;
use std::io;

pub fn split(out_paths: &[&str]) -> Result<(), Box<dyn Error>> {
    let mut reader = fastq::Reader::new(io::stdin());
    let mut writers = Vec::new();
    for path in out_paths {
        writers.push(fastq::Writer::to_file(path)?);
    }
    let mut record = fastq::Record::new();
    let mut i = 0;
    let mut j = 0;
    loop {
        reader.read(&mut record)?;
        if record.is_empty() {
            return Ok(());
        }
        writers[i].write_record(&record)?;
        i = (i + 1) % writers.len();
        j += 1;
        if j % 1000 == 0 {
            info!("{} records written.", j);
        }
    }
}
