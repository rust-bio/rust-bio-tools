//! Filter reads matching names in a text file into a new FASTQ file.
//!
//! ## Usage:
//!
//! Extract the read with identifier `A` from `test.fastq` into a new file `filtered.fastq`
//! ```bash
//! $ cat ids.txt
//! A
//!
//! $ cat test.fastq
//! @A
//! ACTCTATCTA
//! +
//! !!!!!!!!!!
//! @B
//! CTCTATCTCTA
//! +
//! !!!!!!!!!!!
//!
//! $ rbt fastq-filter ids.txt < test.fastq > filtered.fastq
//!
//! $ cat filtered.fastq
//! @A
//! ACTCTATCTA
//! +
//! !!!!!!!!!!
//! ```
//!
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::iter::FromIterator;

pub fn filter(ids_path: &str) -> Result<(), Box<dyn Error>> {
    let mut reader = fastq::Reader::new(io::stdin());
    let mut writer = fastq::Writer::new(io::stdout());
    let f = File::open(ids_path)?;
    let f = BufReader::new(f);
    let ids =
        HashSet::<String>::from_iter(f.lines().filter_map(Result::ok).collect::<Vec<String>>());

    let mut record = fastq::Record::new();

    loop {
        reader.read(&mut record)?;
        if record.is_empty() {
            return Ok(());
        }
        if !ids.contains(record.id()) {
            writer.write_record(&record)?;
        }
    }
}
