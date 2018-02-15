use std::io::{self, BufRead, BufReader};
use std::fs::File;
use std::error::Error;
use std::collections::HashSet;
use std::iter::FromIterator;
use bio::io::fastq;

pub fn filter(ids_path: &str) -> Result<(), Box<Error>> {
    let mut reader = fastq::Reader::new(io::stdin());
    let mut writer = fastq::Writer::new(io::stdout());
    let f = File::open(ids_path)?;
    let f = BufReader::new(f);
    let ids = HashSet::<String>::from_iter(
        f.lines().filter_map(Result::ok).collect::<Vec<String>>()
    );  

    let mut record = fastq::Record::new();

    loop {
        try!(reader.read(&mut record));
        if record.is_empty() {
            return Ok(());
        }
        if !ids.contains(record.id()) {
            try!(writer.write_record(&record));
        }
    }
}