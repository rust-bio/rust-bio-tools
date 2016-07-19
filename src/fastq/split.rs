
use std::io;
use itertools::Itertools;
use bio::io::fastq;

pub fn split(out_paths: &[&str]) {
    let mut reader = fastq::Reader::new(io::stdin());
    let mut writers = out_paths.iter().map(|path| fastq::Writer::to_file(path).ok().unwrap()).collect_vec();
    let mut record = fastq::Record::new();
    let mut i = 0;
    let mut j = 0;
    loop {
        // TODO better error formatting?
        reader.read(&mut record).unwrap();
        if record.is_empty() {
            return;
        }
        writers[i].write_record(&record).unwrap();
        i = (i + 1) % writers.len();
        j += 1;
        if j % 1000 == 0 {
            info!("{} records written.", j);
        }
    }
}
