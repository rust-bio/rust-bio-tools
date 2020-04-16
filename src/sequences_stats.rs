//! Compute statics on sequences from stdin:
//!   - min: length of shortest sequence
//!   - max: length of longest sequence
//!   - average: average length of sequence
//!   - median: median length of sequence
//!   - nb_reads: number of reads
//!   - nb_bases: number of bases
//!   - n50: N50 of sequences
//!
//! ## Usage:
//!
//! ```
//! $ rbt sequences-stats < A.fasta
//! $ rbt sequences-stats -q < A.fastq
//! ```

use bio::io::{fasta, fastq};

use quick_error::quick_error;

use std::error::Error;
use std::io;

quick_error! {
    #[derive(Debug)]
    pub enum InputError {
        NoSequence {
            description("stdin didn't contains any sequence")
        }

    }
}

pub fn stats(fastq: bool) -> Result<(), Box<dyn Error>> {
    let mut lengths = if fastq {
        fastq_lengths()
    } else {
        fasta_lengths()
    };

    if lengths.is_empty() {
        return Err(Box::new(InputError::NoSequence));
    }
    // Sort lengths one time
    lengths.sort();

    let nb_bases = lengths.iter().sum::<usize>();

    println!(
        "min: {min}
max: {max}
average: {average}
mediane: {mediane}
number of reads: {nb_reads}
number of bases: {nb_bases}
n50: {n50}",
        min = lengths[0],                 // First element is the minimal element
        max = lengths[lengths.len() - 1], // last element is the maximal element
        average = average(&lengths),
        mediane = median(&lengths),
        nb_reads = lengths.len(),
        nb_bases = nb_bases,
        n50 = n50(&lengths, nb_bases),
    );

    Ok(())
}

fn fasta_lengths() -> Vec<usize> {
    let reader = fasta::Reader::new(io::stdin());

    let mut lengths = Vec::new();

    let mut records = reader.records();
    while let Some(Ok(record)) = records.next() {
        lengths.push(record.seq().len());
    }

    lengths
}

pub fn fastq_lengths() -> Vec<usize> {
    let reader = fastq::Reader::new(io::stdin());

    let mut lengths = Vec::new();

    let mut records = reader.records();
    while let Some(Ok(record)) = records.next() {
        lengths.push(record.seq().len());
    }

    lengths
}

fn n50(numbers: &[usize], nb_bases_total: usize) -> usize {
    let mut acc = 0;
    for val in numbers.iter() {
        acc += *val;
        if acc > nb_bases_total / 2 {
            return *val;
        }
    }

    return numbers[numbers.len() - 1];
}

fn average(numbers: &[usize]) -> f64 {
    numbers.iter().sum::<usize>() as f64 / numbers.len() as f64
}

fn median(data: &[usize]) -> f64 {
    match data.len() {
        0 => 0.0,
        1 => data[0] as f64,
        len if len % 2 == 0 => {
            let v1 = data[(len / 2) - 1];
            let v2 = data[len / 2];
            (v1 + v2) as f64 / 2.0
        }
        len => data[len / 2] as f64,
    }
}
