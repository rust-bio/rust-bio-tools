use bio::io::fasta;
use serde::Serialize;
use std::path::Path;

pub fn read_fasta(path: &Path, chrom: String, start: u64, stop: u64) -> Vec<Nucleobase> {
    let mut reader = fasta::IndexedReader::from_file(&path).unwrap();
    let index = fasta::Index::with_fasta_file(&path).unwrap();
    let _sequences = index.sequences();

    let mut seq: Vec<u8> = Vec::new();

    reader.fetch(&chrom, start, stop).unwrap();
    reader.read(&mut seq).unwrap();

    let mut fasta = Vec::new();
    let mut ind = start;
    for a in seq {
        let b = Nucleobase {
            start_position: ind as f64 - 0.5,
            end_position: ind as f64 + 0.5,
            marker_type: a as char,
            row: 0,
        };
        fasta.push(b);
        ind += 1;
    }

    fasta
}

pub fn get_fasta_length(path: &Path) -> u64 {
    let index = fasta::Index::with_fasta_file(&path).unwrap();
    let sequences = index.sequences();

    let length = sequences[0].len;

    length
}

#[derive(Serialize, Clone, Debug, PartialEq)]
pub struct Nucleobase {
    start_position: f64,
    end_position: f64,
    marker_type: char,
    row: u8,
}

impl Nucleobase {
    pub fn get_marker_type(&self) -> char {
        self.marker_type
    }
}
