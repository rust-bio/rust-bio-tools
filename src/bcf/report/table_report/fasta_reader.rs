use crate::common::Region;
use anyhow::Context;
use anyhow::Result;
use bio::io::fasta;
use itertools::Itertools;
use serde::Serialize;
use std::collections::HashMap;
use std::path::Path;

pub fn read_fasta<P: AsRef<Path>>(
    path: P,
    region: &Region,
    compensate_0_basing: bool,
) -> Result<Vec<Nucleobase>> {
    let mut reader = fasta::IndexedReader::from_file(&path).unwrap();
    let index =
        fasta::Index::with_fasta_file(&path).context("error reading index file of input FASTA")?;
    let _sequences = index.sequences();

    let mut seq: Vec<u8> = Vec::new();

    reader.fetch(&region.target, region.start, region.end)?;
    reader.read(&mut seq)?;

    let mut fasta = Vec::new();
    let mut ind = region.start;
    if compensate_0_basing {
        ind += 1;
    }
    for a in seq {
        let base = char::from(a);
        let marker = base.to_uppercase().collect_vec().pop().unwrap();
        let b = Nucleobase {
            position: ind,
            marker_type: marker,
            row: 0,
            repeat: base.is_lowercase(),
        };
        fasta.push(b);
        ind += 1;
    }

    Ok(fasta)
}

pub fn get_fasta_lengths(path: &Path) -> Result<HashMap<String, u64>> {
    let index = fasta::Index::with_fasta_file(&path).context("error reading input FASTA")?;
    let sequences = index.sequences();
    Ok(sequences
        .iter()
        .map(|s| (s.name.to_owned(), s.len))
        .collect())
}

#[derive(Serialize, Clone, Debug, PartialEq)]
pub struct Nucleobase {
    position: u64,
    marker_type: char,
    row: u8,
    repeat: bool,
}

impl Nucleobase {
    pub fn get_marker_type(&self) -> char {
        self.marker_type
    }
}
