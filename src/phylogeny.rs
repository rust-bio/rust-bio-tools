//! Compute a phylogeny from a distance matrix, or the distance between two phylogenies.
//!
//! ## Usage:
//!
//! ```
//! $ rbt phylogeny --method {UPGMA,NeighborJoining} phylip.dist [newick.tree]
//! $ rbt robinsonfoulds newick1.tree newick2.tree [outfile]
//! ```

use anyhow::Error;
use std::{fs::File, io::Write};

use crate::cli::PhylogenyMethod;

pub fn phylogeny(
    input: std::path::PathBuf,
    method: PhylogenyMethod,
    output: Option<std::path::PathBuf>,
) -> Result<(), Error> {
    let dm = bio::io::phylip::from_file(&input)?;
    let phylogeny = (match method {
        PhylogenyMethod::UPGMA => bio::phylogeny::upgma,
        PhylogenyMethod::NeighborJoining => bio::phylogeny::neighbor_joining,
    })(dm);
    let s = bio::io::newick::to_string(&phylogeny)?;
    Ok(match output {
        Some(path) => {
            let mut f = File::create(path)?;
            f.write_all(s.as_bytes())?;
            f.write(b"\n")?;
        }
        None => {
            println!("{}", s);
        }
    })
}

pub fn robinson_foulds_distance(
    newick_1: &std::path::PathBuf,
    newick_2: &std::path::PathBuf,
    output: Option<std::path::PathBuf>,
) -> Result<(), Error> {
    let dist = bio::phylogeny::robinson_foulds_distance(
        &bio::io::newick::from_file(newick_1)?,
        &bio::io::newick::from_file(newick_2)?,
    );
    Ok(match output {
        Some(path) => {
            let mut f = File::create(path)?;
            f.write_all(dist.to_string().as_bytes())?;
            f.write(b"\n")?;
        }
        None => {
            println!("{}", dist.to_string());
        }
    })
}
