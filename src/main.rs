extern crate bio;
#[macro_use]
extern crate clap;
use clap::App;
extern crate itertools;

use itertools::Itertools;

pub mod fastq;

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml).get_matches();

    if let Some(matches) = matches.subcommand_matches("fastq-split") {
        if let Some(chunks) = matches.values_of("chunks") {
            fastq::split::split(&chunks.collect_vec())
        }
    }
}
