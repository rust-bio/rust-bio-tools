extern crate bio;
#[macro_use]
extern crate clap;
use clap::App;
extern crate itertools;
#[macro_use]
extern crate log;
extern crate fern;
extern crate rustc_serialize;
extern crate csv;
extern crate rust_htslib;

use itertools::Itertools;

pub mod fastq;
pub mod bam;

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml).version(env!("CARGO_PKG_VERSION")).get_matches();

    fern::init_global_logger(fern::DispatchConfig {
        format: Box::new(|msg, _, _| msg.to_owned()),
        output: vec![fern::OutputConfig::stderr()],
        level: log::LogLevelFilter::Info
    }, log::LogLevelFilter::Trace).unwrap();

    if let Some(matches) = matches.subcommand_matches("fastq-split") {
        fastq::split::split(&matches.values_of("chunks").unwrap().collect_vec())
    }
    else if let Some(matches) = matches.subcommand_matches("bam-depth") {
        bam::depth::depth(&matches.value_of("bam-path").unwrap(), value_t!(matches, "max-read-length", u32).unwrap_or(1000))
    }
}
