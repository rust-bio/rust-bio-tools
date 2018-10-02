#[macro_use]
extern crate log;
extern crate fern;
#[macro_use]
extern crate clap;
extern crate bio;
extern crate itertools;
extern crate rustc_serialize;
extern crate csv;
extern crate rust_htslib;
#[macro_use]
extern crate quick_error;
#[macro_use]
extern crate custom_derive;
#[macro_use]
extern crate newtype_derive;
extern crate rand;
extern crate cogset;

use std::process;

use clap::App;
use itertools::Itertools;

pub mod fastq;
pub mod bam;
pub mod bcf;

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml)
                      .version(env!("CARGO_PKG_VERSION"))
                      .get_matches();

    fern::Dispatch::new()
                   .format(|out, message, _| out.finish(format_args!("{}", message)))
                   .level(
                       if matches.is_present("verbose") {
                           log::LogLevelFilter::Debug
                       } else {
                           log::LogLevelFilter::Info
                       }
                   )
                   .chain(std::io::stderr())
                   .apply().unwrap();

    if let Some(matches) = matches.subcommand_matches("fastq-split") {
        if let Err(e) = fastq::split::split(
            &matches.values_of("chunks").unwrap().collect_vec()
        ) {
            error!("{}", e);
            process::exit(1);
        }
    }
    else if let Some(matches) = matches.subcommand_matches("fastq-filter") {
        if let Err(e) = fastq::filter::filter(
            &matches.value_of("ids").unwrap()
        ) {
            error!("{}", e);
            process::exit(1);
        }
    }

    else if let Some(matches) = matches.subcommand_matches("bam-depth") {
        if let Err(e) = bam::depth::depth(
            &matches.value_of("bam-path").unwrap(),
            value_t!(matches, "max-read-length", u32).unwrap_or(1000),
            value_t!(matches, "include-flags", u16).unwrap_or(0),
            value_t!(matches, "exclude-flags", u16).unwrap_or(4 | 256 | 512 | 1024),
            value_t!(matches, "min-mapq", u8).unwrap_or(0)
        ) {
            error!("{}", e);
            process::exit(1);
        }
    } else if let Some(matches) = matches.subcommand_matches("vcf-to-txt") {
        if let Err(e) = bcf::to_txt::to_txt(
            &matches.values_of("info").map(|values| values.collect_vec()).unwrap_or(vec![]),
            &matches.values_of("format").map(|values| values.collect_vec()).unwrap_or(vec![]),
            matches.is_present("genotypes")
        ) {
            error!("{}", e);
            process::exit(1);
        }
    } else if let Some(matches) = matches.subcommand_matches("vcf-match") {
        if let Err(e) = bcf::match_variants::match_variants(
            matches.value_of("vcf").unwrap(),
            value_t!(matches, "max-dist", u32).unwrap_or(20),
            value_t!(matches, "max-len-diff", u32).unwrap_or(10)
        ) {
            error!("{}", e);
            process::exit(1);
        }
    } else if let Some(_) = matches.subcommand_matches("vcf-baf") {
        if let Err(e) = bcf::baf::calculate_baf() {
            error!("{}", e);
            process::exit(1);
        }
    } else if let Some(matches) = matches.subcommand_matches("group-by-umi") {
        if let Err(e) = bam::group_by_umi::group_by_umi(
            matches.value_of("bam-path").unwrap(),
            value_t!(matches, "max-fingerprint-dist", u64).unwrap()
        ) {
            error!("{}", e);
            process::exit(1);
        }
    }
}
