//! Documentation for Rust Bio Tools
use clap::{load_yaml, value_t};
use log::LevelFilter;

use clap::App;
use fern;
use itertools::Itertools;
use std::error::Error;

pub mod bam;
pub mod bcf;
pub mod common;
pub mod fastq;

fn main() -> Result<(), Box<dyn Error>> {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml)
        .version(env!("CARGO_PKG_VERSION"))
        .get_matches();

    fern::Dispatch::new()
        .format(|out, message, _| out.finish(format_args!("{}", message)))
        .level(if matches.is_present("verbose") {
            LevelFilter::Debug
        } else {
            LevelFilter::Info
        })
        .chain(std::io::stderr())
        .apply()
        .unwrap();

    match matches.subcommand() {
        ("fastq-split", Some(matches)) => {
            fastq::split::split(&matches.values_of("chunks").unwrap().collect_vec())
        }
        ("fastq-filter", Some(matches)) => fastq::filter::filter(&matches.value_of("ids").unwrap()),
        ("bam-depth", Some(matches)) => bam::depth::depth(
            &matches.value_of("bam-path").unwrap(),
            value_t!(matches, "max-read-length", u32).unwrap_or(1000),
            value_t!(matches, "include-flags", u16).unwrap_or(0),
            value_t!(matches, "exclude-flags", u16).unwrap_or(4 | 256 | 512 | 1024),
            value_t!(matches, "min-mapq", u8).unwrap_or(0),
        ),
        ("vcf-to-txt", Some(matches)) => bcf::to_txt::to_txt(
            &matches
                .values_of("info")
                .map(|values| values.collect_vec())
                .unwrap_or(vec![]),
            &matches
                .values_of("format")
                .map(|values| values.collect_vec())
                .unwrap_or(vec![]),
            matches.is_present("genotypes"),
        ),
        ("vcf-match", Some(matches)) => bcf::match_variants::match_variants(
            matches.value_of("vcf").unwrap(),
            value_t!(matches, "max-dist", u32).unwrap_or(20),
            value_t!(matches, "max-len-diff", u32).unwrap_or(10),
        ),
        ("vcf-baf", Some(_)) => bcf::baf::calculate_baf(),
        ("vcf-fix-iupac-alleles", Some(_)) => bcf::fix_iupac_alleles::fix_iupac_alleles(),
        ("vcf-annotate-dgidb", Some(matches)) => bcf::annotate_dgidb::annotate_dgidb(
            &matches.value_of("vcf").unwrap(),
            matches.value_of("api-path").unwrap().to_string(),
            &matches.value_of("field").unwrap(),
            value_t!(matches, "genes-per-request", usize).unwrap(),
        ),
        ("call-consensus-reads", Some(matches)) => match matches.subcommand() {
            ("fastq", Some(matches)) => {
                fastq::call_consensus_reads::call_consensus_reads_from_paths(
                    matches.value_of("fq1").unwrap(),
                    matches.value_of("fq2").unwrap(),
                    matches.value_of("consensus-fq1").unwrap(),
                    matches.value_of("consensus-fq2").unwrap(),
                    matches.value_of("consensus-fq3"),
                    value_t!(matches, "umi-len", usize).unwrap(),
                    value_t!(matches, "max-seq-dist", usize).unwrap(),
                    value_t!(matches, "max-umi-dist", usize).unwrap(),
                    matches.is_present("umi-on-reverse"),
                    matches.is_present("verbose-read-names"),
                    if matches.is_present("insert-size") {
                        Some(value_t!(matches, "insert-size", usize).unwrap())
                    } else {
                        None
                    },
                    if matches.is_present("std-dev") {
                        Some(value_t!(matches, "std-dev", usize).unwrap())
                    } else {
                        None
                    },
                )
            }
            ("bam", Some(matches)) => bam::call_consensus_reads::call_consensus_reads_from_paths(
                matches.value_of("bam").unwrap(),
                matches.value_of("consensus-bam").unwrap(),
                value_t!(matches, "max-seq-dist", usize).unwrap(),
                matches.is_present("verbose-read-names"),
            ),
            _ => unreachable!(),
        },
        // This cannot be reached, since the matches step of
        // clap assures that a valid subcommand is provided
        _ => unreachable!(),
    }
}
