//! Documentation for Rust Bio Tools
use clap::{load_yaml, value_t};
use log::LevelFilter;

use clap::App;
use itertools::Itertools;
use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::path::Path;
use std::str::FromStr;

pub mod bam;
pub mod bcf;
pub mod common;
pub mod fastq;
pub mod sequences_stats;

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
                .unwrap_or_default(),
            &matches
                .values_of("format")
                .map(|values| values.collect_vec())
                .unwrap_or_default(),
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
            match matches.is_present("datasources") {
                true => Some(matches.values_of("datasources").unwrap().collect()),
                false => None,
            },
            value_t!(matches, "genes-per-request", usize).unwrap(),
        ),
        ("vcf-report", Some(matches)) => {
            let mut sample_calls = HashMap::new();
            let mut bam_paths = HashMap::new();
            let output_path = matches.value_of("output-path").unwrap();
            let max_cells = u32::from_str(matches.value_of("max-cells").unwrap()).unwrap();
            bcf::report::embed_js(output_path)?;
            bcf::report::embed_css(output_path)?;
            bcf::report::embed_html(output_path)?;
            let fasta_path = matches.value_of("fasta").unwrap();
            let detail_path = output_path.to_owned() + "/details/";
            fs::create_dir(Path::new(&detail_path))?;
            let vcfs = matches.values_of("vcfs").unwrap();
            let bams = matches.values_of("bams").unwrap();
            for vcf in vcfs {
                let v: Vec<_> = vcf.split('=').collect();
                match sample_calls.insert(v[0].to_owned(), v[1].to_owned()) {
                    None => {}
                    _ => panic!("Found duplicate sample name {}. Please make sure the provided sample names are unique.", v[0].to_owned())
                }
            }
            for bam in bams {
                let b: Vec<_> = bam.split('=').collect();
                match bam_paths.insert(b[0].to_owned(), b[1].to_owned()) {
                    None => {}
                    _ => panic!("Found duplicate sample name {}. Please make sure the provided sample names are unique.", b[0].to_owned())
                }
            }
            for sample in sample_calls.keys().sorted() {
                bcf::report::table_report::table_report(
                    sample_calls.get(sample).unwrap(),
                    fasta_path,
                    bam_paths
                        .get(sample)
                        .unwrap_or_else(|| panic!("No bam provided for sample {}", sample)),
                    output_path,
                    sample,
                )?;
            }

            bcf::report::oncoprint::oncoprint(&sample_calls, output_path, max_cells)
        }
        ("collapse-reads-to-fragments", Some(matches)) => match matches.subcommand() {
            ("fastq", Some(matches)) => {
                fastq::collapse_reads_to_fragments::call_consensus_reads_from_paths(
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
            ("bam", Some(matches)) => {
                bam::collapse_reads_to_fragments::call_consensus_reads_from_paths(
                    matches.value_of("bam").unwrap(),
                    matches.value_of("consensus-bam").unwrap(),
                    value_t!(matches, "max-seq-dist", usize).unwrap(),
                    matches.is_present("verbose-read-names"),
                )
            }
            _ => unreachable!(),
        },
        ("sequence-stats", Some(matches)) => sequences_stats::stats(matches.is_present("fastq")),
        // This cannot be reached, since the matches step of
        // clap assures that a valid subcommand is provided
        _ => unreachable!(),
    }
}
