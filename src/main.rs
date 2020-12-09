//! Documentation for Rust Bio Tools
use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::path::Path;
use std::str::FromStr;

use itertools::Itertools;
use log::LevelFilter;
use rayon::prelude::*;
use structopt::StructOpt;

use cli::Command::*;

pub mod bam;
pub mod bcf;
mod cli;
pub mod common;
pub mod fastq;
pub mod sequences_stats;

fn main() -> Result<(), Box<dyn Error>> {
    let args = cli::Rbt::from_args();

    fern::Dispatch::new()
        .format(|out, message, _| out.finish(format_args!("{}", message)))
        .level(if args.verbose {
            LevelFilter::Debug
        } else {
            LevelFilter::Info
        })
        .chain(std::io::stderr())
        .apply()
        .unwrap();

    match args.cmd {
        FastqSplit { chunks } => {
            fastq::split::split(&chunks.iter().map(|p| p.to_str().unwrap()).collect_vec())?
        }
        FastqFilter { ids } => fastq::filter::filter(&ids).unwrap(),
        BamDepth {
            bam_path,
            max_read_length,
            include_flags,
            exclude_flags,
            min_mapq,
        } => bam::depth::depth(
            &bam_path,
            max_read_length,
            include_flags,
            exclude_flags,
            min_mapq,
        )?,
        VcfToTxt {
            info,
            format,
            genotypes,
        } => bcf::to_txt::to_txt(
            &info.iter().map(|s| s as &str).collect_vec().as_slice(),
            &format.iter().map(|s| s as &str).collect_vec().as_slice(),
            genotypes,
        )?,
        VcfMatch {
            vcf,
            max_dist,
            max_len_diff,
        } => bcf::match_variants::match_variants(vcf, max_dist, max_len_diff)?,
        VcfBaf {} => bcf::baf::calculate_baf()?,
        VcfFixIupacAlleles {} => bcf::fix_iupac_alleles::fix_iupac_alleles()?,
        VcfAnnotateDgidb {
            vcf,
            api_path,
            field,
            datasources,
            genes_per_request,
        } => bcf::annotate_dgidb::annotate_dgidb(
            vcf,
            api_path,
            &*field,
            datasources.as_deref(),
            genes_per_request,
        )?,
        VcfReport {
            fasta,
            vcfs,
            bams,
            cells,
            max_read_depth,
            infos,
            formats,
            custom_js_template,
            custom_js_files,
            tsv,
            threads,
            output_path,
        } => {
            let mut sample_calls = HashMap::new();
            let mut bam_paths = HashMap::new();
            if !Path::new(&output_path).exists() {
                fs::create_dir(Path::new(&output_path)).unwrap_or_else(|_| {
                    panic!(
                        "Couldn't create output directory at {}. Please make sure the path exists.",
                        output_path
                    )
                });
            }
            let js_files_vec = custom_js_files
                .clone()
                .map_or_else(Vec::new, |values| values.into_iter().collect());
            let js_file_names = if let Some(files) = custom_js_files.clone() {
                files
                    .iter()
                    .map(|f| {
                        f.split('/')
                            .collect_vec()
                            .pop()
                            .unwrap_or_else(|| {
                                panic!("Unable to extract file name from path: {:?}", f)
                            })
                            .to_owned()
                    })
                    .collect()
            } else {
                vec![]
            };
            bcf::report::embed_js(&output_path, custom_js_template.as_deref(), js_files_vec)?;
            bcf::report::embed_css(&output_path)?;
            bcf::report::embed_html(&output_path)?;
            let detail_path = output_path.to_owned() + "/details/";
            fs::create_dir(Path::new(&detail_path))?;
            for vcf in vcfs {
                let v: Vec<_> = vcf.split('=').collect();
                match sample_calls.insert(v[0].to_owned(), v[1].to_owned()) {
                    None => {}
                    _ => panic!("Found duplicate sample name {}. Please make sure the provided sample names are unique.", v[0].to_owned())
                }
            }
            for bam in bams {
                let b: Vec<_> = bam.split('=').collect();
                let c: Vec<_> = b[0].split(':').collect();
                let rec = bam_paths.entry(c[0].to_owned()).or_insert_with(Vec::new);
                rec.push((c[1].to_owned(), b[1].to_owned()))
            }

            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()?;

            sample_calls.par_iter().for_each(|(sample, sample_call)| {
                bcf::report::table_report::table_report(
                    sample_call,
                    &fasta,
                    bam_paths
                        .get(sample)
                        .unwrap_or_else(|| panic!("No bam provided for sample {}", sample)),
                    &output_path,
                    sample,
                    infos.clone(),
                    formats.clone(),
                    max_read_depth,
                    js_file_names.clone(),
                )
                .unwrap_or_else(|_| panic!("Failed building table report for sample {}", sample));
            });

            bcf::report::oncoprint::oncoprint(&sample_calls, &output_path, cells, tsv.as_deref())?
        }
        _ => unimplemented!(),
    }
    Ok(())

    // match matches.subcommand() {
    //     ("fastq-split", Some(matches)) => {
    //         fastq::split::split(&matches.values_of("chunks").unwrap().collect_vec())
    //     }
    //     ("fastq-filter", Some(matches)) => fastq::filter::filter(&matches.value_of("ids").unwrap()),
    //     ("bam-depth", Some(matches)) => bam::depth::depth(
    //         &matches.value_of("bam-path").unwrap(),
    //         value_t!(matches, "max-read-length", u32).unwrap_or(1000),
    //         value_t!(matches, "include-flags", u16).unwrap_or(0),
    //         value_t!(matches, "exclude-flags", u16).unwrap_or(4 | 256 | 512 | 1024),
    //         value_t!(matches, "min-mapq", u8).unwrap_or(0),
    //     ),
    //     ("vcf-to-txt", Some(matches)) => bcf::to_txt::to_txt(
    //         &matches
    //             .values_of("info")
    //             .map(|values| values.collect_vec())
    //             .unwrap_or_default(),
    //         &matches
    //             .values_of("format")
    //             .map(|values| values.collect_vec())
    //             .unwrap_or_default(),
    //         matches.is_present("genotypes"),
    //     ),
    //     ("vcf-match", Some(matches)) => bcf::match_variants::match_variants(
    //         matches.value_of("vcf").unwrap(),
    //         value_t!(matches, "max-dist", u32).unwrap_or(20),
    //         value_t!(matches, "max-len-diff", u32).unwrap_or(10),
    //     ),
    //     ("vcf-baf", Some(_)) => bcf::baf::calculate_baf(),
    //     ("vcf-fix-iupac-alleles", Some(_)) => bcf::fix_iupac_alleles::fix_iupac_alleles(),
    //     ("vcf-annotate-dgidb", Some(matches)) => bcf::annotate_dgidb::annotate_dgidb(
    //         &matches.value_of("vcf").unwrap(),
    //         matches.value_of("api-path").unwrap().to_string(),
    //         &matches.value_of("field").unwrap(),
    //         match matches.is_present("datasources") {
    //             true => Some(matches.values_of("datasources").unwrap().collect()),
    //             false => None,
    //         },
    //         value_t!(matches, "genes-per-request", usize).unwrap(),
    //     ),
    //     ("vcf-report", Some(matches)) => {
    //         let mut sample_calls = HashMap::new();
    //         let mut bam_paths = HashMap::new();
    //         let output_path = matches.value_of("output-path").unwrap();
    //         let threads = usize::from_str(matches.value_of("threads").unwrap()).unwrap();
    //         if !Path::new(output_path).exists() {
    //             fs::create_dir(Path::new(output_path)).unwrap_or_else(|_| {
    //                 panic!(
    //                     "Couldn't create output directory at {}. Please make sure the path exists.",
    //                     output_path
    //                 )
    //             });
    //         }
    //         let max_cells = u32::from_str(matches.value_of("max-cells").unwrap()).unwrap();
    //         let max_read_depth =
    //             u32::from_str(matches.value_of("max-read-depth").unwrap()).unwrap();
    //         let custom_js = matches.value_of("custom-js-template");
    //         let js_files = matches.values_of("custom-js-files");
    //         let js_files_vec = js_files
    //             .clone()
    //             .map_or_else(Vec::new, |values| values.collect());
    //         let js_file_names = if let Some(files) = js_files.clone() {
    //             files
    //                 .map(|f| {
    //                     f.split('/')
    //                         .collect_vec()
    //                         .pop()
    //                         .unwrap_or_else(|| {
    //                             panic!("Unable to extract file name from path: {:?}", f)
    //                         })
    //                         .to_owned()
    //                 })
    //                 .collect()
    //         } else {
    //             vec![]
    //         };
    //         let tsv_data = matches.value_of("tsv");
    //         bcf::report::embed_js(output_path, custom_js, js_files_vec)?;
    //         bcf::report::embed_css(output_path)?;
    //         bcf::report::embed_html(output_path)?;
    //         let fasta_path = matches.value_of("fasta").unwrap();
    //         let detail_path = output_path.to_owned() + "/details/";
    //         fs::create_dir(Path::new(&detail_path))?;
    //         let vcfs = matches.values_of("vcfs").unwrap();
    //         let bams = matches.values_of("bams").unwrap();
    //         let infos = matches.values_of("info");
    //         let formats = matches.values_of("format");
    //         for vcf in vcfs {
    //             let v: Vec<_> = vcf.split('=').collect();
    //             match sample_calls.insert(v[0].to_owned(), v[1].to_owned()) {
    //                 None => {}
    //                 _ => panic!("Found duplicate sample name {}. Please make sure the provided sample names are unique.", v[0].to_owned())
    //             }
    //         }
    //         for bam in bams {
    //             let b: Vec<_> = bam.split('=').collect();
    //             let c: Vec<_> = b[0].split(':').collect();
    //             let rec = bam_paths.entry(c[0].to_owned()).or_insert_with(Vec::new);
    //             rec.push((c[1].to_owned(), b[1].to_owned()))
    //         }
    //
    //         rayon::ThreadPoolBuilder::new()
    //             .num_threads(threads)
    //             .build_global()?;
    //
    //         sample_calls.par_iter().for_each(|(sample, sample_call)| {
    //             bcf::report::table_report::table_report(
    //                 sample_call,
    //                 fasta_path,
    //                 bam_paths
    //                     .get(sample)
    //                     .unwrap_or_else(|| panic!("No bam provided for sample {}", sample)),
    //                 output_path,
    //                 sample,
    //                 infos.clone(),
    //                 formats.clone(),
    //                 max_read_depth,
    //                 js_file_names.clone(),
    //             )
    //             .unwrap_or_else(|_| panic!("Failed building table report for sample {}", sample));
    //         });
    //
    //         bcf::report::oncoprint::oncoprint(&sample_calls, output_path, max_cells, tsv_data)
    //     }
    //     ("vcf-split", Some(matches)) => {
    //         bcf::split::split(
    //             matches.value_of("input").unwrap(),
    //             &matches.values_of("output").unwrap().collect::<Vec<_>>(),
    //         )?;
    //         Ok(())
    //     }
    //     ("collapse-reads-to-fragments", Some(matches)) => match matches.subcommand() {
    //         ("fastq", Some(matches)) => {
    //             fastq::collapse_reads_to_fragments::call_consensus_reads_from_paths(
    //                 matches.value_of("fq1").unwrap(),
    //                 matches.value_of("fq2").unwrap(),
    //                 matches.value_of("consensus-fq1").unwrap(),
    //                 matches.value_of("consensus-fq2").unwrap(),
    //                 matches.value_of("consensus-fq3"),
    //                 value_t!(matches, "umi-len", usize).unwrap(),
    //                 value_t!(matches, "max-seq-dist", usize).unwrap(),
    //                 value_t!(matches, "max-umi-dist", usize).unwrap(),
    //                 matches.is_present("umi-on-reverse"),
    //                 matches.is_present("verbose-read-names"),
    //                 if matches.is_present("insert-size") {
    //                     Some(value_t!(matches, "insert-size", usize).unwrap())
    //                 } else {
    //                     None
    //                 },
    //                 if matches.is_present("std-dev") {
    //                     Some(value_t!(matches, "std-dev", usize).unwrap())
    //                 } else {
    //                     None
    //                 },
    //             )
    //         }
    //         ("bam", Some(matches)) => {
    //             bam::collapse_reads_to_fragments::call_consensus_reads_from_paths(
    //                 matches.value_of("bam").unwrap(),
    //                 matches.value_of("consensus-fq1").unwrap(),
    //                 matches.value_of("consensus-fq2").unwrap(),
    //                 matches.value_of("consensus-fq-se").unwrap(),
    //                 matches.value_of("skipped-bam").unwrap(),
    //                 matches.is_present("verbose-read-names"),
    //             )
    //         }
    //         _ => unreachable!(),
    //     },
    //     ("sequence-stats", Some(matches)) => sequences_stats::stats(matches.is_present("fastq")),
    //     // This cannot be reached, since the matches step of
    //     // clap assures that a valid subcommand is provided
    //     _ => unreachable!(),
    // }
}
