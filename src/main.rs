//! Documentation for Rust Bio Tools
use anyhow::{Context, Result};
use itertools::Itertools;
use log::LevelFilter;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use structopt::StructOpt;

use cli::Command::*;

pub mod bam;
pub mod bcf;
mod cli;
pub mod common;
pub mod csv;
pub mod fastq;
pub mod sequences_stats;

fn main() -> Result<()> {
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
            with_filter,
        } => bcf::to_txt::to_txt(
            info.iter().map(|s| s as &str).collect_vec().as_slice(),
            format.iter().map(|s| s as &str).collect_vec().as_slice(),
            genotypes,
            with_filter,
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
        CsvReport {
            csv_path,
            rows_per_page,
            sort_column,
            sort_order,
            separator,
            formatter,
            pin_until,
            output_path,
        } => {
            if !Path::new(&output_path).exists() {
                fs::create_dir_all(Path::new(&output_path))?;
            }
            bcf::report::embed_js(&output_path, false, None, vec![])?;
            bcf::report::embed_css(&output_path, false)?;
            bcf::report::embed_html(&output_path)?;

            let order = match sort_order.as_str() {
                "ascending" => Some(true),
                "descending" => Some(false),
                _ => None,
            };

            csv::report::csv_report(
                &csv_path,
                &output_path,
                rows_per_page as usize,
                separator,
                sort_column.as_deref(),
                order,
                formatter.as_deref(),
                pin_until.as_deref(),
            )?
        }
        PlotBam {
            bam_path,
            reference,
            region,
            max_read_depth,
        } => bam::plot::plot_bam::plot_bam(&bam_path, reference, &region, max_read_depth)?,
        VcfReport {
            fasta,
            vcfs,
            bams,
            cells,
            max_read_depth,
            infos,
            formats,
            plot_info,
            custom_js_template,
            custom_js_files,
            tsv,
            threads,
            annotation_field,
            output_path,
        } => {
            let mut sample_calls = HashMap::new();
            let mut bam_paths = HashMap::new();
            if !Path::new(&output_path).exists() {
                fs::create_dir(Path::new(&output_path)).context(format!(
                    "Couldn't create output directory at {}. Please make sure the path exists.",
                    output_path
                ))?;
            }
            let js_files_vec = custom_js_files
                .clone()
                .map_or_else(Vec::new, |values| values.into_iter().collect());
            let js_file_names = if let Some(files) = custom_js_files {
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
            bcf::report::embed_js(
                &output_path,
                true,
                custom_js_template.as_deref(),
                js_files_vec,
            )?;
            bcf::report::embed_css(&output_path, true)?;
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
                    &annotation_field,
                )
                .unwrap_or_else(|e| {
                    panic!("Failed building table report for sample {}. {}", sample, e)
                });
            });

            bcf::report::oncoprint::oncoprint(
                &sample_calls,
                &output_path,
                cells,
                tsv.as_deref(),
                plot_info,
                &annotation_field,
            )?
        }
        VcfSplit { input, output } => bcf::split::split(input, output.as_ref())?,
        CollapseReadsToFragments { cmd } => match cmd {
            cli::CollapseReadsToFragmentsSubcommand::Fastq {
                fq1,
                fq2,
                consensus_fq1,
                consensus_fq2,
                consensus_fq3,
                umi_len,
                max_seq_dist,
                max_umi_dist,
                umi_on_reverse,
                verbose_read_names,
                insert_size,
                std_dev,
            } => fastq::collapse_reads_to_fragments::call_consensus_reads_from_paths(
                fq1,
                fq2,
                consensus_fq1,
                consensus_fq2,
                consensus_fq3,
                umi_len,
                max_seq_dist,
                max_umi_dist,
                umi_on_reverse,
                verbose_read_names,
                insert_size,
                std_dev,
            )?,
            cli::CollapseReadsToFragmentsSubcommand::Bam {
                bam,
                consensus_fq1,
                consensus_fq2,
                consensus_fq_se,
                skipped_bam,
                verbose_read_names,
            } => bam::collapse_reads_to_fragments::call_consensus_reads_from_paths(
                bam,
                consensus_fq1,
                consensus_fq2,
                consensus_fq_se,
                skipped_bam,
                verbose_read_names,
            )?,
        },
        BamAnonymize {
            bam,
            input_ref,
            output_bam,
            output_ref,
            chr,
            start,
            end,
            keep_only_pairs,
        } => bam::anonymize_reads::anonymize_reads(
            bam,
            input_ref,
            output_bam,
            output_ref,
            chr,
            start - 1..end - 1,
            keep_only_pairs,
        )?,
        SequenceStats { fastq } => sequences_stats::stats(fastq)?,
    }
    Ok(())
}
