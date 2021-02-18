mod alignment_reader;
pub mod create_report_table;
mod fasta_reader;
mod static_reader;

use crate::bcf::report::table_report::create_report_table::make_table_report;
use anyhow::{Context, Result};
use std::fs;
use std::path::Path;

#[allow(clippy::too_many_arguments)]
pub fn table_report(
    vcf: &str,
    fasta: &str,
    bam: &[(String, String)],
    output_path: &str,
    sample: &str,
    info_strings: Option<Vec<String>>,
    format_strings: Option<Vec<String>>,
    max_read_depth: u32,
    js_files: Vec<String>,
) -> Result<()> {
    let detail_path = output_path.to_owned() + "/details/" + sample;
    fs::create_dir(Path::new(&detail_path)).with_context(|| {
        format!(
            "Could not create directory for table report files at location: {}",
            detail_path
        )
    })?;

    let plot_path = detail_path + "/plots/";
    fs::create_dir(Path::new(&plot_path)).with_context(|| {
        format!(
            "Could not create directory for table report plots at location: {}",
            plot_path
        )
    })?;

    make_table_report(
        Path::new(vcf),
        Path::new(fasta),
        bam,
        info_strings,
        format_strings,
        sample.to_owned(),
        output_path,
        max_read_depth,
        js_files,
    )
}
