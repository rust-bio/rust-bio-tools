mod alignment_reader;
pub mod create_report_table;
mod fasta_reader;
mod static_reader;

use crate::bcf::report::table_report::create_report_table::make_table_report;
use clap::Values;
use itertools::__std_iter::FromIterator;
use std::error::Error;
use std::fs;
use std::path::Path;

#[allow(clippy::too_many_arguments)]
pub fn table_report(
    vcf: &str,
    fasta: &str,
    bam: &[(String, String)],
    output_path: &str,
    sample: &str,
    info: Option<Values>,
    format: Option<Values>,
    max_read_depth: u32,
) -> Result<(), Box<dyn Error>> {
    let detail_path = output_path.to_owned() + "/details/" + sample;
    fs::create_dir(Path::new(&detail_path)).unwrap_or_else(|_| {
        panic!(
            "Could not create directory for table report files at location: {:?}",
            detail_path
        )
    });

    let info_strings = if let Some(value) = info {
        let strings = Vec::from_iter(value);
        Some(strings)
    } else {
        None
    };

    let format_strings = if let Some(value) = format {
        let strings = Vec::from_iter(value);
        Some(strings)
    } else {
        None
    };

    Ok(make_table_report(
        Path::new(vcf),
        Path::new(fasta),
        bam,
        info_strings,
        format_strings,
        sample.to_owned(),
        output_path,
        max_read_depth,
    )?)
}
