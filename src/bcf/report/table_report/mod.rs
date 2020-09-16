mod alignment_reader;
pub mod create_report_table;
mod fasta_reader;
mod static_reader;

use crate::bcf::report::table_report::create_report_table::make_table_report;
use chrono::{DateTime, Local};
use clap::Values;
use itertools::__std_iter::FromIterator;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use tera::{Context, Tera};

pub fn table_report(
    vcf: &str,
    fasta: &str,
    bam: &str,
    output_path: &str,
    sample: &str,
    info: Option<Values>,
    format: Option<Values>,
) -> Result<(), Box<dyn Error>> {
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

    let (reports, ann_field_description) = make_table_report(
        Path::new(vcf),
        Path::new(fasta),
        Path::new(bam),
        info_strings,
        format_strings,
    )?;

    let detail_path = output_path.to_owned() + "/details/" + sample;
    fs::create_dir(Path::new(&detail_path))?;

    let local: DateTime<Local> = Local::now();

    for (gene, report_data) in reports {
        let mut templates = Tera::default();
        templates
            .add_raw_template(
                "table_report.html.tera",
                include_str!("report_table.html.tera"),
            )
            .unwrap();
        let mut context = Context::new();
        context.insert("variants", &report_data);
        context.insert("gene", &gene);
        context.insert("description", &ann_field_description);
        context.insert("sample", &sample);
        context.insert("time", &local.format("%a %b %e %T %Y").to_string());
        context.insert("version", &env!("CARGO_PKG_VERSION"));

        let html = templates
            .render("table_report.html.tera", &context)
            .unwrap();
        let filepath = detail_path.clone() + "/" + &gene + ".html";
        let mut file = File::create(filepath)?;
        file.write_all(html.as_bytes())?;
    }

    Ok(())
}
