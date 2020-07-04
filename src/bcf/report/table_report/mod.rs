mod alignment_reader;
mod create_report_table;
mod fasta_reader;
mod static_reader;

use crate::bcf::report::table_report::create_report_table::make_table_report;
use std::error::Error;
use std::io::{stdout, Write};
use std::path::Path;
use tera::{Context, Tera};

pub fn table_report(vcf: &str, fasta: &str, bam: &str) -> Result<(), Box<dyn Error>> {
    let mut templates = Tera::default();
    templates
        .add_raw_template(
            "table_report.html.tera",
            include_str!("report_table.html.tera"),
        )
        .unwrap();
    let mut context = Context::new();
    context.insert(
        "variants",
        &make_table_report(Path::new(vcf), Path::new(fasta), Path::new(bam))?,
    );

    let html = templates
        .render("table_report.html.tera", &context)
        .unwrap();
    stdout().write(html.as_bytes())?;

    Ok(())
}
