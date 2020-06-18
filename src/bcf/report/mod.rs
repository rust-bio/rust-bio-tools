mod alignment_reader;
mod create_report;
mod fasta_reader;
mod static_reader;

use crate::bcf::report::create_report::make_report;
use std::error::Error;
use std::io::{stdout, Write};
use std::path::Path;
use tera::{Context, Tera};

pub fn report(vcf: &str, fasta: &str, bam: &str) -> Result<(), Box<dyn Error>> {
    let mut templates = Tera::default();
    templates
        .add_raw_template("report.html.tera", include_str!("report.html.tera"))
        .unwrap();
    let mut context = Context::new();
    context.insert(
        "variants",
        &make_report(Path::new(vcf), Path::new(fasta), Path::new(bam))?,
    );

    let html = templates.render("report.html.tera", &context).unwrap();
    stdout().write(html.as_bytes())?;

    Ok(())
}
