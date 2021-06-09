use crate::bcf::report::table_report::create_report_table::create_report_data;
use crate::bcf::report::table_report::create_report_table::manipulate_json;
use anyhow::Result;
use chrono::{DateTime, Local};
use itertools::Itertools;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::str::FromStr;
use tera::{Context, Tera};

pub(crate) fn plot_bam(
    bam_path: &str,
    fasta_path: &str,
    region: &str,
    max_read_depth: u32,
    output_path: &str,
) -> Result<()> {
    let splitted_region = region.split(":").collect_vec();
    let chrom = splitted_region[0];
    let span = splitted_region[1].split("-").collect_vec();
    assert_eq!(2, span.len());
    let start = u64::from_str(span[0])?;
    let end = u64::from_str(span[1])?;

    let (content, max_rows) = create_report_data(
        Path::new(fasta_path),
        None,
        Path::new(bam_path),
        chrom.to_owned(),
        start,
        end,
        max_read_depth,
    )?;
    let visualization = manipulate_json(content, start, end, max_rows)?;

    let mut templates = Tera::default();
    templates.add_raw_template("bam_plot.html.tera", include_str!("bam_plot.html.tera"))?;
    let mut context = Context::new();
    let local: DateTime<Local> = Local::now();
    context.insert("time", &local.format("%a %b %e %T %Y").to_string());
    context.insert("version", &env!("CARGO_PKG_VERSION"));
    context.insert("specs", &visualization);
    context.insert("bam", &Path::new(bam_path).iter().last().unwrap().to_str().unwrap());
    context.insert("chrom", &chrom);
    context.insert("start", &start);
    context.insert("end", &end);

    let html = templates.render("bam_plot.html.tera", &context)?;
    let filepath = Path::new(output_path).join("plot.html");
    let mut file = File::create(filepath)?;
    file.write_all(html.as_bytes())?;

    Ok(())
}
