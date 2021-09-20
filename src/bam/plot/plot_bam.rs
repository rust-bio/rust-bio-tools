use crate::bcf::report::table_report::create_report_table::create_report_data;
use crate::bcf::report::table_report::create_report_table::manipulate_json;
use crate::common::Region;
use anyhow::Result;
use chrono::{DateTime, Local};
use itertools::Itertools;
use std::io;
use std::io::Write;
use std::path::Path;
use tera::{Context, Tera};

pub(crate) fn plot_bam<P: AsRef<Path>>(
    bam_paths: &[P],
    fasta_path: P,
    region: &Region,
    max_read_depth: u32,
) -> Result<()> {
    let mut plots = Vec::new();

    let Region { target, start, end } = region.clone();
    for bam_path in bam_paths {
        let content=
            create_report_data(&fasta_path, None, bam_path, region, max_read_depth)?;
        let visualization = manipulate_json(content, start, end)?;

        plots.push(visualization);
    }

    let bams = bam_paths
        .iter()
        .map(|b| b.as_ref().iter().last().unwrap().to_str().unwrap())
        .collect_vec();

    let mut templates = Tera::default();
    templates.add_raw_template("bam_plot.html.tera", include_str!("bam_plot.html.tera"))?;
    let mut context = Context::new();
    let local: DateTime<Local> = Local::now();
    context.insert("time", &local.format("%a %b %e %T %Y").to_string());
    context.insert("version", &env!("CARGO_PKG_VERSION"));
    context.insert("plots", &plots);
    context.insert("bams", &bams);
    context.insert("chrom", &target);
    context.insert("start", &start);
    context.insert("end", &end);

    let html = templates.render("bam_plot.html.tera", &context)?;
    io::stdout().write_all(html.as_bytes())?;

    Ok(())
}
