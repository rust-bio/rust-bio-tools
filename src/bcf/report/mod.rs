use crate::bcf::report::oncoprint::WriteErr;
use anyhow::{Context, Result};
use itertools::Itertools;
use smart_open::smart_open;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;

pub mod oncoprint;
pub mod table_report;

pub fn embed_js(
    output_path: &str,
    vcf_report: bool,
    custom_table_report_js: Option<&str>,
    custom_js_files: Vec<String>,
) -> Result<()> {
    let js_path = output_path.to_owned() + "/js/";
    fs::create_dir(Path::new(&js_path)).context(WriteErr::CantCreateDir {
        dir_path: js_path.to_owned(),
    })?;
    let mut files = vec![
        (
            "bootstrap.bundle.min.js",
            include_str!("js/bootstrap.bundle.min.js"),
        ),
        ("jquery.min.js", include_str!("js/jquery.min.js")),
        ("lz-string.min.js", include_str!("js/lz-string.min.js")),
        (
            "bootstrap-table.min.js",
            include_str!("js/bootstrap-table.min.js"),
        ),
        ("vega.min.js", include_str!("js/vega.min.js")),
        ("vega-lite.min.js", include_str!("js/vega-lite.min.js")),
        ("vega-embed.min.js", include_str!("js/vega-embed.min.js")),
    ];
    let vcf_report_files = vec![
        ("jsonm.min.js", include_str!("js/jsonm.min.js")),
        ("table-report.js", include_str!("js/table-report.js")),
        ("report.js", include_str!("js/report.js")),
        ("gene-report.js", include_str!("js/gene-report.js")),
    ];
    if vcf_report {
        files.extend(vcf_report_files.iter());
        if let Some(path) = custom_table_report_js {
            let custom_file_string = smart_open(path).context("Unable to open custom JS file")?;
            let mut out_file = File::create(js_path.to_owned() + "table-report.js")?;
            out_file.write_all(custom_file_string.as_bytes())?;
        } else {
            files.push(("table-report.js", include_str!("js/table-report.js")))
        }
    }
    for (name, file) in files {
        let mut out_file = File::create(js_path.to_owned() + name)?;
        out_file.write_all(file.as_bytes())?;
    }

    for file in custom_js_files {
        let file_name = file
            .split('/')
            .collect_vec()
            .pop()
            .context(format!("Unable to extract file name from path: {}", file))?;
        let custom_file_string = smart_open(&file).context("Unable to open JS file")?;
        let mut out_file = File::create(js_path.to_owned() + file_name)?;
        out_file.write_all(custom_file_string.as_bytes())?;
    }

    Ok(())
}

pub fn embed_css(output_path: &str, vcf_report: bool) -> Result<()> {
    let css_path = output_path.to_owned() + "/css/";
    fs::create_dir(Path::new(&css_path)).context(WriteErr::CantCreateDir {
        dir_path: css_path.to_owned(),
    })?;
    let mut files = vec![
        ("bootstrap.min.css", include_str!("css/bootstrap.min.css")),
        (
            "bootstrap-table.min.css",
            include_str!("css/bootstrap-table.min.css"),
        ),
    ];
    let vcf_report_files = vec![("oncoprint.css", include_str!("css/oncoprint.css"))];
    let csv_report_files = vec![("csv_report.css", include_str!("../../csv/csv_report.css"))];
    if vcf_report {
        files.extend(vcf_report_files.iter());
    } else {
        files.extend(csv_report_files.iter());
    }
    for (name, file) in files {
        let mut out_file = File::create(css_path.to_owned() + name)?;
        out_file.write_all(file.as_bytes())?;
    }
    Ok(())
}

pub fn embed_html(output_path: &str) -> Result<()> {
    let files = vec![("index.html", include_str!("html/index.html"))];
    for (name, file) in files {
        let out_path = output_path.to_owned() + "/" + name;
        let mut out_file = File::create(&out_path).context(WriteErr::CantCreateDir {
            dir_path: out_path.to_owned(),
        })?;
        out_file.write_all(file.as_bytes())?;
    }
    Ok(())
}
