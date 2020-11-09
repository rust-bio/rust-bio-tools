use itertools::Itertools;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

pub mod oncoprint;
pub mod table_report;

pub fn embed_js(
    output_path: &str,
    custom_table_report_js: Option<&str>,
    custom_js_files: Vec<&str>,
) -> Result<(), Box<dyn Error>> {
    let js_path = output_path.to_owned() + "/js/";
    fs::create_dir(Path::new(&js_path)).unwrap_or_else(|_| {
        panic!(
            "Could not create directory for js files at location: {:?}",
            js_path
        )
    });
    let mut files = vec![
        ("vega.min.js", include_str!("js/vega.min.js")),
        ("vega-lite.min.js", include_str!("js/vega-lite.min.js")),
        ("vega-embed.min.js", include_str!("js/vega-embed.min.js")),
        ("jsonm.min.js", include_str!("js/jsonm.min.js")),
        ("jquery.min.js", include_str!("js/jquery.min.js")),
        (
            "bootstrap-table.min.js",
            include_str!("js/bootstrap-table.min.js"),
        ),
        (
            "bootstrap.bundle.min.js",
            include_str!("js/bootstrap.bundle.min.js"),
        ),
        ("report.js", include_str!("js/report.js")),
        ("gene-report.js", include_str!("js/gene-report.js")),
    ];
    if let Some(path) = custom_table_report_js {
        let mut file_string = String::new();
        let mut custom_file = File::open(path).expect("Unable to open custom JS file");
        custom_file
            .read_to_string(&mut file_string)
            .expect("Unable to read string");
        let mut out_file = File::create(js_path.to_owned() + "table-report.js")?;
        out_file.write_all(file_string.as_bytes())?;
    } else {
        files.push(("table-report.js", include_str!("js/table-report.js")))
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
            .unwrap_or_else(|| panic!("Unable to extract file name from path: {:?}", file));
        let mut file_string = String::new();
        let mut custom_file = File::open(file).expect("Unable to open JS file");
        custom_file
            .read_to_string(&mut file_string)
            .expect("Unable to read string");
        let mut out_file = File::create(js_path.to_owned() + file_name)?;
        out_file.write_all(file_string.as_bytes())?;
    }

    Ok(())
}

pub fn embed_css(output_path: &str) -> Result<(), Box<dyn Error>> {
    let css_path = output_path.to_owned() + "/css/";
    fs::create_dir(Path::new(&css_path)).unwrap_or_else(|_| {
        panic!(
            "Could not create directory for js files at location: {:?}",
            css_path
        )
    });
    let files = vec![
        ("bootstrap.min.css", include_str!("css/bootstrap.min.css")),
        (
            "bootstrap-table.min.css",
            include_str!("css/bootstrap-table.min.css"),
        ),
        ("oncoprint.css", include_str!("css/oncoprint.css")),
    ];
    for (name, file) in files {
        let mut out_file = File::create(css_path.to_owned() + name)?;
        out_file.write_all(file.as_bytes())?;
    }
    Ok(())
}

pub fn embed_html(output_path: &str) -> Result<(), Box<dyn Error>> {
    let files = vec![("index.html", include_str!("html/index.html"))];
    for (name, file) in files {
        let mut out_file = File::create(output_path.to_owned() + "/" + name).unwrap_or_else(|_| {
            panic!(
                "Could not create file for index html at location: {:?}",
                output_path.to_owned() + "/" + name
            )
        });
        out_file.write_all(file.as_bytes())?;
    }
    Ok(())
}
