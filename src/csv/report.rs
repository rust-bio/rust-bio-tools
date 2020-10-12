use chrono::{DateTime, Local};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::iter::FromIterator;
use std::str::FromStr;
use tera::{Context, Tera};
use std::path::Path;
use std::fs;

pub(crate) fn csv_report(
    csv_path: &str,
    output_path: &str,
    rows_per_page: usize,
    separator: &str,
    sort_column: Option<&str>,
    ascending: Option<bool>,
) -> Result<(), Box<dyn Error>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(separator.as_bytes()[0])
        .from_path(csv_path)?;

    let header = rdr.headers()?.clone();
    let titles = Vec::from_iter(header.iter());
    let mut table = Vec::new();
    for res in rdr.records() {
        let row = res?;
        let mut table_entry = HashMap::new();
        for (i, tile) in titles.iter().enumerate() {
            table_entry.insert(tile.to_owned(), row[i].to_owned());
        }
        table.push(table_entry);
    }

    match (sort_column, ascending) {
        (Some(column), Some(true)) => table.sort_by(|a, b| {
            match (
                f32::from_str(a.get(column).unwrap()),
                f32::from_str(b.get(column).unwrap()),
            ) {
                (Ok(float_a), Ok(float_b)) => float_a.partial_cmp(&float_b).unwrap(),
                _ => a.get(column).cmp(&b.get(column)),
            }
        }),
        (Some(column), Some(false)) => table.sort_by(|a, b| {
            match (
                f32::from_str(a.get(column).unwrap()),
                f32::from_str(b.get(column).unwrap()),
            ) {
                (Ok(float_a), Ok(float_b)) => float_b.partial_cmp(&float_a).unwrap(),
                _ => a.get(column).cmp(&b.get(column)),
            }
        }),
        (_, _) => {}
    }

    let pages = if table.len() % rows_per_page == 0 {
        (table.len() / rows_per_page) - 1
    } else {
        table.len() / rows_per_page
    };

    let index_path = output_path.to_owned() + "/indexes/";
    fs::create_dir(Path::new(&index_path)).unwrap_or_else(|_| {
        panic!(
            "Could not create directory for report files at location: {:?}",
            index_path
        )
    });

    for i in 0..pages + 1 {
        let current_table = if i != pages {
            &table[(i * rows_per_page)..((i + 1) * rows_per_page)] // get genes for current page
        } else {
            &table[(i * rows_per_page)..] // get genes for last page
        };

        let page = i + 1;

        let mut templates = Tera::default();
        templates.add_raw_template("csv_report.html.tera", include_str!("csv_report.html.tera"))?;
        let mut context = Context::new();
        context.insert("table", &current_table);
        context.insert("titles", &titles);
        context.insert("current_page", &page);
        context.insert("pages", &(pages + 1));
        let local: DateTime<Local> = Local::now();
        context.insert("time", &local.format("%a %b %e %T %Y").to_string());
        context.insert("version", &env!("CARGO_PKG_VERSION"));

        let html = templates.render("csv_report.html.tera", &context)?;

        let file_path = output_path.to_owned() + "/indexes/index" + &page.to_string() + ".html";
        let mut file = File::create(file_path)?;
        file.write_all(html.as_bytes())?;
    }
    Ok(())
}
