use chrono::{DateTime, Local};
use std::collections::HashMap;
use std::error::Error;
use std::iter::FromIterator;
use tera::{Context, Tera};
use std::fs::File;
use std::io::Write;

pub(crate) fn csv_report(
    csv_path: &str,
    output_path: &str,
    rows_per_page: u32,
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
        (Some(column), Some(true)) => {
            table.sort_by(|a,b| a.get(column).cmp(&b.get(column)))
        }
        (Some(column), Some(false)) => {
            table.sort_by(|a,b| b.get(column).cmp(&a.get(column)))
        }
        (_,_) => {}
    }

    let mut templates = Tera::default();
    templates.add_raw_template("csv_report.html.tera", include_str!("csv_report.html.tera"))?;
    let mut context = Context::new();
    context.insert("table", &table);
    context.insert("titles", &titles);
    let local: DateTime<Local> = Local::now();
    context.insert("time", &local.format("%a %b %e %T %Y").to_string());
    context.insert("version", &env!("CARGO_PKG_VERSION"));

    let html = templates.render("csv_report.html.tera", &context)?;

    let file_path = output_path.to_owned() + "/index.html";
    let mut file = File::create(file_path)?;
    file.write_all(html.as_bytes())?;

    Ok(())
}
