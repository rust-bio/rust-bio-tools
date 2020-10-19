use chrono::{DateTime, Local};
use derive_new::new;
use serde_derive::Serialize;
use serde_json::json;
use simple_excel_writer::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::iter::FromIterator;
use std::path::Path;
use std::str::FromStr;
use tera::{Context, Tera};

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
    let mut numeric = HashMap::new();
    let mut non_numeric = HashMap::new();
    for res in rdr.records() {
        let row = res?;
        let mut table_entry = HashMap::new();
        for (i, tile) in titles.iter().enumerate() {
            table_entry.insert(tile.to_string(), row[i].to_owned());
            match f32::from_str(&row[i]) {
                Ok(_) => {
                    let num = numeric.entry(tile.to_owned()).or_insert_with(|| 0);
                    *num += 1;
                }
                _ => {
                    let no_num = non_numeric.entry(tile.to_owned()).or_insert_with(|| 0);
                    *no_num += 1;
                }
            }
        }
        table.push(table_entry);
    }

    let mut is_numeric = HashMap::new();
    for title in &titles {
        let is_num = match (numeric.get(title), non_numeric.get(title)) {
            (Some(num), Some(no_num)) => num > no_num,
            (None, Some(_)) => false,
            (Some(_), None) => true,
            _ => unreachable!(),
        };
        is_numeric.insert(title.to_owned(), is_num);
    }

    let mut plot_data = HashMap::new();
    for title in &titles {
        let data = match is_numeric.get(title) {
            Some(true) => num_plot(table.clone(), title.to_string()),
            Some(false) => nominal_plot(table.clone(), title.to_string()),
            _ => unreachable!(),
        };
        plot_data.insert(title, data);
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

    let mut wb = Workbook::create(&(output_path.to_owned() + "/report.xlsx"));
    let mut sheet = wb.create_sheet("Report");
    for _ in 1..titles.len() {
        sheet.add_column(Column { width: 50.0 });
    }

    wb.write_sheet(&mut sheet, |sheet_writer| {
        let sw = sheet_writer;
        let mut title_row = Row::new();
        for title in titles.clone() {
            title_row.add_cell(title);
        }
        sw.append_row(title_row)?;
        for row in table.clone() {
            let mut excel_row = Row::new();
            for title in titles.clone() {
                excel_row.add_cell(row.get(title).unwrap().as_str());
            }
            sw.append_row(excel_row)?;
        }
        Ok(())
    })?;

    wb.close()?;

    let pages = if table.len() % rows_per_page == 0 {
        (table.len() / rows_per_page) - 1
    } else {
        table.len() / rows_per_page
    };

    let plot_path = output_path.to_owned() + "/plots/";
    fs::create_dir(Path::new(&plot_path)).unwrap_or_else(|_| {
        panic!(
            "Could not create directory for plot files at location: {:?}",
            plot_path
        )
    });

    for title in &titles {
        let mut templates = Tera::default();
        templates.add_raw_template("plot.html.tera", include_str!("plot.html.tera"))?;
        let mut context = Context::new();
        context.insert("table", &json!(plot_data.get(title)).to_string());
        context.insert("title", &title);
        let html = templates.render("plot.html.tera", &context)?;

        let file_path = plot_path.to_owned() + title + ".html";
        let mut file = File::create(file_path)?;
        file.write_all(html.as_bytes())?;
    }

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

fn num_plot(table: Vec<HashMap<String, String>>, column: String) -> Vec<PlotRecord> {
    let mut values = Vec::new();
    let mut nan = 0;
    for row in table {
        match f32::from_str(row.get(&column).unwrap()) {
            Ok(val) => values.push(val.to_owned()),
            _ => nan += 1,
        }
    }
    let min = values.iter().fold(f32::INFINITY, |a, &b| a.min(b));
    let max = values.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b));
    let bins = 6;
    let step = (max - min) / bins as f32;
    let mut binned_data = HashMap::new();
    for val in values {
        for i in 0..bins {
            let lower_bound = min + i as f32 * step;
            let upper_bound = lower_bound + step;
            let bin_name = lower_bound.to_string() + " - " + &upper_bound.to_string();
            let entry = binned_data.entry(bin_name).or_insert_with(|| 0);
            if ((i < (bins - 1) && val < upper_bound) || (i < bins && val <= upper_bound))
                && val >= lower_bound
            {
                *entry += 1;
            }
        }
    }
    if nan > 0 {
        binned_data.insert(String::from("NaN"), nan);
    }

    let mut plot_data = Vec::new();
    for (k, v) in binned_data {
        let plot_record = PlotRecord { key: k, value: v };
        plot_data.push(plot_record);
    }
    plot_data
}

fn nominal_plot(table: Vec<HashMap<String, String>>, column: String) -> Vec<PlotRecord> {
    let mut values = Vec::new();
    for row in table {
        let val = row.get(&column).unwrap();
        values.push(val.to_owned());
    }
    let mut count_values = HashMap::new();
    for v in values {
        let entry = count_values.entry(v.to_owned()).or_insert_with(|| 0);
        *entry += 1;
    }

    let mut plot_data = Vec::new();
    for (k, v) in count_values {
        let plot_record = PlotRecord { key: k, value: v };
        plot_data.push(plot_record);
    }

    if plot_data.len() > 10 {
        plot_data.sort_by(|a, b| b.value.cmp(&a.value));
        plot_data = plot_data.into_iter().take(10).collect();
    }

    plot_data
}

#[derive(new, Serialize, Debug)]
struct PlotRecord {
    key: String,
    value: u32,
}
