use chrono::{DateTime, Local};
use derive_new::new;
use itertools::Itertools;
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
    let mut num_plot_data = HashMap::new();

    for title in &titles {
        match is_numeric.get(title) {
            Some(true) => {
                let plot = num_plot(table.clone(), title.to_string());
                num_plot_data.insert(title, plot);
            }
            Some(false) => {
                let plot = nominal_plot(table.clone(), title.to_string());
                plot_data.insert(title, plot);
            }
            _ => unreachable!(),
        };
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
        templates.add_raw_template("plot.js.tera", include_str!("plot.js.tera"))?;
        let mut context = Context::new();
        match is_numeric.get(title) {
            Some(true) => {
                context.insert(
                    "table",
                    &json!(num_plot_data.get(title).unwrap()).to_string(),
                );
                context.insert("num", &true);
            }
            Some(false) => {
                context.insert("table", &json!(plot_data.get(title).unwrap()).to_string());
                context.insert("num", &false);
            }
            _ => unreachable!(),
        }
        context.insert("title", &title);
        let js = templates.render("plot.js.tera", &context)?;

        let file_path = plot_path.to_owned() + title + ".js";
        let mut file = File::create(file_path)?;
        file.write_all(js.as_bytes())?;
    }

    let index_path = output_path.to_owned() + "/indexes/";
    fs::create_dir(Path::new(&index_path)).unwrap_or_else(|_| {
        panic!(
            "Could not create directory for report files at location: {:?}",
            index_path
        )
    });

    let prefixes = make_prefixes(table.clone(), titles.clone(), rows_per_page);

    let prefix_path = output_path.to_owned() + "/prefixes/";
    fs::create_dir(Path::new(&prefix_path)).unwrap_or_else(|_| {
        panic!(
            "Could not create directory for report files at location: {:?}",
            prefix_path
        )
    });

    for title in &titles {
        if let Some(prefix_table) = prefixes.get(title.to_owned()) {
            let mut templates = Tera::default();
            templates.add_raw_template(
                "prefix_table.html.tera",
                include_str!("prefix_table.html.tera"),
            )?;
            let mut context = Context::new();
            context.insert("title", title);
            context.insert("table", prefix_table);
            let html = templates.render("prefix_table.html.tera", &context)?;

            let file_path = output_path.to_owned() + "/prefixes/" + title + ".html";
            let mut file = File::create(file_path)?;
            file.write_all(html.as_bytes())?;

            let title_path = prefix_path.to_owned() + "/" + title + "/";
            fs::create_dir(Path::new(&title_path)).unwrap_or_else(|_| {
                panic!(
                    "Could not create directory for report files at location: {:?}",
                    title_path
                )
            });

            for (prefix, values) in prefix_table {
                let mut templates = Tera::default();
                templates.add_raw_template(
                    "lookup_table.html.tera",
                    include_str!("lookup_table.html.tera"),
                )?;
                let mut context = Context::new();
                context.insert("title", title);
                context.insert("values", values);
                let html = templates.render("lookup_table.html.tera", &context)?;

                let file_path = title_path.to_owned() + prefix + ".html";
                let mut file = File::create(file_path)?;
                file.write_all(html.as_bytes())?;
            }
        }
    }

    for (i, current_table) in table.chunks(rows_per_page).enumerate() {
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

fn num_plot(table: Vec<HashMap<String, String>>, column: String) -> Vec<BinnedPlotRecord> {
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
    let bins = 20;
    let step = (max - min) / bins as f32;
    let mut binned_data = HashMap::new();
    let mut bin_borders = HashMap::new();
    for val in values {
        for i in 0..bins {
            let lower_bound = min + i as f32 * step;
            let upper_bound = lower_bound + step;
            let bin_name = String::from("bin") + &i.to_string();
            bin_borders.insert(bin_name.to_owned(), (lower_bound, upper_bound));
            let entry = binned_data.entry(bin_name.to_owned()).or_insert_with(|| 0);
            if ((i < (bins - 1) && val < upper_bound) || (i < bins && val <= upper_bound))
                && val >= lower_bound
            {
                *entry += 1;
            }
        }
    }
    if nan > 0 {
        bin_borders.insert(
            String::from("bin") + &bins.to_string(),
            (f32::NAN, f32::NAN),
        );
        binned_data.insert(String::from("bin") + &bins.to_string(), nan);
    }
    let mut plot_data = Vec::new();
    for (name, v) in binned_data {
        let (lower_bound, upper_bound) = bin_borders.get(&name).unwrap();
        let plot_record = BinnedPlotRecord {
            bin_start: *lower_bound,
            value: v,
            bin_end: *upper_bound,
        };
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

fn make_prefixes(
    table: Vec<HashMap<String, String>>,
    titles: Vec<&str>,
    rows_per_page: usize,
) -> HashMap<String, HashMap<String, Vec<(String, usize)>>> {
    let mut title_map = HashMap::new();
    for (i, partial_table) in table.chunks(rows_per_page).enumerate() {
        let page = i + 1;
        let prefix_len = 3;
        for row in partial_table {
            for key in &titles {
                let value = &row[key.to_owned()];
                let entry = value.split_whitespace().take(1).collect_vec()[0];
                if entry.len() >= prefix_len {
                    let prefix = entry.chars().take(prefix_len).collect::<String>();
                    let prefix_map = title_map
                        .entry(key.to_string())
                        .or_insert_with(HashMap::new);
                    let values = prefix_map.entry(prefix).or_insert_with(Vec::new);
                    values.push((value.to_owned(), page));
                }
            }
        }
        // write stuff to output map with page like so: HashMap<column_title, HashMap<prefix, Vec<(value, page)>>>
    }
    title_map
}

#[derive(new, Serialize, Debug, Clone)]
struct PlotRecord {
    key: String,
    value: u32,
}

#[derive(new, Serialize, Debug, Clone)]
struct BinnedPlotRecord {
    bin_start: f32,
    bin_end: f32,
    value: u32,
}
