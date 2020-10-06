use std::error::Error;

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

    let header = rdr.headers()?;
    unimplemented!()
}
