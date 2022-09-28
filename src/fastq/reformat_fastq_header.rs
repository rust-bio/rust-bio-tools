use anyhow::Result;
use bio::io::fastq;
use handlebars::Handlebars;
use regex::Regex;
use std::collections::HashMap;
use std::io;

pub fn reformat_header(desc_regex: &str, desc_format: &str) -> Result<()> {
    let regex = Regex::new(desc_regex)?;
    let capture_groups = parse_format_fields(desc_regex)?;
    process_fastq(&regex, desc_format, &capture_groups)?;
    Ok(())
}

fn parse_format_fields(desc_regex: &str) -> Result<Vec<String>> {
    let re = Regex::new("<([A-Za-z]+)>")?;
    let caps = re.captures(desc_regex).unwrap();
    let capture_fields: Vec<String> = (1..caps.len())
        .map(|i| caps.get(i).unwrap().as_str().to_string())
        .collect();
    Ok(capture_fields)
}

fn process_fastq(desc_regex: &Regex, desc_format: &str, capture_groups: &[String]) -> Result<()> {
    let reader = fastq::Reader::new(io::stdin());
    let mut writer = fastq::Writer::new(io::stdout());
    for result in reader.records() {
        let reg = Handlebars::new();
        let record = result?;
        let description_opt = record.desc();
        let description_extended = match description_opt {
            Some(description) => {
                let mut description_extended = description.to_owned();
                for caps in desc_regex.captures_iter(description) {
                    let mut field_replacements = HashMap::new();
                    for field in capture_groups {
                        field_replacements.insert(field, caps.name(field).unwrap().as_str());
                    }
                    let rendered_entry = reg.render_template(desc_format, &field_replacements)?;
                    //TODO That's dirty
                    description_extended.push(' ');
                    description_extended.push_str(&rendered_entry);
                }
                Some(description_extended)
            }
            None => None,
        };
        let record_out = fastq::Record::with_attrs(
            record.id(),
            description_extended.as_deref(),
            record.seq(),
            record.qual(),
        );
        writer.write_record(&record_out)?;
    }
    Ok(())
}
