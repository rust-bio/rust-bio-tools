use anyhow::Result;
use bio::io::fastq;
use flate2::bufread::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use regex::Regex;
use std::fs;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

fn writer<P: AsRef<Path>>(path: P) -> Result<fastq::Writer<Box<dyn std::io::Write>>> {
    let w: Box<dyn Write> = if path.as_ref().extension().unwrap() == "gz" {
        Box::new(
            fs::File::create(&path)
                .map(BufWriter::new)
                .map(|w| GzEncoder::new(w, Compression::default()))?,
        )
    } else {
        Box::new(fs::File::create(&path).map(BufWriter::new)?)
    };
    Ok(fastq::Writer::new(w))
}

fn reader<P: AsRef<Path>>(path: P) -> Result<fastq::Reader<BufReader<Box<dyn std::io::Read>>>> {
    let r: Box<dyn Read> = if path.as_ref().extension().unwrap() == "gz" {
        Box::new(
            fs::File::open(&path)
                .map(BufReader::new)
                .map(MultiGzDecoder::new)?,
        )
    } else {
        Box::new(fs::File::open(&path).map(BufReader::new)?)
    };
    Ok(fastq::Reader::new(r))
}

//TODO That is really bad for just adding a suffix to the file name
fn build_fastq_out(fastq_in: &PathBuf) -> Result<PathBuf> {
    let suffix = "reheadered.fastq.gz";
    let mut fastq_out = if fastq_in.extension().unwrap() == "gz" {
        fastq_in.to_str().unwrap().strip_suffix("fastq.gz").unwrap()
    } else {
        fastq_in.to_str().unwrap().strip_suffix("fastq").unwrap()
    }
    .to_owned();
    fastq_out.push_str(suffix);
    Ok(Path::new(&fastq_out).to_path_buf())
}

pub fn reformat_header(fastqs: Vec<PathBuf>, desc_regex: &str, desc_format: &str) -> Result<()> {
    let regex = Regex::new(desc_regex)?;
    let capture_groups = parse_format_fields(desc_regex)?;
    for fastq in fastqs {
        process_fastq(&fastq, &regex, desc_format, &capture_groups)?;
    }
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

fn process_fastq(
    fastq_in: &PathBuf,
    desc_regex: &Regex,
    desc_format: &str,
    capture_groups: &[String],
) -> Result<()> {
    let reader = reader(&fastq_in)?;
    let fastq_out = build_fastq_out(fastq_in)?;
    let mut writer = writer(fastq_out)?;
    for result in reader.records() {
        let record = result?;
        let description_opt = record.desc();
        let description_formatted = match description_opt {
            Some(description) => {
                let caps = desc_regex.captures(description).unwrap();
                for caps in desc_regex.captures_iter(description) {
                    dbg!(caps);
                    //let desc_updated = format!(desc_format.to_string(), umi=caps.name("umi").unwrap());
                }
                None
            }
            None => None,
        };
        let record_out = fastq::Record::with_attrs(
            record.id(),
            description_formatted,
            record.seq(),
            record.qual(),
        );
        writer.write_record(&record_out)?;
    }
    Ok(())
}
