//! Create a variant table from a VCF file.
//!
//! ## Usage:
//! ```bash
//! $ rbt vcf-to-txt --genotypes --fmt S --info T X SOMATIC < tests/test.vcf > tests/variant-table.txt
//! ```
//!
use anyhow::{bail, Result};
use derive_new::new;
use itertools::Itertools;
use rust_htslib::bcf;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::Read;
use std::io;
use std::io::Write;
use std::str;
use thiserror::Error;

#[derive(new)]
pub struct Writer {
    inner: io::BufWriter<io::Stdout>,
    #[new(value = "0")]
    field_count: usize,
}

impl Writer {
    fn write_integer(&mut self, value: i32) -> Result<()> {
        let fmt = if value.is_missing() {
            "".to_owned()
        } else {
            format!("{}", value)
        };
        self.write_field(fmt.as_bytes())
    }

    fn write_float(&mut self, value: f32) -> Result<()> {
        let fmt = if value.is_missing() {
            "".to_owned()
        } else {
            format!("{}", value)
        };
        self.write_field(fmt.as_bytes())
    }

    fn write_flag(&mut self, value: bool) -> Result<()> {
        self.write_field(format!("{}", value).as_bytes())
    }

    fn write_field(&mut self, value: &[u8]) -> Result<()> {
        if self.field_count > 0 {
            self.inner.write_all(b"\t")?;
        }
        self.inner.write_all(value)?;
        self.field_count += 1;
        Ok(())
    }

    fn newline(&mut self) -> Result<()> {
        self.inner.write_all(b"\n")?;
        self.field_count = 0;
        Ok(())
    }
}

const HEADER_COMMON: &[u8] = b"VARIANT";

pub fn to_txt(
    info_tags: &[&str],
    format_tags: &[&str],
    show_genotypes: bool,
    show_filter: bool,
) -> Result<()> {
    let mut reader = bcf::Reader::from_stdin()?;
    let mut writer = Writer::new(io::BufWriter::new(io::stdout()));

    let mut common_n = 5 + info_tags.len();
    if show_filter {
        common_n += 1
    }
    writer.write_field(HEADER_COMMON)?;
    for _ in 1..common_n {
        writer.write_field(HEADER_COMMON)?;
    }
    let show_samples = show_genotypes || !format_tags.is_empty();
    if show_samples {
        for sample in reader.header().samples() {
            writer.write_field(sample)?;
            for _ in 1..format_tags.len() + show_genotypes as usize {
                writer.write_field(sample)?;
            }
        }
    }
    writer.newline()?;
    writer.write_field(b"CHROM")?;
    writer.write_field(b"POS")?;
    writer.write_field(b"REF")?;
    writer.write_field(b"ALT")?;
    writer.write_field(b"QUAL")?;
    if show_filter {
        writer.write_field(b"FILTER")?;
    }
    for name in info_tags {
        writer.write_field(name.as_bytes())?;
    }
    if show_samples {
        for _ in 0..reader.header().sample_count() {
            if show_genotypes {
                writer.write_field(b"GT")?;
            }
            for name in format_tags {
                writer.write_field(name.as_bytes())?;
            }
        }
    }
    writer.newline()?;
    let mut rec = reader.empty_record();
    loop {
        match reader.read(&mut rec) {
            Some(Ok(())) => (),
            None => break,
            Some(Err(e)) => bail!(e),
        };
        let alleles = rec
            .alleles()
            .into_iter()
            .map(|a| a.to_owned())
            .collect_vec();
        for (i, allele) in alleles[1..].iter().enumerate() {
            writer.write_field(reader.header().rid2name(rec.rid().unwrap())?)?;
            writer.write_integer(rec.pos() as i32 + 1)?;
            writer.write_field(&alleles[0])?;
            writer.write_field(allele)?;
            match rec.qual() {
                q if q.is_missing() => writer.write_field(b"")?,
                q => writer.write_float(q)?,
            }

            if show_filter {
                if rec.has_filter(".".as_bytes()) {
                    writer.write_field(b"")?
                } else if rec.has_filter("PASS".as_bytes()) {
                    writer.write_field(b"PASS")?
                } else {
                    let mut filters = Vec::new();
                    for (i, filter) in rec.filters().enumerate() {
                        if i != 0 {
                            filters.push(b';');
                        }
                        filters.extend_from_slice(&reader.header().id_to_name(filter));
                    }
                    writer.write_field(&filters)?;
                }
            }

            for name in info_tags {
                let _name = name.as_bytes();
                if let Ok((tag_type, tag_length)) = rec.header().info_type(_name) {
                    let get_idx = || match tag_length {
                        bcf::header::TagLength::Fixed(_) => Ok(0),
                        bcf::header::TagLength::AltAlleles => Ok(i),
                        bcf::header::TagLength::Alleles => Ok(i + 1),
                        bcf::header::TagLength::Variable => Ok(0),
                        _ => Err(Box::new(ParseError::UnsupportedTagLength)),
                    };

                    match tag_type {
                        bcf::header::TagType::Flag => {
                            writer.write_flag(rec.info(_name).flag()?)?;
                        }
                        bcf::header::TagType::Integer => {
                            let i = get_idx()?;
                            if let Some(values) = rec.info(_name).integer()? {
                                writer.write_integer(values[i])?;
                            } else {
                                writer.write_field(b"")?;
                            }
                        }
                        bcf::header::TagType::Float => {
                            let i = get_idx()?;
                            if let Some(values) = rec.info(_name).float()? {
                                writer.write_float(values[i])?;
                            } else {
                                writer.write_field(b"")?;
                            }
                        }
                        bcf::header::TagType::String => {
                            let i = get_idx()?;
                            if let Some(values) = rec.info(_name).string()? {
                                writer.write_field(values[i])?;
                            } else {
                                writer.write_field(b"")?;
                            }
                        }
                    }
                } else {
                    // tag undefined, write NA
                    writer.write_field(b"")?;
                }
            }

            let genotypes = if show_genotypes {
                let genotypes = rec.genotypes()?;

                Some(
                    (0..reader.header().sample_count() as usize)
                        .map(|s| format!("{}", genotypes.get(s)))
                        .collect_vec(),
                )
            } else {
                None
            };

            for s in 0..reader.header().sample_count() as usize {
                if let Some(ref genotypes) = genotypes {
                    writer.write_field(genotypes[s].as_bytes())?;
                }
                for name in format_tags {
                    let _name = name.as_bytes();
                    if let Ok((tag_type, tag_length)) = reader.header().format_type(_name) {
                        let i = match tag_length {
                            bcf::header::TagLength::Fixed(_) => 0,
                            bcf::header::TagLength::AltAlleles => i,
                            bcf::header::TagLength::Alleles => i + 1,
                            _ => bail!(ParseError::UnsupportedTagLength),
                        };

                        match tag_type {
                            bcf::header::TagType::Flag => {
                                panic!("Unable to find FORMAT \"{0}\" in the input file! Is \"{0}\" an INFO tag?", name);
                            }
                            bcf::header::TagType::Integer => {
                                writer.write_field(
                                    format!("{}", rec.format(_name).integer()?[s][i]).as_bytes(),
                                )?;
                            }
                            bcf::header::TagType::Float => {
                                writer.write_field(
                                    format!("{}", rec.format(_name).float()?[s][i]).as_bytes(),
                                )?;
                            }
                            bcf::header::TagType::String => {
                                writer.write_field(rec.format(_name).string()?[s])?;
                            }
                        }
                    } else {
                        // tag undefined, write NA
                        writer.write_field(b"")?;
                    }
                }
            }
            writer.newline()?;
        }
    }

    Ok(())
}

#[derive(Error, Debug)]
pub enum ParseError {
    #[error("currently, only R, A, and 1 are supported multiplicities of tags")]
    UnsupportedTagLength,
}
