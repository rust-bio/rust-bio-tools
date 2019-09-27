//! Create a variant table from a VCF file.
//!
//! ## Usage:
//! ```bash
//! $ rbt vcf-to-txt --genotypes --fmt S --info T X SOMATIC < tests/test.vcf > tests/variant-table.txt
//! ```
//!
use crate::errors;
use derive_new::new;
use itertools::Itertools;
use rust_htslib::bcf;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::Read;
use snafu::ResultExt;
use std::io;
use std::io::Write;
use std::str;

#[derive(new)]
pub struct Writer {
    inner: io::BufWriter<io::Stdout>,
    #[new(value = "0")]
    field_count: usize,
}

impl Writer {
    fn write_integer(&mut self, value: i32) -> errors::Result<()> {
        let fmt = if value.is_missing() {
            "".to_owned()
        } else {
            format!("{}", value)
        };
        self.write_field(fmt.as_bytes())
    }

    fn write_float(&mut self, value: f32) -> errors::Result<()> {
        let fmt = if value.is_missing() {
            "".to_owned()
        } else {
            format!("{}", value)
        };
        self.write_field(fmt.as_bytes())
    }

    fn write_flag(&mut self, value: bool) -> errors::Result<()> {
        self.write_field(format!("{}", value).as_bytes())
    }

    fn write_field(&mut self, value: &[u8]) -> errors::Result<()> {
        if self.field_count > 0 {
            self.inner
                .write(b"\t")
                .context(errors::BCFStdIOWriteError)?;
        }
        self.inner
            .write(value)
            .context(errors::BCFStdIOWriteError)?;
        self.field_count += 1;
        Ok(())
    }

    fn newline(&mut self) -> errors::Result<()> {
        self.inner
            .write(b"\n")
            .context(errors::BCFStdIOWriteError)?;
        self.field_count = 0;
        Ok(())
    }
}

const HEADER_COMMON: &'static [u8] = b"VARIANT";

pub fn to_txt(
    info_tags: &[&str],
    format_tags: &[&str],
    show_genotypes: bool,
) -> errors::Result<()> {
    let mut reader = bcf::Reader::from_stdin().context(errors::BCFReaderStdinError)?;
    let mut writer = Writer::new(io::BufWriter::new(io::stdout()));

    let common_n = 5 + info_tags.len();
    r#try!(writer.write_field(HEADER_COMMON));
    for _ in 1..common_n {
        r#try!(writer.write_field(HEADER_COMMON));
    }
    let show_samples = show_genotypes || !format_tags.is_empty();
    if show_samples {
        for sample in reader.header().samples() {
            r#try!(writer.write_field(sample));
            for _ in 1..format_tags.len() + show_genotypes as usize {
                r#try!(writer.write_field(sample));
            }
        }
    }
    r#try!(writer.newline());
    r#try!(writer.write_field(b"CHROM"));
    r#try!(writer.write_field(b"POS"));
    r#try!(writer.write_field(b"REF"));
    r#try!(writer.write_field(b"ALT"));
    r#try!(writer.write_field(b"QUAL"));
    for name in info_tags {
        r#try!(writer.write_field(name.as_bytes()));
    }
    if show_samples {
        for _ in 0..reader.header().sample_count() {
            if show_genotypes {
                r#try!(writer.write_field(b"GT"));
            }
            for name in format_tags {
                r#try!(writer.write_field(name.as_bytes()));
            }
        }
    }
    r#try!(writer.newline());
    let mut rec = reader.empty_record();
    let mut record_idx: usize = 0;
    loop {
        match reader.read(&mut rec) {
            Err(e) => {
                return Err(e).context(errors::BCFReadError {
                    header: format!("{:?}", reader.header()),
                });
            }
            Ok(false) => {
                // reached end of file
                break;
            }
            Ok(true) => {

        let alleles = rec
            .alleles()
            .into_iter()
            .map(|a| a.to_owned())
            .collect_vec();
        for (i, allele) in alleles[1..].iter().enumerate() {
            r#try!(
                writer.write_field(reader.header().rid2name(rec.rid().unwrap()).context(
                    errors::BCFReadIdError {
                        rid: rec.rid().unwrap(),
                        header: format!("{:?}", reader.header()),
                    }
                )?)
            );
            r#try!(writer.write_integer(rec.pos() as i32 + 1));
            r#try!(writer.write_field(&alleles[0]));
            r#try!(writer.write_field(allele));
            match rec.qual() {
                q if q.is_missing() => r#try!(writer.write_field(b"")),
                q => r#try!(writer.write_float(q)),
            }

            for name in info_tags {
                let _name = name.as_bytes();
                if let Ok((tag_type, tag_length)) = rec.header().info_type(_name) {
                    let get_idx = || match tag_length {
                        bcf::header::TagLength::Fixed(_) => Ok(0),
                        bcf::header::TagLength::AltAlleles => Ok(i),
                        bcf::header::TagLength::Alleles => Ok(i + 1),
                        bcf::header::TagLength::Variable => Ok(0),
                        _ => Err(errors::Error::UnsupportedTagLengthError {
                            tag_length: format!("{:?}", tag_length),
                        }),
                    };

                    match tag_type {
                        bcf::header::TagType::Flag => {
                            r#try!(writer.write_flag(
                                rec.info(_name)
                                    .flag()
                                    .context(errors::BCFInfoReadError { record_idx: i })?
                            ));
                        }
                        bcf::header::TagType::Integer => {
                            let i = r#try!(get_idx());
                            if let Some(values) = rec
                                .info(_name)
                                .integer()
                                .context(errors::BCFInfoReadError { record_idx })?
                            {
                                r#try!(writer.write_integer(values[i]));
                            } else {
                                r#try!(writer.write_field(b""));
                            }
                        }
                        bcf::header::TagType::Float => {
                            let i = r#try!(get_idx());
                            if let Some(values) = rec
                                .info(_name)
                                .float()
                                .context(errors::BCFInfoReadError { record_idx })?
                            {
                                r#try!(writer.write_float(values[i]));
                            } else {
                                r#try!(writer.write_field(b""));
                            }
                        }
                        bcf::header::TagType::String => {
                            let i = r#try!(get_idx());
                            if let Some(values) = rec
                                .info(_name)
                                .string()
                                .context(errors::BCFInfoReadError { record_idx })?
                            {
                                r#try!(writer.write_field(values[i]));
                            } else {
                                r#try!(writer.write_field(b""));
                            }
                        }
                    }
                } else {
                    // tag undefined, write NA
                    r#try!(writer.write_field(b""));
                }
            }

            let genotypes = if show_genotypes {
                let genotypes = rec.genotypes().context(errors::BCFFormatReadError)?;
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
                    r#try!(writer.write_field(genotypes[s].as_bytes()));
                }
                for name in format_tags {
                    let _name = name.as_bytes();
                    if let Ok((tag_type, tag_length)) = reader.header().format_type(_name) {
                        let i = match tag_length {
                            bcf::header::TagLength::Fixed(_) => 0,
                            bcf::header::TagLength::AltAlleles => i,
                            bcf::header::TagLength::Alleles => i + 1,
                            _ => Err(errors::Error::UnsupportedTagLengthError {
                                tag_length: format!("{:?}", tag_length),
                            })?,
                        };

                        match tag_type {
                            bcf::header::TagType::Flag => Err(errors::Error::FlagTypeError)?,
                            bcf::header::TagType::Integer => {
                                r#try!(writer.write_field(
                                    format!(
                                        "{}",
                                        rec.format(_name)
                                            .integer()
                                            .context(errors::BCFFormatReadError)?[s][i]
                                    )
                                    .as_bytes()
                                ));
                            }
                            bcf::header::TagType::Float => {
                                r#try!(writer.write_field(
                                    format!(
                                        "{}",
                                        rec.format(_name)
                                            .float()
                                            .context(errors::BCFFormatReadError)?[s][i]
                                    )
                                    .as_bytes()
                                ));
                            }
                            bcf::header::TagType::String => {
                                r#try!(writer.write_field(
                                    rec.format(_name)
                                        .string()
                                        .context(errors::BCFFormatReadError)?[s]
                                ));
                            }
                        }
                    } else {
                        // tag undefined, write NA
                        r#try!(writer.write_field(b""));
                    }
                }
            }
            r#try!(writer.newline());
        }
        record_idx += 1;
            }
        }
    }
    Ok(())
}
