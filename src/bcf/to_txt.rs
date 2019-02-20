use std::error::Error;
use std::io;
use std::io::Write;
use std::str;

use itertools::Itertools;
use rust_htslib::bcf;
use rust_htslib::bcf::record::Numeric;

pub struct Writer {
    inner: io::BufWriter<io::Stdout>,
    field_count: usize,
}

impl Writer {
    fn new(inner: io::BufWriter<io::Stdout>) -> Self {
        Writer {
            inner: inner,
            field_count: 0,
        }
    }

    fn write_integer(&mut self, value: i32) -> Result<(), Box<dyn Error>> {
        let fmt = if value.is_missing() {
            "".to_owned()
        } else {
            format!("{}", value)
        };
        self.write_field(fmt.as_bytes())
    }

    fn write_float(&mut self, value: f32) -> Result<(), Box<dyn Error>> {
        let fmt = if value.is_missing() {
            "".to_owned()
        } else {
            format!("{}", value)
        };
        self.write_field(fmt.as_bytes())
    }

    fn write_flag(&mut self, value: bool) -> Result<(), Box<dyn Error>> {
        self.write_field(format!("{}", value).as_bytes())
    }

    fn write_field(&mut self, value: &[u8]) -> Result<(), Box<dyn Error>> {
        if self.field_count > 0 {
            r#try!(self.inner.write(b"\t"));
        }
        r#try!(self.inner.write(value));
        self.field_count += 1;
        Ok(())
    }

    fn newline(&mut self) -> Result<(), Box<dyn Error>> {
        r#try!(self.inner.write(b"\n"));
        self.field_count = 0;
        Ok(())
    }
}

const HEADER_COMMON: &'static [u8] = b"VARIANT";

pub fn to_txt(
    info_tags: &[&str],
    format_tags: &[&str],
    show_genotypes: bool,
) -> Result<(), Box<dyn Error>> {
    let mut reader = bcf::Reader::from_stdin()?;
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

    let mut record = bcf::Record::new();
    loop {
        if let Err(e) = reader.read(&mut record) {
            if e.is_eof() {
                break;
            } else {
                return Err(Box::new(e));
            }
        }

        let alleles = record
            .alleles()
            .into_iter()
            .map(|a| a.to_owned())
            .collect_vec();
        for (i, allele) in alleles[1..].iter().enumerate() {
            r#try!(writer.write_field(reader.header().rid2name(record.rid().unwrap())));
            r#try!(writer.write_integer(record.pos() as i32 + 1));
            r#try!(writer.write_field(&alleles[0]));
            r#try!(writer.write_field(allele));
            match record.qual() {
                q if q.is_missing() => r#try!(writer.write_field(b"")),
                q => r#try!(writer.write_float(q)),
            }

            for name in info_tags {
                let _name = name.as_bytes();
                if let Ok((tag_type, tag_length)) = reader.header().info_type(_name) {
                    let get_idx = || match tag_length {
                        bcf::header::TagLength::Fixed => Ok(0),
                        bcf::header::TagLength::AltAlleles => Ok(i),
                        bcf::header::TagLength::Alleles => Ok(i + 1),
                        bcf::header::TagLength::Variable => Ok(0),
                        _ => Err(Box::new(ParseError::UnsupportedTagLength)),
                    };

                    match tag_type {
                        bcf::header::TagType::Flag => {
                            r#try!(writer.write_flag(r#try!(record.info(_name).flag())));
                        }
                        bcf::header::TagType::Integer => {
                            let i = r#try!(get_idx());
                            if let Some(values) = r#try!(record.info(_name).integer()) {
                                r#try!(writer.write_integer(values[i]));
                            } else {
                                r#try!(writer.write_field(b""));
                            }
                        }
                        bcf::header::TagType::Float => {
                            let i = r#try!(get_idx());
                            if let Some(values) = r#try!(record.info(_name).float()) {
                                r#try!(writer.write_float(values[i]));
                            } else {
                                r#try!(writer.write_field(b""));
                            }
                        }
                        bcf::header::TagType::String => {
                            let i = r#try!(get_idx());
                            if let Some(values) = r#try!(record.info(_name).string()) {
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
                let genotypes = r#try!(record.genotypes());
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
                            bcf::header::TagLength::Fixed => 0,
                            bcf::header::TagLength::AltAlleles => i,
                            bcf::header::TagLength::Alleles => i + 1,
                            _ => return Err(Box::new(ParseError::UnsupportedTagLength)),
                        };

                        match tag_type {
                            bcf::header::TagType::Flag => {
                                panic!("there is no flag type for format");
                            }
                            bcf::header::TagType::Integer => {
                                r#try!(
                                    writer.write_field(
                                        format!("{}", r#try!(record.format(_name).integer())[s][i])
                                            .as_bytes()
                                    )
                                );
                            }
                            bcf::header::TagType::Float => {
                                r#try!(
                                    writer.write_field(
                                        format!("{}", r#try!(record.format(_name).float())[s][i])
                                            .as_bytes()
                                    )
                                );
                            }
                            bcf::header::TagType::String => {
                                r#try!(writer.write_field(r#try!(record.format(_name).string())[s]));
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
    }

    Ok(())
}

quick_error! {
    #[derive(Debug)]
    pub enum ParseError {
        UnsupportedTagLength {
            description("currently, only R, A, and 1 are supported multiplicities of tags")
        }

    }
}
