use std::error::Error;
use std::io;
use std::str;
use std::io::Write;

use itertools::Itertools;
use rust_htslib::bcf;
use rust_htslib::bcf::record::Numeric;


pub struct Writer {
    inner: io::BufWriter<io::Stdout>,
    field_count: usize
}


impl Writer {
    fn new(inner: io::BufWriter<io::Stdout>) -> Self {
        Writer {
            inner: inner,
            field_count: 0
        }
    }

    fn write_integer(&mut self, value: i32) -> Result<(), Box<Error>> {
        self.write_field(format!("{}", value).as_bytes())
    }

    fn write_float(&mut self, value: f32) -> Result<(), Box<Error>> {
        self.write_field(format!("{}", value).as_bytes())
    }

    fn write_flag(&mut self, value: bool) -> Result<(), Box<Error>> {
        self.write_field(format!("{}", value).as_bytes())
    }

    fn write_field(&mut self, value: &[u8]) -> Result<(), Box<Error>> {
        if self.field_count > 0 {
            try!(self.inner.write(b"\t"));
        }
        try!(self.inner.write(value));
        self.field_count += 1;
        Ok(())
    }

    fn newline(&mut self) -> Result<(), Box<Error>> {
        try!(self.inner.write(b"\n"));
        self.field_count = 0;
        Ok(())
    }
}


const HEADER_COMMON: &'static [u8] = b"VARIANT";


pub fn to_txt(
    info_tags: &[&str],
    format_tags: &[&str],
    show_genotypes: bool
) -> Result<(), Box<Error>> {
    let reader = try!(bcf::Reader::from_path(&"-"));
    let mut writer = Writer::new(io::BufWriter::new(io::stdout()));

    let common_n = 5 + info_tags.len();
    try!(writer.write_field(HEADER_COMMON));
    for _ in 1..common_n {
        try!(writer.write_field(HEADER_COMMON));
    }
    for sample in reader.header.samples() {
        try!(writer.write_field(sample));
        for _ in 1..format_tags.len() + show_genotypes as usize {
            try!(writer.write_field(sample));
        }
    }
    try!(writer.newline());
    try!(writer.write_field(b"CHROM"));
    try!(writer.write_field(b"POS"));
    try!(writer.write_field(b"REF"));
    try!(writer.write_field(b"ALT"));
    try!(writer.write_field(b"QUAL"));
    for name in info_tags {
        try!(writer.write_field(name.as_bytes()));
    }
    for _ in 0..reader.header.sample_count() {
        if show_genotypes {
            try!(writer.write_field(b"GT"));
        }
        for name in format_tags {
            try!(writer.write_field(name.as_bytes()));
        }
    }
    try!(writer.newline());

    let mut record = bcf::Record::new();
    loop {
        if let Err(e) = reader.read(&mut record) {
            if e.is_eof() {
                break;
            } else {
                return Err(Box::new(e));
            }
        }

        let alleles = record.alleles().into_iter().map(|a| a.to_owned()).collect_vec();
        for (i, allele) in alleles[1..].iter().enumerate() {
            try!(writer.write_field(reader.header.rid2name(record.rid().unwrap())));
            try!(writer.write_integer(record.pos() as i32 + 1));
            try!(writer.write_field(&alleles[0]));
            try!(writer.write_field(allele));
            match record.qual() {
                q if q.is_missing() => try!(writer.write_field(b"")),
                q => try!(writer.write_float(q))
            }

            for name in info_tags {
                let _name = name.as_bytes();
                if let Ok((tag_type, tag_length)) = reader.header.info_type(_name) {
                    let get_idx = || {
                        match tag_length {
                            bcf::header::TagLength::Fixed => {
                                Ok(0)
                            },
                            bcf::header::TagLength::AltAlleles => {
                                Ok(i)
                            },
                            bcf::header::TagLength::Variable => {
                                Ok(i)
                            },
                            _ => Err(Box::new(ParseError::UnsupportedTagLength))
                        }
                    };

                    match tag_type {
                        bcf::header::TagType::Flag => {
                            try!(writer.write_flag(try!(record.info(_name).flag())));
                        },
                        bcf::header::TagType::Integer => {
                            let i = try!(get_idx());
                            if let Some(values) = try!(record.info(_name).integer()) {
                                try!(writer.write_integer(values[i]));
                            } else {
                                try!(writer.write_field(b""));
                            }
                        },
                        bcf::header::TagType::Float => {
                            let i = try!(get_idx());
                            if let Some(values) = try!(record.info(_name).float()) {
                                try!(writer.write_float(values[i]));
                            } else {
                                try!(writer.write_field(b""));
                            }
                        },
                        bcf::header::TagType::String => {
                            let i = try!(get_idx());
                            if let Some(values) = try!(record.info(_name).string()) {
                                try!(writer.write_field(values[i]));
                            } else {
                                try!(writer.write_field(b""));
                            }
                        }
                    }
                } else {
                    // tag undefined, write NA
                    try!(writer.write_field(b""));
                }
            }

            let genotypes = if show_genotypes {
                let genotypes = try!(record.genotypes());
                Some(
                    (0..reader.header.sample_count() as usize).map(|s| {
                        format!("{}", genotypes.get(s))
                    }).collect_vec()
                )
            } else {
                None
            };

            for s in 0..reader.header.sample_count() as usize {
                if let Some(ref genotypes) = genotypes {
                    try!(writer.write_field(genotypes[s].as_bytes()));
                }
                for name in format_tags {
                    let _name = name.as_bytes();
                    if let Ok((tag_type, tag_length)) = reader.header.format_type(_name) {
                        let i = match tag_length {
                            bcf::header::TagLength::Fixed => {
                                0
                            },
                            bcf::header::TagLength::AltAlleles => {
                                i
                            },
                            _ => return Err(Box::new(ParseError::UnsupportedTagLength))
                        };

                        match tag_type {
                            bcf::header::TagType::Flag => {
                                panic!("there is no flag type for format");
                            },
                            bcf::header::TagType::Integer => {
                                try!(writer.write_field(format!("{}", try!(record.format(_name).integer())[s][i]).as_bytes()));
                            },
                            bcf::header::TagType::Float => {
                                try!(writer.write_field(format!("{}", try!(record.format(_name).float())[s][i]).as_bytes()));
                            },
                            bcf::header::TagType::String => {
                                try!(writer.write_field(try!(record.format(_name).string())[s]));
                            }
                        }
                    } else {
                        // tag undefined, write NA
                        try!(writer.write_field(b""));
                    }
                }
            }
            try!(writer.newline());
        }
    }

    Ok(())
}


quick_error! {
    #[derive(Debug)]
    pub enum ParseError {
        UnsupportedTagLength {
            description("currently, only A and 1 are supported multiplicities of tags")
        }

    }
}
