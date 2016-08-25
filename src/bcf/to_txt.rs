use std::error::Error;
use std::io;
use std::str;
use std::io::Write;

use itertools::Itertools;
use rust_htslib::bcf;


pub struct Writer(io::BufWriter<io::Stdout>);


impl Writer {
    fn write_field(&mut self, value: &[u8]) -> Result<(), Box<Error>> {
        let &mut Writer(ref mut w) = self;
        try!(w.write(value));
        try!(w.write(b"\t"));

        Ok(())
    }
    fn newline(&mut self) -> Result<(), Box<Error>> {
        let &mut Writer(ref mut w) = self;
        try!(w.write(b"\n"));

        Ok(())
    }
}


pub fn to_txt(
    info_tags: &[&str],
    format_tags: &[&str],
) -> Result<(), Box<Error>> {
    let reader = try!(bcf::Reader::new(&"-"));
    let mut writer = Writer(io::BufWriter::new(io::stdout()));

    let common_n = 5 + info_tags.len();
    try!(writer.write_field(b"COMMON"));
    for _ in 1..common_n {
        try!(writer.write_field(b""));
    }
    for sample in reader.header.samples() {
        try!(writer.write_field(sample));
        for _ in 1..format_tags.len() {
            try!(writer.write_field(b""));
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
            try!(writer.write_field(format!("{}", record.pos()).as_bytes()));
            try!(writer.write_field(&alleles[0]));
            try!(writer.write_field(allele));
            try!(writer.write_field(format!("{}", record.qual()).as_bytes()));

            for name in info_tags {
                let _name = name.as_bytes();
                let (tag_type, tag_length) = try!(reader.header.info_type(_name));
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
                        try!(writer.write_field(format!("{}", try!(record.info(_name).flag())).as_bytes()));
                    },
                    bcf::header::TagType::Integer => {
                        try!(writer.write_field(format!("{}", try!(record.info(_name).integer())[i]).as_bytes()));
                    },
                    bcf::header::TagType::Float => {
                        try!(writer.write_field(format!("{}", try!(record.info(_name).float())[i]).as_bytes()));
                    },
                    bcf::header::TagType::String => {
                        try!(writer.write_field(try!(record.info(_name).string())[i]));
                    }
                }
            }

            for s in 0..reader.header.sample_count() as usize {
                for name in format_tags {
                    let _name = name.as_bytes();
                    let (tag_type, tag_length) = try!(reader.header.format_type(_name));
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
