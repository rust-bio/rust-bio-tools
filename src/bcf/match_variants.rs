use std::error::Error;
use std::collections::{VecDeque, vec_deque};
use std::str;
use std::cmp;
use std::u32;
use std::mem;

use itertools::Itertools;
use rust_htslib::bcf;


pub fn match_variants(matchbcf: &str, max_dist: u32, max_len_diff: u32) -> Result<(), Box<Error>> {
    let inbcf = try!(bcf::Reader::new(&"-"));
    let mut header = bcf::Header::with_template(&inbcf.header);

    header.push_record(
        format!("##INFO=<ID=MATCHING,Number=A,Type=Integer,\
        Description=\"For each alternative allele, 1 if it matches a variant in another VCF/BCF \
        and 0 otherwise. For indels, matching is fuzzy: distance of centres <= {}, difference of \
        lengths <= {}\">", max_dist, max_len_diff).as_bytes()
    );
    header.push_record(
        format!("##rust-bio-tools={}", env!("CARGO_PKG_VERSION")).as_bytes()
    );
    header.push_record(
        b"##rust-bio-tools-subcommand=vcf-match"
    );
    let mut outbcf = try!(bcf::Writer::new(&"-", &header, false, false));
    let mut buffer = RecordBuffer::new(try!(bcf::Reader::new(&matchbcf)), max_dist);

    let mut rec = bcf::Record::new();
    loop {
        if let Err(e) = inbcf.read(&mut rec) {
            if e.is_eof() {
                break;
            }
            return Err(Box::new(e));
        }
        outbcf.translate(&mut rec);

        if let Some(rid) = rec.rid() {
            let chrom = inbcf.header.rid2name(rid);
            let pos = rec.pos();
            // move buffer to pos
            try!(buffer.fill(chrom, pos));

            let var = try!(Variant::new(&mut rec));
            let matching = var.alleles.iter().map(|a| {
                for v in buffer.iter() {
                    if var.matches(v, a, max_dist, max_len_diff) {
                        return 1;
                    }
                }
                0
            }).collect_vec();

            try!(rec.push_info_integer(b"MATCHING", &matching));
            try!(outbcf.write(&rec));
        } else {
            try!(outbcf.write(&rec));
        }
    }

    Ok(())
}


#[derive(Debug)]
pub struct Variant {
    rid: u32,
    pos: u32,
    alleles: Vec<VariantType>
}


impl Variant {

    pub fn new(rec: &mut bcf::Record) -> Result<Self, Box<Error>> {
        let pos = rec.pos();

        let svlen = if let Ok(Some(svlen)) = rec.info(b"SVLEN").integer() {
            Some(svlen.to_owned())
        } else { None };
        let svtype = if let Ok(Some(svtype)) = rec.info(b"SVTYPE").string() {
            Some(svtype.into_iter().map(|s| s.to_owned()).collect_vec())
        } else { None };

        let mut _alleles = Vec::new();
        let alleles = rec.alleles();
        let refallele = alleles[0];
        for (i, &a) in alleles[1..].iter().enumerate() {
            let vartype = match (&svlen, &svtype) {

                (&Some(ref svlen), &Some(ref svtype)) => {
                    let svlen = svlen.get(i).unwrap_or(&svlen[0]).abs();
                    let svtype = svtype.get(i).unwrap_or(&svtype[0]);
                    match &svtype[..] {
                        b"INS" => VariantType::Insertion(svlen as u32),
                        b"DEL" => VariantType::Deletion(svlen as u32),
                        t => return Err(Box::new(MatchError::UnsupportedVariant(try!(str::from_utf8(t)).to_owned())))
                    }
                },

                (&Some(ref svlen), _) if a[0] == b'<' => {
                    let svlen = svlen.get(i).unwrap_or(&svlen[0]).abs();
                    match a {
                        b"<DEL>" => VariantType::Deletion(svlen as u32),
                        b"<INS>" => VariantType::Insertion(svlen as u32),
                        a => return Err(Box::new(MatchError::UnsupportedVariant(try!(str::from_utf8(a)).to_owned())))
                    }
                },

                _ => {
                    if a.len() < refallele.len() {
                        VariantType::Deletion((refallele.len() - a.len()) as u32)
                    } else if a.len() > refallele.len() {
                        VariantType::Insertion((a.len() - refallele.len()) as u32)
                    } else if a.len() == 1 {
                        VariantType::SNV(a[0])
                    } else {
                        return Err(Box::new(MatchError::UnsupportedVariant("complex".to_owned())))
                    }
                }
            };
            _alleles.push(vartype);
        }
        Ok(Variant {
            rid: rec.rid().unwrap(),
            pos: pos,
            alleles: _alleles
        })
    }

    pub fn centerpoint(&self, allele: &VariantType) -> u32 {
        match allele {
            &VariantType::SNV(_) => self.pos,
            &VariantType::Insertion(_) => self.pos,
            &VariantType::Deletion(len) => (self.pos as f64 + len as f64 / 2.0) as u32
        }
    }

    pub fn matches(&self, other: &Variant, allele: &VariantType, max_dist: u32, max_len_diff: u32) -> bool {
        for b in &other.alleles {
            let dist = (self.centerpoint(allele) as i32 - other.centerpoint(b) as i32).abs() as u32;
            match (allele, b) {
                (&VariantType::SNV(a), &VariantType::SNV(b)) => {
                    if a == b && dist == 0 {
                        return true;
                    }
                },
                (&VariantType::Insertion(l1), &VariantType::Insertion(l2)) |
                (&VariantType::Deletion(l1), &VariantType::Deletion(l2)) => {
                    if (l1 as i32 - l2 as i32).abs() as u32 <= max_len_diff && dist <= max_dist {
                        return true;
                    }
                },
                _ => continue
            }
        }
        false
    }
}


#[derive(Debug)]
pub enum VariantType {
    SNV(u8),
    Insertion(u32),
    Deletion(u32)
}


quick_error! {
    #[derive(Debug)]
    pub enum MatchError {
        UnsupportedVariant(vartype: String) {
            description("unsupported variant")
            display("variant type {} is not supported", vartype)
        }
    }
}


pub struct RecordBuffer {
    reader: bcf::Reader,
    ringbuffer: VecDeque<Variant>,
    ringbuffer2: VecDeque<Variant>,
    window: u32
}

impl RecordBuffer {
    pub fn new(reader: bcf::Reader, window: u32) -> Self {
        RecordBuffer {
            reader: reader,
            ringbuffer: VecDeque::new(),
            ringbuffer2: VecDeque::new(),
            window: window
        }
    }

    pub fn rid(&self) -> Option<u32> {
        self.ringbuffer.back().map(|var| var.rid)
    }

    pub fn fill(&mut self, chrom: &[u8], pos: u32) -> Result<(), Box<Error>> {
        let window_start = cmp::max(pos as i32 - self.window as i32 - 1, 0) as u32;
        let window_end = pos + self.window;
        let rid = try!(self.reader.header.name2rid(chrom));

        if let Some(last_rid) = self.rid() {
            if rid != last_rid {
                // swap with buffer for next rid
                mem::swap(&mut self.ringbuffer2, &mut self.ringbuffer);
                // clear second buffer
                self.ringbuffer2.clear();
            } else {
                // remove records too far left of from wrong rid
                let to_remove = self.ringbuffer.iter().take_while(|var| var.pos < window_start || var.rid != rid).count();
                self.ringbuffer.drain(..to_remove);
            }
        }

        // extend to the right
        let mut rec = bcf::Record::new();
        loop {
            if let Err(e) = self.reader.read(&mut rec) {
                if e.is_eof() {
                    break;
                }
                return Err(Box::new(e));
            }
            let pos = rec.pos();
            if let Some(rec_rid) = rec.rid() {
                if rec_rid == rid {
                    if pos > window_end {
                        // Record is beyond our window. Store it anyways but stop.
                        self.ringbuffer.push_back(try!(Variant::new(&mut rec)));
                        break;
                    } else if pos >= window_start {
                        // Record is within out window.
                        self.ringbuffer.push_back(try!(Variant::new(&mut rec)));
                    }
                } else if rec_rid > rid {
                    // record comes from next rid. Store it in second buffer but stop filling.
                    self.ringbuffer2.push_back(try!(Variant::new(&mut rec)));
                    break;
                } else {
                    // Record comes from previous rid. Ignore it.
                    continue;
                }
            } else {
                // skip records without proper rid
                continue;
            }
        }

        Ok(())
    }

    pub fn iter(&self) -> vec_deque::Iter<Variant> {
        self.ringbuffer.iter()
    }
}
