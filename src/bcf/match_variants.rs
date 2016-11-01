use std::error::Error;
use std::collections::{VecDeque, vec_deque};
use std::str;
use std::cmp;
use std::mem;

use itertools::Itertools;
use rust_htslib::bcf;


pub fn match_variants(matchbcf: &str, max_dist: u32, max_len_diff: u32) -> Result<(), Box<Error>> {
    let inbcf = try!(bcf::Reader::new(&"-"));
    let mut header = bcf::Header::with_template(&inbcf.header);

    header.push_record(
        format!("##INFO=<ID=MATCHING,Number=A,Type=Integer,\
        Description=\"For each alternative allele, -1 if it does not match a variant in another VCF/BCF. \
        If it matches a variant, an id i>=0 is points to the i-th variant in the VCF/BCF (counting each \
        alternative allele separately). For indels, matching is fuzzy: distance of centres <= {}, difference of \
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
    let mut i = 0;
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

            let var = try!(Variant::new(&mut rec, &mut i));
            let matching = var.alleles.iter().map(|a| {
                for v in buffer.iter() {
                    if let Some(id) = var.matches(v, a, max_dist, max_len_diff) {
                        return id as i32;
                    }
                }
                -1
            }).collect_vec();

            try!(rec.push_info_integer(b"MATCHING", &matching));
            try!(outbcf.write(&rec));
        } else {
            try!(outbcf.write(&rec));
        }

        if (i) % 1000 == 0 {
            info!("{} variants written.", i);
        }
    }
    info!("{} variants written.", i);

    Ok(())
}


#[derive(Debug)]
pub struct Variant {
    id: u32,
    rid: u32,
    pos: u32,
    alleles: Vec<VariantType>
}


impl Variant {

    pub fn new(rec: &mut bcf::Record, id: &mut u32) -> Result<Self, Box<Error>> {
        let pos = rec.pos();

        let svlen = if let Ok(Some(svlen)) = rec.info(b"SVLEN").integer() {
            Some(svlen[0].abs() as u32)
        } else { None };
        let svtype = if let Ok(Some(svtype)) = rec.info(b"SVTYPE").string() {
            Some(svtype[0].to_owned())
        } else { None };
        let end = if let Ok(Some(end)) = rec.info(b"END").integer() {
            Some(end[0] as u32)
        } else { None };
        let inslen = if let Ok(Some(inslen)) = rec.info(b"INSLEN").integer() {
            Some(inslen[0] as u32)
        } else { None };
        let alleles = rec.alleles();
        let refallele = alleles[0];

        let _alleles: Vec<VariantType> = if let Some(svtype) = svtype {
            vec![
                if svtype == b"INS" {
                    let svlen = match (svlen, inslen) {
                        (Some(svlen), _)     => svlen,
                        (None, Some(inslen)) => inslen,
                        _ => {
                            return Err(Box::new(MatchError::MissingTag("SVLEN or INSLEN".to_owned())));
                        }
                    };
                    VariantType::Insertion(svlen)
                } else if svtype == b"DEL" {
                    let svlen = match(svlen, end) {
                        (Some(svlen), _)  => svlen,
                        (None, Some(end)) => end - 1 - pos,
                        _ => {
                            return Err(Box::new(MatchError::MissingTag("SVLEN or END".to_owned())));
                        }
                    };
                    VariantType::Deletion(svlen)
                } else {
                    warn!("Unsupported variant {}", try!(str::from_utf8(&svtype)));
                    VariantType::Unsupported
                }
            ]
        } else {
            let mut _alleles = Vec::with_capacity(alleles.len() - 1);
            for a in &alleles[1..] {
                _alleles.push(
                    if a.len() < refallele.len() {
                        VariantType::Deletion((refallele.len() - a.len()) as u32)
                    } else if a.len() > refallele.len() {
                        VariantType::Insertion((a.len() - refallele.len()) as u32)
                    } else if a.len() == 1 {
                        VariantType::SNV(a[0])
                    } else {
                        warn!("Unsupported variant {} -> {}", try!(str::from_utf8(refallele)), try!(str::from_utf8(a)));
                        VariantType::Unsupported
                    }
                );
            }
            _alleles
        };

        let var = Variant {
            id: *id,
            rid: rec.rid().unwrap(),
            pos: pos,
            alleles: _alleles
        };
        *id += alleles.len() as u32 - 1;
        Ok(var)
    }

    pub fn centerpoint(&self, allele: &VariantType) -> u32 {
        match allele {
            &VariantType::SNV(_) => self.pos,
            &VariantType::Insertion(_) => self.pos,
            &VariantType::Deletion(len) => (self.pos as f64 + len as f64 / 2.0) as u32,
            &VariantType::Unsupported => panic!("Unsupported variant.")
        }
    }

    pub fn matches(&self, other: &Variant, allele: &VariantType, max_dist: u32, max_len_diff: u32) -> Option<u32> {
        if allele.is_unsupported() {
            return None;
        }
        for (j, b) in other.alleles.iter().enumerate() {
            if b.is_unsupported() {
                continue;
            }
            let dist = (self.centerpoint(allele) as i32 - other.centerpoint(b) as i32).abs() as u32;
            match (allele, b) {
                (&VariantType::SNV(a), &VariantType::SNV(b)) => {
                    if a == b && dist == 0 {
                        return Some(other.id(j));
                    }
                },
                (&VariantType::Insertion(l1), &VariantType::Insertion(l2)) |
                (&VariantType::Deletion(l1), &VariantType::Deletion(l2)) => {
                    if (l1 as i32 - l2 as i32).abs() as u32 <= max_len_diff && dist <= max_dist {
                        return Some(other.id(j));
                    }
                },
                // TODO: for now, ignore complex variants
                _ => continue
            }
        }
        None
    }

    pub fn id(&self, allele: usize) -> u32 {
        self.id + allele as u32
    }
}


#[derive(Debug)]
pub enum VariantType {
    SNV(u8),
    Insertion(u32),
    Deletion(u32),
    Unsupported
}


impl VariantType {
    pub fn is_unsupported(&self) -> bool {
        match self {
            &VariantType::Unsupported => true,
            _ => false
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum MatchError {
        MissingTag(tag: String) {
            description("missing tag")
            display("missing tag {}", tag)
        }
    }
}


pub struct RecordBuffer {
    reader: bcf::Reader,
    ringbuffer: VecDeque<Variant>,
    ringbuffer2: VecDeque<Variant>,
    window: u32,
    i: u32
}

impl RecordBuffer {
    pub fn new(reader: bcf::Reader, window: u32) -> Self {
        RecordBuffer {
            reader: reader,
            ringbuffer: VecDeque::new(),
            ringbuffer2: VecDeque::new(),
            window: window,
            i: 0
        }
    }

    pub fn last_rid(&self) -> Option<u32> {
        self.ringbuffer.back().map(|var| var.rid)
    }

    pub fn next_rid(&self) -> Option<u32> {
        self.ringbuffer2.back().map(|var| var.rid)
    }

    fn swap_buffers(&mut self) {
        // swap with buffer for next rid
        mem::swap(&mut self.ringbuffer2, &mut self.ringbuffer);
        // clear second buffer
        self.ringbuffer2.clear();
    }

    fn drain_left(&mut self, rid: u32, window_start: u32) {
        // remove records too far left or from wrong rid
        let to_remove = self.ringbuffer.iter().take_while(|var| var.pos < window_start || var.rid != rid).count();
        self.ringbuffer.drain(..to_remove);
    }

    pub fn fill(&mut self, chrom: &[u8], pos: u32) -> Result<(), Box<Error>> {
        let window_start = cmp::max(pos as i32 - self.window as i32 - 1, 0) as u32;
        let window_end = pos + self.window;
        let rid = try!(self.reader.header.name2rid(chrom));

        match (self.last_rid(), self.next_rid()) {
            (Some(last_rid), _) => {
                if last_rid != rid {
                    self.swap_buffers();
                } else {
                    self.drain_left(rid, window_start);
                }
            },
            (_, Some(_)) => {
                self.swap_buffers();
                self.drain_left(rid, window_start);
            },
            _ => ()
        }

        if !self.ringbuffer2.is_empty() {
            // We have already read beyond the current rid. Hence we can't extend to the right for
            // this rid.
            return Ok(())
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
            let alt_count = rec.alleles().len() as u32 - 1;
            if let Some(rec_rid) = rec.rid() {
                if rec_rid == rid {
                    if pos > window_end {
                        // Record is beyond our window. Store it anyways but stop.
                        self.ringbuffer.push_back(try!(Variant::new(&mut rec, &mut self.i)));
                        break;
                    } else if pos >= window_start {
                        // Record is within our window.
                        self.ringbuffer.push_back(try!(Variant::new(&mut rec, &mut self.i)));
                    } else {
                        // Record is upstream of our window, ignore it
                        self.i += alt_count;
                        continue
                    }
                } else if rec_rid > rid {
                    // record comes from next rid. Store it in second buffer but stop filling.
                    self.ringbuffer2.push_back(try!(Variant::new(&mut rec, &mut self.i)));
                    break;
                } else {
                    // Record comes from previous rid. Ignore it.
                    self.i += alt_count;
                    continue;
                }
            } else {
                // skip records without proper rid
                self.i += alt_count;
                continue;
            }
        }

        Ok(())
    }

    pub fn iter(&self) -> vec_deque::Iter<Variant> {
        self.ringbuffer.iter()
    }
}
