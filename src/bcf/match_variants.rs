//! Annotate for each variant in a VCF/BCF at STDIN whether it is contained in a given second VCF/BCF.
//!
//! The matching is fuzzy for indels and exact for SNVs.
//! Results are printed as BCF to STDOUT, with an additional INFO tag MATCHING.
//! The two vcfs do not have to be sorted.
//!
//! ## Usage:
//! ```bash
//! rbt vcf-match -d 50 -l 20 tests/test3.vcf < tests/test2.vcf > tests/matching.bcf
//! ```
//!
use itertools::Itertools;
use log::{info, warn};
use quick_error::quick_error;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use std::collections::{btree_map, BTreeMap, HashMap};
use std::error::Error;
use std::str;

pub struct VarIndex {
    inner: HashMap<Vec<u8>, BTreeMap<u32, Vec<Variant>>>,
    max_dist: u32,
}

impl VarIndex {
    pub fn new(mut reader: bcf::Reader, max_dist: u32) -> Result<Self, Box<dyn Error>> {
        let mut inner = HashMap::new();
        let mut i = 0;
        let mut rec = reader.empty_record();
        loop {
            if let Err(e) = reader.read(&mut rec) {
                if e.is_eof() {
                    break;
                }
                return Err(Box::new(e));
            }
            if let Some(rid) = rec.rid() {
                let chrom = reader.header().rid2name(rid)?;
                let recs = inner.entry(chrom.to_owned()).or_insert(BTreeMap::new());
                recs.entry(rec.pos())
                    .or_insert_with(|| Vec::new())
                    .push(Variant::new(&mut rec, &mut i)?);
            //recs.insert(rec.pos(), Variant::new(&mut rec, &mut i)?);
            } else {
                // skip records without rid
                let alt_count = rec.alleles().len() as u32 - 1;
                i += alt_count;
            }
        }

        Ok(VarIndex { inner, max_dist })
    }

    pub fn range(&self, chrom: &[u8], pos: u32) -> Option<btree_map::Range<'_, u32, Vec<Variant>>> {
        self.inner
            .get(chrom)
            .map(|recs| recs.range(pos.saturating_sub(self.max_dist)..pos + self.max_dist))
    }
}

pub fn match_variants(
    matchbcf: &str,
    max_dist: u32,
    max_len_diff: u32,
) -> Result<(), Box<dyn Error>> {
    let mut inbcf = bcf::Reader::from_stdin()?;
    let mut header = bcf::Header::from_template(inbcf.header());

    header.push_record(
        format!("##INFO=<ID=MATCHING,Number=A,Type=Integer,\
        Description=\"For each alternative allele, -1 if it does not match a variant in another VCF/BCF. \
        If it matches a variant, an id i>=0 points to the i-th variant in the VCF/BCF (counting each \
        alternative allele separately). For indels, matching is fuzzy: distance of centres <= {}, difference of \
        lengths <= {}\">", max_dist, max_len_diff).as_bytes()
    );
    let mut outbcf = bcf::Writer::from_path(&"-", &header, false, false)?;
    let index = VarIndex::new(bcf::Reader::from_path(matchbcf)?, max_dist)?;

    let mut rec = inbcf.empty_record();
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
            let chrom = inbcf.header().rid2name(rid)?;
            let pos = rec.pos();

            let var = Variant::new(&mut rec, &mut i)?;
            let matching = var
                .alleles
                .iter()
                .map(|a| {
                    if let Some(range) = index.range(chrom, pos) {
                        for v in Itertools::flatten(range.map(|(_, idx_vars)| idx_vars)) {
                            if let Some(id) = var.matches(v, a, max_dist, max_len_diff) {
                                return id as i32;
                            }
                        }
                    }
                    -1
                })
                .collect_vec();

            rec.push_info_integer(b"MATCHING", &matching)?;
            outbcf.write(&rec)?;
        } else {
            outbcf.write(&rec)?;
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
    alleles: Vec<VariantType>,
}

impl Variant {
    pub fn new(rec: &mut bcf::Record, id: &mut u32) -> Result<Self, Box<dyn Error>> {
        let pos = rec.pos();

        let svlens = if let Ok(Some(svlens)) = rec.info(b"SVLEN").integer() {
            Some(svlens.into_iter().map(|l| l.abs() as u32).collect_vec())
        } else {
            None
        };
        let svtype = if let Ok(Some(svtype)) = rec.info(b"SVTYPE").string() {
            Some(svtype[0].to_owned())
        } else {
            None
        };
        let end = if let Ok(Some(end)) = rec.info(b"END").integer() {
            Some(end[0] as u32)
        } else {
            None
        };
        let inslen = if let Ok(Some(inslen)) = rec.info(b"INSLEN").integer() {
            Some(inslen[0] as u32)
        } else {
            None
        };
        let alleles = rec.alleles();
        let refallele = alleles[0];

        let _alleles: Vec<VariantType> = if let Some(svtype) = svtype {
            vec![if svtype == b"INS" {
                match (svlens, inslen) {
                    (Some(svlens), _) => VariantType::Insertion(svlens[0]),
                    (None, Some(inslen)) => VariantType::Insertion(inslen),
                    _ => {
                        warn!("Unsupported variant INS without SVLEN or INSLEN");
                        VariantType::Unsupported
                    }
                }
            } else if svtype == b"DEL" {
                let svlen = match (svlens, end) {
                    (Some(svlens), _) => svlens[0],
                    (None, Some(end)) => end - 1 - pos,
                    _ => {
                        return Err(Box::new(MatchError::MissingTag("SVLEN or END".to_owned())));
                    }
                };
                VariantType::Deletion(svlen)
            } else {
                warn!("Unsupported variant {}", str::from_utf8(&svtype)?);
                VariantType::Unsupported
            }]
        } else {
            let mut _alleles = Vec::with_capacity(alleles.len() - 1);
            for (i, a) in alleles[1..].iter().enumerate() {
                _alleles.push(if a == b"<DEL>" {
                    if let Some(ref svlens) = svlens {
                        VariantType::Deletion(svlens[i])
                    } else {
                        return Err(Box::new(MatchError::MissingTag("SVLEN".to_owned())));
                    }
                } else if a.len() < refallele.len() {
                    VariantType::Deletion((refallele.len() - a.len()) as u32)
                } else if a.len() > refallele.len() {
                    VariantType::Insertion((a.len() - refallele.len()) as u32)
                } else if a.len() == 1 {
                    VariantType::SNV(a[0])
                } else {
                    warn!(
                        "Unsupported variant {} -> {}",
                        str::from_utf8(refallele)?,
                        str::from_utf8(a)?
                    );
                    VariantType::Unsupported
                });
            }
            _alleles
        };
        let var = Variant {
            id: *id,
            rid: rec.rid().unwrap(),
            pos,
            alleles: _alleles,
        };
        *id += alleles.len() as u32 - 1;
        Ok(var)
    }

    pub fn centerpoint(&self, allele: &VariantType) -> u32 {
        match allele {
            &VariantType::SNV(_) => self.pos,
            &VariantType::Insertion(_) => self.pos,
            &VariantType::Deletion(len) => (self.pos as f64 + len as f64 / 2.0) as u32,
            &VariantType::Unsupported => panic!("Unsupported variant."),
        }
    }

    pub fn matches(
        &self,
        other: &Variant,
        allele: &VariantType,
        max_dist: u32,
        max_len_diff: u32,
    ) -> Option<u32> {
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
                }
                (&VariantType::Insertion(l1), &VariantType::Insertion(l2))
                | (&VariantType::Deletion(l1), &VariantType::Deletion(l2)) => {
                    if (l1 as i32 - l2 as i32).abs() as u32 <= max_len_diff && dist <= max_dist {
                        return Some(other.id(j));
                    }
                }
                // TODO: for now, ignore complex variants
                _ => continue,
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
    Unsupported,
}

impl VariantType {
    pub fn is_unsupported(&self) -> bool {
        match self {
            &VariantType::Unsupported => true,
            _ => false,
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
