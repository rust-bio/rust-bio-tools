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
use crate::errors;
use itertools::Itertools;
use log::{info, warn};
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use snafu::ResultExt;
use std::collections::{btree_map, BTreeMap, HashMap};
use std::str;

pub struct VarIndex {
    inner: HashMap<Vec<u8>, BTreeMap<u32, Vec<Variant>>>,
    max_dist: u32,
}

impl VarIndex {
    pub fn new(mut reader: bcf::Reader, max_dist: u32) -> errors::Result<Self> {
        let mut inner = HashMap::new();
        let mut i = 0;
        let mut rec = reader.empty_record();
        loop {
            if let Err(e) = reader.read(&mut rec) {
                if e.is_eof() {
                    break;
                }
                return Err(e).context(errors::BCFReadError {
                    header: format!("{:?}", reader.header()),
                });
            }
            if let Some(rid) = rec.rid() {
                let chrom = reader
                    .header()
                    .rid2name(rid)
                    .context(errors::BCFReadIdError {
                        rid,
                        header: format!("{:?}", reader.header()),
                    })?;
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

pub fn match_variants(matchbcf: &str, max_dist: u32, max_len_diff: u32) -> errors::Result<()> {
    let mut inbcf = bcf::Reader::from_stdin().context(errors::BCFReaderStdinError {})?;
    let mut header = bcf::Header::from_template(inbcf.header());

    header.push_record(
        format!("##INFO=<ID=MATCHING,Number=A,Type=Integer,\
        Description=\"For each alternative allele, -1 if it does not match a variant in another VCF/BCF. \
        If it matches a variant, an id i>=0 points to the i-th variant in the VCF/BCF (counting each \
        alternative allele separately). For indels, matching is fuzzy: distance of centres <= {}, difference of \
        lengths <= {}\">", max_dist, max_len_diff).as_bytes()
    );
    let mut outbcf =
        bcf::Writer::from_stdout(&header, false, false).context(errors::BCFWriterStdoutError {
            header: format!("{:?}", header),
        })?;
    let index = VarIndex::new(
        bcf::Reader::from_path(matchbcf)
            .context(errors::BCFReaderFromPathError { path: matchbcf })?,
        max_dist,
    )?;

    let mut rec = inbcf.empty_record();
    let mut i = 0;
    loop {
        if let Err(e) = inbcf.read(&mut rec) {
            if e.is_eof() {
                break;
            }
            return Err(e).context(errors::BCFReadError {
                header: format!("{:?}", header),
            });
        }
        outbcf.translate(&mut rec);

        if let Some(rid) = rec.rid() {
            let chrom = inbcf
                .header()
                .rid2name(rid)
                .context(errors::BCFReadIdError {
                    rid,
                    header: format!("{:?}", inbcf.header()),
                })?;
            let pos = rec.pos();

            let var = Variant::new(&mut rec, &mut i)?;
            let matching = var
                .alleles
                .iter()
                .map(|a| {
                    if let Some(range) = index.range(chrom, pos) {
                        for v in Itertools::flatten(range.map(|(_, idx_vars)| idx_vars)) {
                            if let Some(id) = var.matches(v, a, max_dist, max_len_diff)? {
                                return Ok(id as i32);
                            }
                        }
                    }
                    Ok(-1)
                })
                .collect::<errors::Result<Vec<i32>>>()?;

            rec.push_info_integer(b"MATCHING", &matching)
                .context(errors::BCFTagWriteError {
                    data: format!("{:?}", matching),
                    fd: String::from("MATCHING"),
                })?;
            outbcf.write(&rec).context(errors::BCFWriteError)?;
        } else {
            outbcf.write(&rec).context(errors::BCFWriteError)?;
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
    pub fn new(rec: &mut bcf::Record, id: &mut u32) -> errors::Result<Self> {
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
                        return Err(errors::Error::MissingTagError {
                            tag: "SVLEN or END".to_string(),
                        });
                    }
                };
                VariantType::Deletion(svlen)
            } else {
                warn!(
                    "Unsupported variant {}",
                    str::from_utf8(&svtype).context(errors::StdStrUtf8Error)?
                );
                VariantType::Unsupported
            }]
        } else {
            let mut _alleles = Vec::with_capacity(alleles.len() - 1);
            for (i, a) in alleles[1..].iter().enumerate() {
                _alleles.push(if a == b"<DEL>" {
                    if let Some(ref svlens) = svlens {
                        VariantType::Deletion(svlens[i])
                    } else {
                        return Err(errors::Error::MissingTagError {
                            tag: "SVLEN".to_string(),
                        });
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
                        str::from_utf8(refallele).context(errors::StdStrUtf8Error)?,
                        str::from_utf8(a).context(errors::StdStrUtf8Error)?
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

    pub fn centerpoint(&self, allele: &VariantType) -> errors::Result<u32> {
        match allele {
            &VariantType::SNV(_) => Ok(self.pos),
            &VariantType::Insertion(_) => Ok(self.pos),
            &VariantType::Deletion(len) => Ok((self.pos as f64 + len as f64 / 2.0) as u32),
            &VariantType::Unsupported => Err(errors::Error::UnsupportedVariantError),
        }
    }

    pub fn matches(
        &self,
        other: &Variant,
        allele: &VariantType,
        max_dist: u32,
        max_len_diff: u32,
    ) -> errors::Result<Option<u32>> {
        if allele.is_unsupported() {
            return Ok(None);
        }
        for (j, b) in other.alleles.iter().enumerate() {
            if b.is_unsupported() {
                continue;
            }
            let dist =
                (self.centerpoint(allele)? as i32 - other.centerpoint(b)? as i32).abs() as u32;
            match (allele, b) {
                (&VariantType::SNV(a), &VariantType::SNV(b)) => {
                    if a == b && dist == 0 {
                        return Ok(Some(other.id(j)));
                    }
                }
                (&VariantType::Insertion(l1), &VariantType::Insertion(l2))
                | (&VariantType::Deletion(l1), &VariantType::Deletion(l2)) => {
                    if (l1 as i32 - l2 as i32).abs() as u32 <= max_len_diff && dist <= max_dist {
                        return Ok(Some(other.id(j)));
                    }
                }
                // TODO: for now, ignore complex variants
                _ => continue,
            }
        }
        Ok(None)
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
