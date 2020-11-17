use std::collections::HashMap;
use std::iter;
use std::path::Path;

use anyhow::Context;
use anyhow::Result;
use itertools::Itertools;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;

pub fn split<P: AsRef<Path>>(input_bcf: P, output_bcfs: &[P]) -> Result<()> {
    let info = BCFInfo::new(&input_bcf).context("error reading input VCF/BCF")?;
    let mut reader = bcf::Reader::from_path(input_bcf).context("error reading input VCF/BCF")?;
    let header = bcf::Header::from_template(reader.header());
    let mut bnd_cache = HashMap::new();

    let chunk_size = info.n_records / output_bcfs.len() as u64;

    let mut i = 0;
    for (chunk, out_path) in output_bcfs.iter().enumerate() {
        let chunk = chunk as u64;
        let mut writer = bcf::Writer::from_path(out_path, &header, false, bcf::Format::BCF)?;
        let mut written = 0;
        let write_err = || {
            format!(
                "error writing record to {}",
                out_path.as_ref().as_os_str().to_str().unwrap()
            )
        };

        loop {
            let mut rec = reader.empty_record();
            let read_err = || format!("error reading {}-th record of input VCF/BCF", i + 1);

            if !reader.read(&mut rec).with_context(read_err)? {
                // EOF
                break;
            }

            let towrite = if is_bnd(&mut rec).with_context(read_err)? {
                if let Some(group) = BreakendGroup::from(&mut rec) {
                    if let Some(end) = info.end(&group) {
                        // BND is part of a group.
                        if i == end {
                            // BND is last BND of group.
                            // Write out all previous BNDs of group.
                            if let Some(bnds) =
                                bnd_cache.get(info.supergroups.get(&group).unwrap_or(&group))
                            {
                                for rec in bnds {
                                    writer.write(rec).with_context(write_err)?;
                                    written += 1;
                                }
                            }
                            Some(rec)
                        } else {
                            // Cache BND record for later.
                            let entry = bnd_cache.entry(group).or_insert_with(Vec::new);
                            entry.push(rec);
                            None
                        }
                    } else {
                        Some(rec)
                    }
                } else {
                    Some(rec)
                }
            } else {
                Some(rec)
            };

            if let Some(towrite) = towrite {
                writer.write(&towrite).with_context(write_err)?;
                written += 1;
            }

            i += 1;

            if chunk < output_bcfs.len() as u64 - 1 && written >= chunk_size {
                // go on with next chunk
                break;
            }
        }
    }

    Ok(())
}

#[derive(Eq, PartialEq, Hash, Clone, Debug)]
enum BreakendGroup {
    Event(Vec<u8>),
    Mates(Vec<Vec<u8>>),
}

impl BreakendGroup {
    fn subgroups<'a>(&'a self) -> Box<dyn Iterator<Item = Self> + 'a> {
        if let BreakendGroup::Mates(mates) = self {
            if mates.len() > 2 {
                return Box::new(
                    (2..mates.len())
                        .map(move |k| {
                            if let BreakendGroup::Mates(mates) = self {
                                mates.iter().cloned().combinations(k)
                            } else {
                                unreachable!();
                            }
                        })
                        .flatten()
                        .map(BreakendGroup::Mates),
                );
            }
        }
        Box::new(iter::empty())
    }

    fn len(&self) -> usize {
        match self {
            BreakendGroup::Event(_) => 1,
            BreakendGroup::Mates(mates) => mates.len(),
        }
    }

    fn from(rec: &mut bcf::Record) -> Option<Self> {
        if let Some(event) = event(rec) {
            Some(BreakendGroup::Event(event))
        } else if let Some(mateids) = mateids(rec) {
            let mut mates: Vec<_> = mateids
                .into_iter()
                .map(|mateid| mateid.to_owned())
                .collect();
            let id = rec.id();
            mates.push(id);
            mates.sort();
            Some(BreakendGroup::Mates(mates))
        } else {
            None
        }
    }
}

#[derive(Debug)]
struct BCFInfo {
    n_records: u64,
    bnd_ends: HashMap<BreakendGroup, u64>,
    supergroups: HashMap<BreakendGroup, BreakendGroup>,
}

impl BCFInfo {
    fn new<P: AsRef<Path>>(input_bcf: P) -> Result<Self> {
        let mut reader = bcf::Reader::from_path(input_bcf)?;
        let mut bnd_ends = HashMap::new();
        let mut supergroups: HashMap<BreakendGroup, BreakendGroup> = HashMap::new();

        let mut i = 0;
        let mut set_end = |group: BreakendGroup, i: u64| {
            let end = bnd_ends.entry(group).or_insert(0);
            *end = i;
        };
        for res in reader.records() {
            let mut rec = res?;

            if is_bnd(&mut rec)? {
                if let Some(group) = BreakendGroup::from(&mut rec) {
                    if let Some(repr) = supergroups.get(&group) {
                        // this group is part of a supergroup
                        set_end(repr.clone(), i)
                    } else {
                        // register all subgroups
                        for subgroup in group.subgroups() {
                            if let Some(repr) = supergroups.get(&subgroup) {
                                if repr.len() < group.len() {
                                    supergroups.insert(subgroup, group.clone());
                                }
                            }
                        }
                        set_end(group, i)
                    }
                }
            }
            i += 1;
        }

        Ok(BCFInfo {
            n_records: i,
            bnd_ends,
            supergroups,
        })
    }

    fn end(&self, group: &BreakendGroup) -> Option<u64> {
        self.bnd_ends
            .get(self.supergroups.get(group).unwrap_or(group))
            .copied()
    }
}

fn is_bnd(record: &mut bcf::Record) -> Result<bool> {
    Ok(record.info(b"SVTYPE").string().map_or(false, |entries| {
        entries.map_or(false, |entries| entries[0] == b"BND")
    }))
}

fn event(record: &mut bcf::Record) -> Option<Vec<u8>> {
    if let Ok(Some(event)) = record.info(b"EVENT").string() {
        Some(event[0].to_owned())
    } else {
        None
    }
}

fn mateids(record: &mut bcf::Record) -> Option<Vec<&[u8]>> {
    record.info(b"MATEID").string().unwrap_or(None)
}
