use std::collections::HashMap;
use std::path::Path;

use anyhow::Context;
use anyhow::Result;
use itertools::Itertools;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;

pub fn split<P: AsRef<Path>>(input_bcf: P, output_bcfs: &[P]) -> Result<()> {
    let n_records = bcf::Reader::from_path(input_bcf.as_ref())
        .context("error reading input VCF/BCF")?
        .records()
        .fold(0_u64, |count, _| count + 1);
    let mut reader = bcf::Reader::from_path(input_bcf).context("error reading input VCF/BCF")?;
    let header = bcf::Header::from_template(reader.header());
    let mut bnd_cache = HashMap::new();

    let chunk_size = n_records / output_bcfs.len() as u64;

    let mut writers = output_bcfs
        .iter()
        .map(|path| {
            bcf::Writer::from_path(path, &header, false, bcf::Format::Bcf)
                .context("error creating output VCF/BCF")
        })
        .collect::<Result<Vec<_>>>()?;

    for (rec, i) in reader.records().zip(0..) {
        let rec = rec?;

        let mut chunk = i / (chunk_size + 1);
        if rec.is_bnd() {
            if let Some(group) = BreakendGroup::from(&rec) {
                let event_chunk = match group {
                    BreakendGroup::Event(id) => bnd_cache.entry(id).or_insert(chunk),
                    BreakendGroup::Mates(ids) => {
                        let ids = ids.clone();
                        bnd_cache.entry(ids.concat()).or_insert(chunk)
                    }
                };
                chunk = *event_chunk;
            }
        };
        let writer = &mut writers[chunk as usize];
        writer.write(&rec)?;
    }

    Ok(())
}

#[derive(Eq, PartialEq, Hash, Clone, Debug)]
enum BreakendGroup {
    Event(Vec<u8>),
    Mates(Vec<Vec<u8>>),
}

impl BreakendGroup {
    fn from(rec: &bcf::Record) -> Option<Self> {
        if let Some(event) = rec.event() {
            Some(BreakendGroup::Event(event))
        } else if let Some(mut mates) = rec.mateids() {
            let id = rec.id();
            mates.push(id);
            mates.sort();
            Some(BreakendGroup::Mates(mates))
        } else {
            None
        }
    }
}

type Id = Vec<u8>;

trait BndRecord {
    fn is_bnd(&self) -> bool;
    fn event(&self) -> Option<Id>;
    fn mateids(&self) -> Option<Vec<Id>>;
}

impl BndRecord for bcf::Record {
    fn is_bnd(&self) -> bool {
        self.info(b"SVTYPE").string().map_or(false, |entries| {
            entries.map_or(false, |entries| entries[0] == b"BND")
        })
    }

    fn event(&self) -> Option<Id> {
        if let Ok(Some(event)) = self.info(b"EVENT").string() {
            Some(event[0].to_owned())
        } else {
            None
        }
    }

    fn mateids(&self) -> Option<Vec<Id>> {
        match self.info(b"MATEID").string() {
            Ok(Some(s)) => Some(s.clone().into_iter().map(|v| v.to_vec()).collect_vec()),
            _ => None,
        }
    }
}
