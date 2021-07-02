use crate::bcf::report::table_report::alignment_reader::{
    make_nucleobases, read_indexed_bam, AlignmentMatch, AlignmentNucleobase,
};
use crate::bcf::report::table_report::create_report_table::VariantType;
use crate::common::Region;
use anyhow::Result;
use rand::rngs::StdRng;
use rand::seq::IteratorRandom;
use rand_core::SeedableRng;
use serde::Serialize;
use std::collections::{BTreeMap, HashSet};
use std::path::Path;

#[derive(Serialize, Clone, Debug)]
pub struct StaticAlignmentMatch {
    #[serde(flatten)]
    alignment: AlignmentMatch,
    row: u16,
}

#[derive(Serialize, Clone)]
pub struct StaticAlignmentNucleobase {
    #[serde(flatten)]
    nucleobase: AlignmentNucleobase,
    row: u16,
}

#[derive(Serialize, Clone)]
pub struct Variant {
    pub(crate) marker_type: String,
    pub(crate) reference: String,
    pub(crate) alternatives: Option<String>,
    pub(crate) start_position: f64,
    pub(crate) end_position: f64,
    pub(crate) row: i8,
    pub(crate) var_type: VariantType,
}

fn calc_rows(
    reads: Vec<AlignmentNucleobase>,
    matches: Vec<AlignmentMatch>,
    max_read_depth: u32,
    variant: Option<&Variant>,
) -> (
    Vec<StaticAlignmentNucleobase>,
    Vec<StaticAlignmentMatch>,
    usize,
) {
    let mut row_ends = vec![0; 10000];

    let mut read_names: BTreeMap<String, u16> = BTreeMap::new();

    let mut reads_wr: Vec<StaticAlignmentNucleobase> = Vec::new();
    let mut matches_wr: Vec<StaticAlignmentMatch> = Vec::new();

    let mut max_row = 0;

    for r in matches {
        let overlaps = if variant.is_some() {
            r.read_start < variant.unwrap().start_position as u32
                && r.read_end > variant.unwrap().end_position as u32
        } else {
            true
        };

        if overlaps {
            let mut row: u16 = 0;

            if read_names.contains_key(&r.name) {
                row = *read_names.get(&r.name).unwrap();
            } else {
                for (i, _) in row_ends.iter().enumerate().take(10000).skip(1) {
                    if r.read_start > row_ends[i] {
                        if i > max_row {
                            max_row = i;
                        }
                        row = i as u16;
                        row_ends[i] = r.read_end;
                        read_names.insert(r.name.clone(), i as u16);
                        break;
                    }
                }
            }

            let base = StaticAlignmentMatch {
                alignment: r.clone(),
                row,
            };

            matches_wr.push(base);
        }
    }

    for r in reads {
        let overlaps = if variant.is_some() {
            r.read_start < variant.unwrap().start_position as u32
                && r.read_end > variant.unwrap().end_position as u32
        } else {
            true
        };

        if overlaps {
            let mut row: u16 = 0;

            if read_names.contains_key(&r.name) {
                row = *read_names.get(&r.name).unwrap();
            } else {
                for (i, _) in row_ends.iter().enumerate().take(10000).skip(1) {
                    if r.read_start > row_ends[i] {
                        if i > max_row {
                            max_row = i;
                        }
                        row = i as u16;
                        row_ends[i] = r.read_end;
                        read_names.insert(r.name.clone(), i as u16);
                        break;
                    }
                }
            }

            let base = StaticAlignmentNucleobase {
                nucleobase: r.clone(),
                row,
            };

            reads_wr.push(base);
        }
    }

    if max_row > max_read_depth as usize {
        let mut rng = StdRng::seed_from_u64(42);
        let random_rows: HashSet<_> = (0..max_row as u32)
            .choose_multiple(&mut rng, max_read_depth as usize)
            .into_iter()
            .collect();
        reads_wr = reads_wr
            .into_iter()
            .filter(|b| random_rows.contains(&(b.row as u32)))
            .collect();
        matches_wr = matches_wr
            .into_iter()
            .filter(|b| random_rows.contains(&(b.row as u32)))
            .collect();
        max_row = max_read_depth as usize;
    }

    (reads_wr, matches_wr, max_row)
}

pub fn get_static_reads<P: AsRef<Path>>(
    path: P,
    fasta_path: P,
    region: &Region,
    max_read_depth: u32,
    variant: Option<&Variant>,
) -> Result<(
    Vec<StaticAlignmentNucleobase>,
    Vec<StaticAlignmentMatch>,
    usize,
)> {
    let alignments = read_indexed_bam(path, region)?;
    let (msm, m) = make_nucleobases(fasta_path, region, alignments)?;
    Ok(calc_rows(msm, m, max_read_depth, variant))
}
