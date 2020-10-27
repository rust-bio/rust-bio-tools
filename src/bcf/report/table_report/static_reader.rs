use crate::bcf::report::table_report::alignment_reader::{
    make_nucleobases, read_indexed_bam, AlignmentMatch, AlignmentNucleobase,
};
use crate::bcf::report::table_report::create_report_table::VariantType;
use serde::Serialize;
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Serialize, Clone, Debug)]
pub struct StaticAlignmentMatch {
    #[serde(flatten)]
    alignment: AlignmentMatch,
    row: u8,
}

#[derive(Serialize, Clone)]
pub struct StaticAlignmentNucleobase {
    #[serde(flatten)]
    nucleobase: AlignmentNucleobase,
    row: u8,
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
) -> (Vec<StaticAlignmentNucleobase>, Vec<StaticAlignmentMatch>) {
    let mut row_ends = vec![0; 1000];

    let mut read_names: BTreeMap<String, u8> = BTreeMap::new();

    let mut reads_wr: Vec<StaticAlignmentNucleobase> = Vec::new();
    let mut matches_wr: Vec<StaticAlignmentMatch> = Vec::new();

    for r in matches {
        let mut row: u8 = 0;

        if read_names.contains_key(&r.name) {
            row = *read_names.get(&r.name).unwrap();
        } else {
            for (i, _) in row_ends.iter().enumerate().take(1000).skip(1) {
                if r.read_start > row_ends[i] {
                    row = i as u8;
                    row_ends[i] = r.read_end;
                    read_names.insert(r.name.clone(), i as u8);
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

    for r in reads {
        let mut row: u8 = 0;

        if read_names.contains_key(&r.name) {
            row = *read_names.get(&r.name).unwrap();
        } else {
            for (i, _) in row_ends.iter().enumerate().take(1000).skip(1) {
                if r.read_start > row_ends[i] {
                    row = i as u8;
                    row_ends[i] = r.read_end;
                    read_names.insert(r.name.clone(), i as u8);
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

    (reads_wr, matches_wr)
}

pub fn get_static_reads(
    path: &Path,
    fasta_path: &Path,
    chrom: String,
    from: u64,
    to: u64,
) -> (Vec<StaticAlignmentNucleobase>, Vec<StaticAlignmentMatch>) {
    let alignments = read_indexed_bam(path, chrom.clone(), from, to);
    let (msm, m) = make_nucleobases(fasta_path, chrom, alignments, from, to);
    calc_rows(msm, m)
}
