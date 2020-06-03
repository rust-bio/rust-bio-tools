use crate::bcf::report::alignment_reader::{
    make_nucleobases, read_indexed_bam, AlignmentMatch, AlignmentNucleobase, Marker,
};
use crate::bcf::report::create_report::VariantType;
use serde::Serialize;
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Serialize, Clone, Debug)]
pub struct StaticAlignmentMatch {
    marker_type: Marker,
    start_position: f64,
    end_position: f64,
    flags: BTreeMap<u16, &'static str>,
    name: String,
    row: u8,
}

#[derive(Serialize, Clone)]
pub struct StaticAlignmentNucleobase {
    marker_type: Marker,
    bases: String,
    start_position: f64,
    end_position: f64,
    flags: BTreeMap<u16, &'static str>,
    name: String,
    row: u8,
}

#[derive(Serialize, Clone)]
pub struct StaticVariant {
    pub(crate) marker_type: String,
    pub(crate) reference: String,
    pub(crate) alternatives: Option<String>,
    pub(crate) start_position: f64,
    pub(crate) end_position: f64,
    pub(crate) row: i8,
    pub(crate) var_type: VariantType,
}

pub fn decode_static_flags(flag_vec: Vec<u16>) -> BTreeMap<u16, &'static str> {
    let mut string_map = BTreeMap::new();

    const FLAG_1: &'static str = "template having multiple segments in sequencing";
    const FLAG_2: &'static str = "each segment properly aligned according to the aligner";
    const FLAG_3: &'static str = "segment unmapped";
    const FLAG_4: &'static str = "next segment in the template unmapped";
    const FLAG_5: &'static str = "SEQ being reverse complemented";
    const FLAG_6: &'static str =
        "SEQ of the next segment in the template being reverse complemented";
    const FLAG_7: &'static str = "the first segment in the template ";
    const FLAG_8: &'static str = "the last segment in the template";
    const FLAG_9: &'static str = "secondary alignment";
    const FLAG_10: &'static str = "not passing filters, such as platform/vendor quality controls";
    const FLAG_11: &'static str = "PCR or optical duplicate";
    const FLAG_12: &'static str = "vega lite lines";

    let mut flags_map = BTreeMap::new();
    flags_map.insert(0x1, FLAG_1);
    flags_map.insert(0x2, FLAG_2);
    flags_map.insert(0x4, FLAG_3);
    flags_map.insert(0x8, FLAG_4);
    flags_map.insert(0x10, FLAG_5);
    flags_map.insert(0x20, FLAG_6);
    flags_map.insert(0x40, FLAG_7);
    flags_map.insert(0x80, FLAG_8);
    flags_map.insert(0x100, FLAG_9);
    flags_map.insert(0x200, FLAG_10);
    flags_map.insert(0x400, FLAG_11);
    flags_map.insert(0x800, FLAG_12);

    for (flag, text) in flags_map {
        for f in flag_vec.clone() {
            if (flag & f) == flag {
                string_map.insert(flag, text);
            }
        }
    }

    string_map
}

fn calc_rows(
    reads: Vec<AlignmentNucleobase>,
    matches: Vec<AlignmentMatch>,
) -> (Vec<StaticAlignmentNucleobase>, Vec<StaticAlignmentMatch>) {
    let mut row_ends = vec![0; 30];

    let mut read_names: BTreeMap<String, u8> = BTreeMap::new();

    let mut reads_wr: Vec<StaticAlignmentNucleobase> = Vec::new();
    let mut matches_wr: Vec<StaticAlignmentMatch> = Vec::new();

    for r in matches {
        let mut row: u8 = 0;

        if read_names.contains_key(&r.name) {
            row = *read_names.get(&r.name).unwrap();
        } else {
            for i in 1..30 {
                if r.read_start > row_ends[i] {
                    row = i as u8;
                    row_ends[i] = r.read_end;
                    read_names.insert(r.name.clone(), i as u8);
                    break;
                }
            }
        }

        let f = decode_static_flags(r.flags);

        let base = StaticAlignmentMatch {
            marker_type: r.marker_type,
            start_position: r.start_position,
            end_position: r.end_position,
            flags: f,
            name: r.name,
            row: row,
        };

        matches_wr.push(base);
    }

    for r in reads {
        let mut row: u8 = 0;

        if read_names.contains_key(&r.name) {
            row = *read_names.get(&r.name).unwrap();
        } else {
            for i in 1..30 {
                if r.read_start > row_ends[i] {
                    row = i as u8;
                    row_ends[i] = r.read_end;
                    read_names.insert(r.name.clone(), i as u8);
                    break;
                }
            }
        }

        let f = decode_static_flags(r.flags);

        let base = StaticAlignmentNucleobase {
            marker_type: r.marker_type,
            bases: r.bases,
            start_position: r.start_position,
            end_position: r.end_position,
            flags: f,
            name: r.name,
            row: row,
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
    let static_bases = calc_rows(msm, m);

    static_bases
}
