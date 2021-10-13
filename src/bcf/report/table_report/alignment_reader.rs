extern crate rust_htslib;

use self::rust_htslib::bam::FetchDefinition;
use crate::bcf::report::table_report::fasta_reader::read_fasta;
use crate::common::Region;
use anyhow::Result;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::{bam, bam::Read};
use serde::Serialize;
use std::path::Path;

#[derive(Serialize, Clone, Debug, PartialEq)]
pub enum Marker {
    A,
    T,
    C,
    G,
    N,
    Deletion,
    Insertion,
    Match,
    Pairing,
}

#[derive(Clone, Debug)]
pub struct Alignment {
    sequence: String,
    pos: i64,
    length: u16,
    flags: Vec<u16>,
    name: String,
    cigar: CigarStringView,
    paired: bool,
    mate_pos: i64,
    tid: i32,
    mate_tid: i32,
    mapq: u8,
}

#[derive(Serialize, Clone, Debug, PartialEq)]
pub struct AlignmentNucleobase {
    pub marker_type: Marker,
    pub bases: String,
    pub start_position: f64,
    pub end_position: f64,
    pub flags: Vec<u16>,
    pub name: String,
    pub read_start: u32,
    pub read_end: u32,
    pub mapq: u8,
    pub cigar: String,
}

#[derive(Serialize, Clone, Debug, PartialEq)]
pub struct AlignmentMatch {
    pub marker_type: Marker,
    pub start_position: f64,
    pub end_position: f64,
    pub flags: Vec<u16>,
    pub name: String,
    pub read_start: u32,
    pub read_end: u32,
    pub mapq: u8,
    pub cigar: String,
}

pub fn decode_flags(code: u16) -> Vec<u16> {
    let flags_map = vec![
        0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 0x100, 0x200, 0x400, 0x800,
    ];
    let mut read_map = Vec::new();

    for flag in flags_map {
        if (flag & code) == flag {
            read_map.push(flag);
        }
    }

    read_map
}

pub fn read_indexed_bam<P: AsRef<Path>>(path: P, region: &Region) -> Result<Vec<Alignment>> {
    let mut bam = bam::IndexedReader::from_path(&path)?;
    let chrom = &region.target;
    let (from, to) = (region.start, region.end);
    let tid = bam.header().tid(chrom.as_bytes()).unwrap() as i32;

    let mut alignments: Vec<Alignment> = Vec::new();

    bam.fetch(FetchDefinition::Region(tid, from as i64, to as i64))?;

    for r in bam.records() {
        let rec = r?;
        let a = make_alignment(&rec);
        alignments.push(a);
    }

    Ok(alignments)
}

fn make_alignment(record: &bam::Record) -> Alignment {
    let has_pair = record.is_paired();

    let mate_pos = record.mpos();

    let tid = record.tid();
    let mtid = record.mtid();

    //Cigar String
    let cigstring = record.cigar();

    //Position
    let pos = record.pos();

    //Länge
    let le = record.seq().len() as u16;

    //Sequenz
    let seq = record.seq().as_bytes();
    let sequenz = String::from_utf8(seq).unwrap();

    //Flags
    let flgs = record.flags();
    let flag_string = decode_flags(flgs);

    //Name
    let n = record.qname().to_owned();
    let name = String::from_utf8(n).unwrap();

    Alignment {
        sequence: sequenz,
        pos,
        length: le,
        cigar: cigstring,
        flags: flag_string,
        name,
        paired: has_pair,
        mate_pos,
        tid,
        mate_tid: mtid,
        mapq: record.mapq(),
    }
}

pub fn make_nucleobases<P: AsRef<Path>>(
    fasta_path: P,
    region: &Region,
    snippets: Vec<Alignment>,
) -> Result<(Vec<AlignmentNucleobase>, Vec<AlignmentMatch>)> {
    let mut bases: Vec<AlignmentNucleobase> = Vec::new();
    let mut matches: Vec<AlignmentMatch> = Vec::new();

    let ref_bases = read_fasta(fasta_path, region, false)?;
    let (from, to) = (region.start, region.end);

    for snippet in snippets {
        let mut cigar_offset: i64 = 0;
        let mut read_offset: i64 = 0;
        let base_string = snippet.sequence.clone();
        let char_vec: Vec<char> = base_string.chars().collect();

        let mut soft_clip_begin = true;

        let temp_snippet = snippet.clone();

        if temp_snippet.paired
            && (temp_snippet.pos + temp_snippet.length as i64) < temp_snippet.mate_pos
            && temp_snippet.tid == temp_snippet.mate_tid
        {
            let pairing = AlignmentMatch {
                marker_type: Marker::Pairing,
                start_position: (temp_snippet.pos + temp_snippet.length as i64) as f64 - 0.5,
                end_position: temp_snippet.mate_pos as f64 - 0.5,
                flags: temp_snippet.flags.clone(),
                name: temp_snippet.name.clone(),
                read_start: temp_snippet.pos as u32,
                read_end: (temp_snippet.mate_pos + 100) as u32,
                mapq: snippet.mapq,
                cigar: snippet.cigar.to_string(),
            };

            matches.push(pairing);
        }

        for c in snippet.cigar.iter() {
            let mut match_count = 0;
            let mut match_start = 0;
            let mut match_ending = false;

            match c {
                rust_htslib::bam::record::Cigar::Match(c) => {
                    for _i in 0..rust_htslib::bam::record::Cigar::Match(*c).len() {
                        let snip = snippet.clone();
                        let b = char_vec[cigar_offset as usize];

                        if snip.pos + read_offset >= from as i64
                            && snip.pos + read_offset < to as i64
                        {
                            let ref_index = snip.pos + read_offset - from as i64;
                            let ref_base = &ref_bases[ref_index as usize];

                            if ref_base.get_marker_type() == b {
                                // Create long rule while bases match
                                if match_count == 0 {
                                    match_start = snip.pos as i64 + read_offset;
                                }
                                match_count += 1;
                                match_ending = true;

                            //m = Marker::Match; // Match with reference fasta
                            } else {
                                let (mtch, base) = make_markers(
                                    snip.clone(),
                                    b,
                                    read_offset,
                                    match_start,
                                    match_count,
                                );
                                if let Some(m) = mtch {
                                    matches.push(m)
                                }
                                bases.push(base);

                                match_count = 0;
                                match_start = 0;

                                match_ending = false;
                            }
                        }
                        cigar_offset += 1;
                        read_offset += 1;
                    }

                    if match_ending {
                        let mtch =
                            end_mismatch_detection(snippet.clone(), match_start, match_count);
                        matches.push(mtch);
                    }

                    soft_clip_begin = false;
                }
                rust_htslib::bam::record::Cigar::Ins(c) => {
                    let snip = snippet.clone();
                    let p: f64 = snip.pos as f64 + read_offset as f64 - 0.5;
                    let m: Marker = Marker::Insertion;
                    let rs;
                    let re;

                    let mut b = String::from("");
                    for _i in 0..rust_htslib::bam::record::Cigar::Ins(*c).len() {
                        let char = char_vec[cigar_offset as usize];
                        b.push(char);

                        cigar_offset += 1;
                    }

                    if snip.paired && snip.tid == snip.mate_tid {
                        if snip.pos < snip.mate_pos {
                            re = snip.mate_pos + 100;
                            rs = snip.pos;
                        } else {
                            rs = snip.mate_pos;
                            re = snip.pos + snip.length as i64;
                        }
                    } else {
                        rs = snip.pos;
                        re = snip.pos + snip.length as i64;
                    }

                    let base = AlignmentNucleobase {
                        marker_type: m,
                        bases: b,
                        start_position: p as f64 + 0.5,
                        end_position: p as f64 + 1.5,
                        flags: snip.flags,
                        name: snip.name,
                        read_start: rs as u32,
                        read_end: re as u32,
                        mapq: snippet.mapq,
                        cigar: snippet.cigar.to_string(),
                    };

                    if from as f64 <= (base.start_position + 0.5)
                        && (base.start_position + 0.5) <= to as f64
                    {
                        bases.push(base);
                    }

                    soft_clip_begin = false;
                }
                rust_htslib::bam::record::Cigar::Del(c) => {
                    for _i in 0..rust_htslib::bam::record::Cigar::Del(*c).len() {
                        let snip = snippet.clone();
                        let marker = Marker::Deletion;
                        let position = snip.pos as i64 + read_offset;
                        let flags = snip.flags;
                        let name = snip.name;
                        let read_start;
                        let read_end;
                        let empty_bases = String::from("");

                        if snip.paired && snip.tid == snip.mate_tid {
                            if snip.pos < snip.mate_pos {
                                read_end = snip.mate_pos + 100;
                                read_start = snip.pos;
                            } else {
                                read_start = snip.mate_pos;
                                read_end = snip.pos + snip.length as i64;
                            }
                        } else {
                            read_start = snip.pos;
                            read_end = snip.pos + snip.length as i64;
                        }

                        let base = AlignmentNucleobase {
                            marker_type: marker,
                            bases: empty_bases,
                            start_position: position as f64 + 0.5,
                            end_position: position as f64 + 1.5,
                            flags,
                            name,
                            read_start: read_start as u32,
                            read_end: read_end as u32,
                            mapq: snippet.mapq,
                            cigar: snippet.cigar.to_string(),
                        };

                        read_offset += 1;

                        if from as f64 <= (base.start_position + 0.5)
                            && (base.start_position + 0.5) <= to as f64
                        {
                            bases.push(base);
                        }
                    }

                    soft_clip_begin = false;
                }
                rust_htslib::bam::record::Cigar::SoftClip(c) => {
                    for _i in 0..rust_htslib::bam::record::Cigar::SoftClip(*c).len() {
                        let snip = snippet.clone();
                        let b = char_vec[cigar_offset as usize];

                        if snip.pos + read_offset >= from as i64
                            && snip.pos + read_offset < to as i64
                        {
                            let ref_index = snip.pos + read_offset - from as i64;
                            let ref_base = &ref_bases[ref_index as usize];

                            if ref_base.get_marker_type() == b {
                                // Create long rule while bases match
                                if match_count == 0 {
                                    match_start = snip.pos as i64 + read_offset;
                                }
                                match_count += 1;
                                match_ending = true;
                            } else {
                                let (mtch, base) = make_markers(
                                    snip.clone(),
                                    b,
                                    read_offset,
                                    match_start,
                                    match_count,
                                );
                                if let Some(m) = mtch {
                                    matches.push(m)
                                }
                                bases.push(base);

                                match_count = 0;
                                match_start = 0;

                                match_ending = false;
                            }
                        }

                        cigar_offset += 1;
                        if !soft_clip_begin {
                            read_offset += 1;
                        }
                    }

                    if match_ending {
                        let mtch =
                            end_mismatch_detection(snippet.clone(), match_start, match_count);
                        matches.push(mtch);
                    }

                    soft_clip_begin = false;
                }
                _ => {
                    soft_clip_begin = false;
                }
            }
        }
    }
    Ok((bases, matches))
}

fn make_markers(
    snip: Alignment,
    base: char,
    read_offset: i64,
    match_start: i64,
    match_count: i64,
) -> (Option<AlignmentMatch>, AlignmentNucleobase) {
    let marker: Marker;

    match base {
        // Mismatch
        'A' => marker = Marker::A,
        'T' => marker = Marker::T,
        'C' => marker = Marker::C,
        'N' => marker = Marker::N,
        'G' => marker = Marker::G,
        _ => marker = Marker::Deletion,
    }

    let position = snip.pos as i64 + read_offset;
    let flags = snip.flags;
    let name = snip.name;

    let read_start: i64;
    let read_end: i64;

    if snip.paired && snip.tid == snip.mate_tid {
        if snip.pos < snip.mate_pos {
            read_end = snip.mate_pos + 100;
            read_start = snip.pos;
        } else {
            read_start = snip.mate_pos;
            read_end = snip.pos + snip.length as i64;
        }
    } else {
        read_start = snip.pos;
        read_end = snip.pos + snip.length as i64;
    }

    let mut mtch = None;

    if match_count > 0 {
        // First mismatch detection must lead to new creation of all previous matches
        mtch = Some(AlignmentMatch {
            marker_type: Marker::Match,
            start_position: match_start as f64 + 0.5,
            end_position: (match_start + match_count - 1) as f64 + 1.5,
            flags: flags.clone(),
            name: name.clone(),
            read_start: read_start as u32,
            read_end: read_end as u32,
            mapq: snip.mapq,
            cigar: snip.cigar.to_string(),
        });
    }

    let base = AlignmentNucleobase {
        marker_type: marker,
        bases: base.to_string(),
        start_position: position as f64 + 0.5,
        end_position: position as f64 + 1.5,
        flags,
        name,
        read_start: read_start as u32,
        read_end: read_end as u32,
        mapq: snip.mapq,
        cigar: snip.cigar.to_string(),
    };
    (mtch, base)
}

fn end_mismatch_detection(snip: Alignment, match_start: i64, match_count: i64) -> AlignmentMatch {
    let f = snip.flags;
    let n = snip.name;

    let rs: i64;
    let re: i64;

    if snip.paired && snip.tid == snip.mate_tid {
        if snip.pos < snip.mate_pos {
            re = snip.mate_pos + 100;
            rs = snip.pos;
        } else {
            rs = snip.mate_pos;
            re = snip.pos + snip.length as i64;
        }
    } else {
        rs = snip.pos;
        re = snip.pos + snip.length as i64;
    }

    AlignmentMatch {
        marker_type: Marker::Match,
        start_position: match_start as f64 + 0.5,
        end_position: (match_start + match_count - 1) as f64 + 1.5,
        flags: f,
        name: n,
        read_start: rs as u32,
        read_end: re as u32,
        mapq: snip.mapq,
        cigar: snip.cigar.to_string(),
    }
}
