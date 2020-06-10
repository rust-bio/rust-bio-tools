extern crate rust_htslib;

use crate::bcf::report::fasta_reader::read_fasta;
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
}

pub fn decode_flags(code: u16) -> Vec<u16> {
    let mut flags_map = Vec::new();
    flags_map.push(0x1);
    flags_map.push(0x2);
    flags_map.push(0x4);
    flags_map.push(0x8);
    flags_map.push(0x10);
    flags_map.push(0x20);
    flags_map.push(0x40);
    flags_map.push(0x80);
    flags_map.push(0x100);
    flags_map.push(0x200);
    flags_map.push(0x400);
    flags_map.push(0x800);

    let mut read_map = Vec::new();

    for flag in flags_map {
        if (flag & code.clone()) == flag {
            read_map.push(flag);
        }
    }

    read_map
}

pub fn read_indexed_bam(path: &Path, chrom: String, from: u64, to: u64) -> Vec<Alignment> {
    let mut bam = bam::IndexedReader::from_path(&path).unwrap();
    let tid = bam.header().tid(chrom.as_bytes()).unwrap();

    let mut alignments: Vec<Alignment> = Vec::new();

    bam.fetch(tid, from, to).unwrap();

    for r in bam.records() {
        let rec = r.unwrap();

        let a = make_alignment(&rec);

        alignments.push(a);
    }

    alignments
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

    let read = Alignment {
        sequence: sequenz,
        pos: pos,
        length: le,
        cigar: cigstring,
        flags: flag_string,
        name: name,
        paired: has_pair,
        mate_pos: mate_pos,
        tid: tid,
        mate_tid: mtid,
    };

    read
}

pub fn make_nucleobases(
    fasta_path: &Path,
    chrom: String,
    snippets: Vec<Alignment>,
    from: u64,
    to: u64,
) -> (Vec<AlignmentNucleobase>, Vec<AlignmentMatch>) {
    let mut bases: Vec<AlignmentNucleobase> = Vec::new();
    let mut matches: Vec<AlignmentMatch> = Vec::new();

    let ref_bases = read_fasta(fasta_path, chrom, from, to);

    for s in snippets {
        let mut cigar_offset: i64 = 0;
        let mut read_offset: i64 = 0;
        let base_string = s.sequence.clone();
        let char_vec: Vec<char> = base_string.chars().collect();

        let mut soft_clip_begin = true;

        let p = s.clone();

        if p.paired && (p.pos + p.length as i64) < p.mate_pos && p.tid == p.mate_tid {
            let pairing = AlignmentMatch {
                marker_type: Marker::Pairing,
                start_position: (p.pos + p.length as i64) as f64 - 0.5,
                end_position: p.mate_pos as f64 - 0.5,
                flags: p.flags.clone(),
                name: p.name.clone(),
                read_start: p.pos.clone() as u32,
                read_end: (p.mate_pos.clone() + 100) as u32,
            };

            matches.push(pairing);
        }

        for c in s.cigar.iter() {
            let mut match_count = 0;
            let mut match_start = 0;
            let mut match_ending = false;

            match c {
                rust_htslib::bam::record::Cigar::Match(c) => {
                    for _i in 0..rust_htslib::bam::record::Cigar::Match(*c).len() {
                        let snip = s.clone();
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
                                match mtch {
                                    Some(m) => matches.push(m),
                                    _ => {}
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
                        let mtch = end_mismatch_detection(s.clone(), match_start, match_count);
                        matches.push(mtch);
                    }

                    soft_clip_begin = false;
                }
                rust_htslib::bam::record::Cigar::Ins(c) => {
                    let snip = s.clone();
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
                        start_position: p.clone() as f64 - 0.5,
                        end_position: p as f64 + 0.5,
                        flags: snip.flags,
                        name: snip.name,
                        read_start: rs as u32,
                        read_end: re as u32,
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
                        let snip = s.clone();
                        let m = Marker::Deletion;
                        let p = snip.pos as i64 + read_offset;
                        let f = snip.flags;
                        let n = snip.name;
                        let rs;
                        let re;
                        let b = String::from("");

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
                            start_position: p.clone() as f64 - 0.5,
                            end_position: p as f64 + 0.5,
                            flags: f,
                            name: n,
                            read_start: rs as u32,
                            read_end: re as u32,
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
                        let snip = s.clone();
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
                                match mtch {
                                    Some(m) => matches.push(m),
                                    _ => {}
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
                        let mtch = end_mismatch_detection(s.clone(), match_start, match_count);
                        matches.push(mtch);
                    }

                    soft_clip_begin = false;
                }
                rust_htslib::bam::record::Cigar::HardClip(c) => {
                    for _i in 0..rust_htslib::bam::record::Cigar::HardClip(*c).len() {
                        cigar_offset += 1;
                    }

                    soft_clip_begin = false;
                }
                _ => {
                    soft_clip_begin = false;
                }
            }
        }
    }
    (bases, matches)
}

fn make_markers(
    snip: Alignment,
    b: char,
    read_offset: i64,
    match_start: i64,
    match_count: i64,
) -> (Option<AlignmentMatch>, AlignmentNucleobase) {
    let m: Marker;

    match b {
        // Mismatch
        'A' => m = Marker::A,
        'T' => m = Marker::T,
        'C' => m = Marker::C,
        'N' => m = Marker::N,
        'G' => m = Marker::G,
        _ => m = Marker::Deletion,
    }

    let p = snip.pos as i64 + read_offset;
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

    let mut mtch = None;

    if match_count > 0 {
        // First mismatch detection must lead to new creation of all previous matches
        mtch = Some(AlignmentMatch {
            marker_type: Marker::Match,
            start_position: match_start as f64 - 0.5,
            end_position: (match_start + match_count - 1) as f64 + 0.5,
            flags: f.clone(),
            name: n.clone(),
            read_start: rs.clone() as u32,
            read_end: re.clone() as u32,
        });
    }

    let base = AlignmentNucleobase {
        marker_type: m,
        bases: b.to_string(),
        start_position: p.clone() as f64 - 0.5,
        end_position: p as f64 + 0.5,
        flags: f,
        name: n,
        read_start: rs as u32,
        read_end: re as u32,
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

    let mtch = AlignmentMatch {
        marker_type: Marker::Match,
        start_position: match_start as f64 - 0.5,
        end_position: (match_start + match_count - 1) as f64 + 0.5,
        flags: f.clone(),
        name: n.clone(),
        read_start: rs.clone() as u32,
        read_end: re.clone() as u32,
    };

    mtch
}
