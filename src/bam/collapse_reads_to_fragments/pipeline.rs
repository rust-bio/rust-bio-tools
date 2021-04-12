use super::calc_consensus::{CalcNonOverlappingConsensus, CalcOverlappingConsensus};
use bio::io::fastq;
use derive_new::new;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::io;
use std::ops::Deref;
use uuid::Uuid;

#[derive(new)]
pub struct CallConsensusRead<W: io::Write> {
    bam_reader: bam::Reader,
    fq1_writer: fastq::Writer<W>,
    fq2_writer: fastq::Writer<W>,
    fq_se_writer: fastq::Writer<W>,
    bam_skipped_writer: bam::Writer,
    verbose_read_names: bool,
}

type Position = i64;
type GroupIDs = HashSet<GroupID>;
type RecordIDs = Vec<RecordID>;

#[derive(Hash, PartialEq, Eq)]
pub enum RecordID {
    Regular(Vec<u8>),
    Splitted(Vec<u8>),
}

#[derive(Hash, PartialEq, Eq, Clone)]
pub enum GroupID {
    Regular(i64),
    Splitted(i64),
}

impl<W: io::Write> CallConsensusRead<W> {
    pub fn call_consensus_reads(&mut self) -> Result<(), Box<dyn Error>> {
        let mut group_end_idx: BTreeMap<Position, GroupIDs> = BTreeMap::new();
        let mut duplicate_groups: HashMap<GroupID, RecordIDs> = HashMap::new();
        let mut record_storage: HashMap<RecordID, RecordStorage> = HashMap::new();

        for (i, result) in self.bam_reader.records().enumerate() {
            let mut record = result?;
            if !record.is_unmapped() {
                //Process completed duplicate groups
                calc_consensus_complete_groups(
                    &mut group_end_idx,
                    &mut duplicate_groups,
                    Some(&record.pos()),
                    &mut record_storage,
                    &mut self.fq1_writer,
                    &mut self.fq2_writer,
                    &mut self.fq_se_writer,
                    &mut self.bam_skipped_writer,
                    self.verbose_read_names,
                )?;
                group_end_idx = group_end_idx.split_off(&record.pos()); //Remove processed indexes
            }
            if record.is_unmapped() || record.is_mate_unmapped() || (record.tid() != record.mtid())
            {
                self.bam_skipped_writer.write(&record)?;
                continue;
            }
            if record.is_supplementary() {
                //TODO Supplementary Alignment
                continue;
            }
            record.cache_cigar();
            let duplicate_id_option = record.aux(b"DI");
            let record_id = record.qname();
            //Check if record has duplicate ID
            match duplicate_id_option {
                //Case: duplicate ID exists
                Some(duplicate_id) => {
                    match record_storage.get_mut(&RecordID::Regular(record_id.to_owned())) {
                        //Case: Right record
                        Some(storage_entry) => {
                            //For right record save end position and duplicate group ID
                            let record_end_pos = record.cigar_cached().unwrap().end_pos() - 1;
                            let group_id = match storage_entry {
                                RecordStorage::PairedRecords { ref mut r2_rec, .. } => {
                                    let group_id = duplicate_id.integer();
                                    r2_rec.get_or_insert(IndexedRecord {
                                        rec: record,
                                        rec_id: i,
                                    });
                                    GroupID::Regular(group_id)
                                }
                                // This arm is reached if a mate is mapped to another chromosome.
                                // In that case a new duplicate and record ID is required
                                //TODO Handle split reads correctly (see TODO below)
                                RecordStorage::SingleRecord { .. } => {
                                    let group_id = duplicate_id.integer();
                                    duplicate_groups
                                        .entry(GroupID::Splitted(group_id))
                                        .or_insert_with(Vec::new)
                                        .push(RecordID::Splitted(record_id.to_owned()));
                                    record_storage.insert(
                                        RecordID::Splitted(record_id.to_owned()),
                                        RecordStorage::SingleRecord {
                                            rec: IndexedRecord {
                                                rec: record,
                                                rec_id: i,
                                            },
                                        },
                                    );
                                    GroupID::Splitted(group_id)
                                }
                            };
                            group_end_idx
                                .entry(record_end_pos)
                                .or_insert_with(HashSet::new)
                                .insert(group_id);
                        }
                        //Case: Left record or record w/o mate
                        None => {
                            //TODO Records mapped to different chromosomes should be saved as single record with a Spliited recordID
                            duplicate_groups
                                .entry(GroupID::Regular(duplicate_id.integer()))
                                .or_insert_with(Vec::new)
                                .push(RecordID::Regular(record_id.to_owned()));
                            if !record.is_paired() {
                                //If right or single record save end position and duplicate group ID
                                group_end_idx
                                    .entry(record.cigar_cached().unwrap().end_pos() - 1)
                                    .or_insert_with(HashSet::new)
                                    .insert(GroupID::Regular(duplicate_id.integer()));
                                record_storage.insert(
                                    RecordID::Regular(record_id.to_owned()),
                                    RecordStorage::SingleRecord {
                                        rec: IndexedRecord {
                                            rec: record,
                                            rec_id: i,
                                        },
                                    },
                                );
                            } else {
                                record_storage.insert(
                                    RecordID::Regular(record_id.to_owned()),
                                    RecordStorage::PairedRecords {
                                        r1_rec: IndexedRecord {
                                            rec: record,
                                            rec_id: i,
                                        },
                                        r2_rec: None,
                                    },
                                );
                            }
                        }
                    }
                }
                //Duplicate ID not existing
                //Record is writen to bam file if it or its mate is unmapped
                //If record is right mate consensus is calculated
                //Else record is added to hashMap
                None => {
                    match record_storage.get_mut(&RecordID::Regular(record_id.to_owned())) {
                        //Case: Left record
                        None => {
                            record_storage.insert(
                                RecordID::Regular(record_id.to_owned()),
                                RecordStorage::PairedRecords {
                                    r1_rec: IndexedRecord {
                                        rec: record,
                                        rec_id: i,
                                    },
                                    r2_rec: None,
                                },
                            );
                        }
                        //Case: Left record already stored
                        Some(_record_pair) => {
                            let (rec_id, l_rec) = match record_storage
                                .remove(&RecordID::Regular(record_id.to_owned()))
                                .unwrap()
                            {
                                RecordStorage::PairedRecords { r1_rec, .. } => {
                                    (r1_rec.rec_id, r1_rec.into_rec())
                                }
                                RecordStorage::SingleRecord { .. } => unreachable!(),
                            };
                            if cigar_has_softclips(&l_rec) || cigar_has_softclips(&record) {
                                self.bam_skipped_writer.write(&l_rec)?;
                                self.bam_skipped_writer.write(&record)?;
                            } else {
                                let alignment_vectors = calc_read_alignments(&l_rec, &record);
                                match alignment_vectors {
                                    Some((r1_alignment, r2_alignment)) => {
                                        let uuid = &Uuid::new_v4().to_hyphenated().to_string();

                                        self.fq_se_writer.write_record(
                                            &CalcOverlappingConsensus::new(
                                                &[l_rec],
                                                &[record],
                                                &r1_alignment,
                                                &r2_alignment,
                                                &[rec_id, i],
                                                uuid,
                                                self.verbose_read_names,
                                            )
                                            .calc_consensus()
                                            .0,
                                        )?;
                                    }
                                    None => {
                                        self.bam_skipped_writer.write(&l_rec)?;
                                        self.bam_skipped_writer.write(&record)?;
                                    }
                                };
                            }
                        }
                    }
                }
            }
        }
        //Process remaining groups
        calc_consensus_complete_groups(
            &mut group_end_idx,
            &mut duplicate_groups,
            None,
            &mut record_storage,
            &mut self.fq1_writer,
            &mut self.fq2_writer,
            &mut self.fq_se_writer,
            &mut self.bam_skipped_writer,
            self.verbose_read_names,
        )?;
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
pub fn calc_consensus_complete_groups<'a, W: io::Write>(
    group_end_idx: &mut BTreeMap<Position, GroupIDs>,
    duplicate_groups: &mut HashMap<GroupID, RecordIDs>,
    end_pos: Option<&i64>,
    record_storage: &mut HashMap<RecordID, RecordStorage>,
    fq1_writer: &'a mut fastq::Writer<W>,
    fq2_writer: &'a mut fastq::Writer<W>,
    fq_se_writer: &'a mut fastq::Writer<W>,
    bam_skipped_writer: &'a mut bam::Writer,
    verbose_read_names: bool,
) -> Result<(), Box<dyn Error>> {
    let group_ids: HashSet<GroupID> = group_end_idx
        .range(
            ..end_pos.unwrap_or(
                &(group_end_idx
                    .iter()
                    .next_back()
                    .map_or(0, |(entry, _)| *entry)
                    + 1),
            ),
        )
        .flat_map(|(_, group_ids)| group_ids.clone())
        .collect();
    for group_id in group_ids {
        let cigar_groups = group_reads_by_cigar(
            duplicate_groups.remove(&group_id).unwrap(),
            record_storage,
            bam_skipped_writer,
        )?;
        for cigar_group in cigar_groups.values() {
            match cigar_group {
                CigarGroup::PairedRecords {
                    r1_recs,
                    r2_recs,
                    r1_seqids,
                    r2_seqids,
                } => {
                    let alignment_vectors = calc_read_alignments(&r1_recs[0], &r2_recs[0]);
                    match alignment_vectors {
                        Some((r1_alignment, r2_alignment)) => {
                            let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                            let mut seqids = r1_seqids.clone();
                            seqids.append(&mut r2_seqids.clone());
                            fq_se_writer.write_record(
                                &CalcOverlappingConsensus::new(
                                    &r1_recs,
                                    &r2_recs,
                                    &r1_alignment,
                                    &r2_alignment,
                                    &seqids,
                                    uuid,
                                    verbose_read_names,
                                )
                                .calc_consensus()
                                .0,
                            )?;
                        }
                        None => {
                            // If reads do not overlap or CIGAR in overlapping region differs R1 and R2 are handled sepperatly
                            if r1_recs.len() > 1 {
                                let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                                fq1_writer.write_record(
                                    &CalcNonOverlappingConsensus::new(&r1_recs, &r1_seqids, uuid)
                                        .calc_consensus()
                                        .0,
                                )?;
                                fq2_writer.write_record(
                                    &CalcNonOverlappingConsensus::new(&r2_recs, &r2_seqids, uuid)
                                        .calc_consensus()
                                        .0,
                                )?;
                            } else {
                                bam_skipped_writer.write(&r1_recs[0])?;
                                bam_skipped_writer.write(&r2_recs[0])?;
                            }
                        }
                    };
                }
                CigarGroup::SingleRecords { recs, seqids } => match recs.len().cmp(&1) {
                    Ordering::Greater => {
                        let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                        fq_se_writer.write_record(
                            &CalcNonOverlappingConsensus::new(&recs, &seqids, uuid)
                                .calc_consensus()
                                .0,
                        )?;
                    }
                    _ => {
                        bam_skipped_writer.write(&recs[0])?;
                    }
                },
            }
        }
    }
    Ok(())
}

fn group_reads_by_cigar(
    record_ids: Vec<RecordID>,
    record_storage: &mut HashMap<RecordID, RecordStorage>,
    bam_skipped_writer: &mut bam::Writer,
) -> Result<HashMap<Cigar, CigarGroup>, Box<dyn Error>> {
    let mut cigar_groups: HashMap<Cigar, CigarGroup> = HashMap::new();
    for rec_id in record_ids {
        match record_storage.remove(&rec_id).unwrap() {
            RecordStorage::PairedRecords { r1_rec, r2_rec } => {
                let r1_rec_id = r1_rec.rec_id;
                let r1_rec_entry = r1_rec.into_rec();
                let r2_rec_unwrapped = r2_rec.unwrap();
                let r2_rec_id = r2_rec_unwrapped.rec_id;
                let r2_rec_entry = r2_rec_unwrapped.into_rec();

                if cigar_has_softclips(&r1_rec_entry) || cigar_has_softclips(&r2_rec_entry) {
                    bam_skipped_writer.write(&r1_rec_entry)?;
                    bam_skipped_writer.write(&r2_rec_entry)?;
                } else {
                    let cigar_tuple = Cigar::Tuple {
                        r1_cigar: r1_rec_entry.raw_cigar().to_vec(),
                        r2_cigar: r2_rec_entry.raw_cigar().to_vec(),
                    };
                    if !cigar_groups.contains_key(&cigar_tuple) {
                        cigar_groups.insert(
                            cigar_tuple.clone(),
                            CigarGroup::PairedRecords {
                                r1_recs: Vec::new(),
                                r2_recs: Vec::new(),
                                r1_seqids: Vec::new(),
                                r2_seqids: Vec::new(),
                            },
                        );
                    }
                    match cigar_groups.get_mut(&cigar_tuple) {
                        Some(CigarGroup::PairedRecords {
                            r1_recs,
                            r2_recs,
                            r1_seqids,
                            r2_seqids,
                        }) => {
                            r1_recs.push(r1_rec_entry);
                            r2_recs.push(r2_rec_entry);
                            r1_seqids.push(r1_rec_id);
                            r2_seqids.push(r2_rec_id);
                        }
                        _ => unreachable!(),
                    }
                }
            }
            RecordStorage::SingleRecord { rec } => {
                let rec_id = rec.rec_id;
                let rec_entry = rec.into_rec();
                if cigar_has_softclips(&rec_entry) {
                    bam_skipped_writer.write(&rec_entry)?;
                } else {
                    let cigar_single = Cigar::Single {
                        cigar: rec_entry.raw_cigar().to_vec(),
                    };
                    if !cigar_groups.contains_key(&cigar_single) {
                        cigar_groups.insert(
                            cigar_single.clone(),
                            CigarGroup::SingleRecords {
                                recs: Vec::new(),
                                seqids: Vec::new(),
                            },
                        );
                    }
                    match cigar_groups.get_mut(&cigar_single) {
                        Some(CigarGroup::SingleRecords { recs, seqids }) => {
                            recs.push(rec_entry);
                            seqids.push(rec_id);
                        }
                        _ => unreachable!(),
                    }
                }
            }
        };
    }
    Ok(cigar_groups)
}

fn calc_read_alignments(
    r1_rec: &bam::Record,
    r2_rec: &bam::Record,
) -> Option<(Vec<bool>, Vec<bool>)> {
    let r1_start = r1_rec.pos();
    let r1_end = r1_rec.cigar_cached().unwrap().end_pos();
    let r2_start = r2_rec.pos();
    let r2_end = r1_rec.cigar_cached().unwrap().end_pos();

    if r1_start <= r2_start {
        //Check if reads overlap
        if r1_end >= r2_start {
            let offset = r2_start - r1_start;
            calc_alignment_vectors(offset, r1_rec, r2_rec)
        } else {
            //Reads do not overlap
            None
        }
    } else {
        //R2 starts before R1
        if r2_end >= r1_start {
            let offset = r1_start - r2_start;
            calc_alignment_vectors(offset, r2_rec, r1_rec)
        } else {
            None
        }
    }
}

fn calc_alignment_vectors(
    mut offset: i64,
    r1_rec: &bam::Record,
    r2_rec: &bam::Record,
) -> Option<(Vec<bool>, Vec<bool>)> {
    let mut r1_vec = Vec::new();
    let mut r2_vec = Vec::new();
    let mut r1_cigarstring = r1_rec
        .cigar_cached()
        .unwrap()
        .iter()
        .flat_map(|cigar| vec![cigar.char(); cigar.len() as usize])
        .collect::<Vec<char>>()
        .into_iter();
    let mut r2_cigarstring = r2_rec
        .cigar_cached()
        .unwrap()
        .iter()
        .flat_map(|cigar| vec![cigar.char(); cigar.len() as usize])
        .collect::<Vec<char>>()
        .into_iter();
    let mut r1_cigar = r1_cigarstring.next();
    let mut r2_cigar = match offset == 0 {
        true => r2_cigarstring.next(),
        false => None,
    };
    loop {
        if r2_cigar == None {
            match r1_cigar {
                None => break,
                Some('M') | Some('X') | Some('=') | Some('D') | Some('N') => {
                    offset -= 1;
                }
                Some('S') => unreachable!(),
                Some(_) => {}
            }
            match_single_cigar(&r1_cigar, &mut r1_vec, &mut r2_vec);
            r1_cigar = r1_cigarstring.next();
            if offset == 0 {
                r2_cigar = r2_cigarstring.next();
            }
        } else if r1_cigar == None {
            match_single_cigar(&r2_cigar, &mut r2_vec, &mut r1_vec);
            r2_cigar = r2_cigarstring.next();
        } else {
            if r1_cigar != r2_cigar {
                return None;
            }
            match (r1_cigar, r2_cigar) {
                (Some('M'), Some('M'))
                | (Some('X'), Some('X'))
                | (Some('='), Some('='))
                | (Some('I'), Some('I')) => {
                    r1_vec.push(true);
                    r2_vec.push(true);
                    r1_cigar = r1_cigarstring.next();
                    r2_cigar = r2_cigarstring.next();
                }
                (Some('D'), Some('D')) | (Some('H'), Some('H')) => {
                    r1_cigar = r1_cigarstring.next();
                    r2_cigar = r2_cigarstring.next();
                }
                (None, None) | (None, Some(_)) | (Some(_), None) | (Some(_), Some(_)) => {
                    unreachable!()
                }
            };
        }
    }
    Some((r1_vec, r2_vec))
}

/* fn calc_alignment_vectors(
    mut offset: i64,
    r1_rec: &bam::Record,
    r2_rec: &bam::Record,
) -> (Vec<bool>, Vec<bool>) {
    let mut r1_vec = Vec::new();
    let mut r2_vec = Vec::new();
    let mut r1_cigarstring = r1_rec
        .cigar_cached()
        .unwrap()
        .iter()
        .flat_map(|cigar| vec![cigar.char(); cigar.len() as usize])
        .collect::<Vec<char>>()
        .into_iter();
    let mut r2_cigarstring = r2_rec
        .cigar_cached()
        .unwrap()
        .iter()
        .flat_map(|cigar| vec![cigar.char(); cigar.len() as usize])
        .collect::<Vec<char>>()
        .into_iter();
    let mut r1_cigar = r1_cigarstring.next();
    let mut r2_cigar = match offset == 0 {
        true => r2_cigarstring.next(),
        false => None,
    };
    loop {
        if r2_cigar == None {
            match r1_cigar {
                None => break,
                Some('M') | Some('S') | Some('X') | Some('=') | Some('D') | Some('N') => {
                    offset -= 1;
                }
                Some(_) => {}
            }
            match_single_cigar(&r1_cigar, &mut r1_vec, &mut r2_vec);
            r1_cigar = r1_cigarstring.next();
            if offset == 0 {
                r2_cigar = r2_cigarstring.next();
            }
        } else if r1_cigar == None {
            match_single_cigar(&r2_cigar, &mut r2_vec, &mut r1_vec);
            r2_cigar = r2_cigarstring.next();
        } else {
            match (r1_cigar, r2_cigar) {
                (Some('M'), Some('M'))
                | (Some('S'), Some('S'))
                | (Some('X'), Some('X'))
                | (Some('='), Some('='))
                | (Some('I'), Some('I'))
                | (Some('M'), Some('S'))
                | (Some('M'), Some('X'))
                | (Some('M'), Some('='))
                | (Some('S'), Some('M'))
                | (Some('S'), Some('X'))
                | (Some('S'), Some('='))
                | (Some('X'), Some('M'))
                | (Some('X'), Some('S'))
                | (Some('X'), Some('='))
                | (Some('='), Some('M'))
                | (Some('='), Some('S'))
                | (Some('='), Some('X')) => {
                    r1_vec.push(true);
                    r2_vec.push(true);
                    r1_cigar = r1_cigarstring.next();
                    r2_cigar = r2_cigarstring.next();
                }
                (Some('D'), Some('D'))
                | (Some('H'), Some('H'))
                | (Some('N'), Some('N'))
                | (Some('P'), Some('P'))
                | (Some('D'), Some('H'))
                | (Some('D'), Some('N'))
                | (Some('D'), Some('P'))
                | (Some('H'), Some('D'))
                | (Some('H'), Some('N'))
                | (Some('H'), Some('P'))
                | (Some('N'), Some('D'))
                | (Some('N'), Some('P'))
                | (Some('N'), Some('H'))
                | (Some('P'), Some('D'))
                | (Some('P'), Some('N'))
                | (Some('P'), Some('H')) => {
                    r1_cigar = r1_cigarstring.next();
                    r2_cigar = r2_cigarstring.next();
                }
                (Some('I'), Some(_)) => {
                    r1_vec.push(true);
                    r1_cigar = r1_cigarstring.next();
                    r2_vec.push(false);
                }
                (Some(_), Some('I')) => {
                    r1_vec.push(false);
                    r2_vec.push(true);
                    r2_cigar = r2_cigarstring.next();
                }
                (Some('M'), Some('D'))
                | (Some('M'), Some('N'))
                | (Some('M'), Some('H'))
                | (Some('M'), Some('P'))
                | (Some('S'), Some('D'))
                | (Some('S'), Some('N'))
                | (Some('S'), Some('H'))
                | (Some('S'), Some('P'))
                | (Some('X'), Some('D'))
                | (Some('X'), Some('N'))
                | (Some('X'), Some('H'))
                | (Some('X'), Some('P'))
                | (Some('='), Some('D'))
                | (Some('='), Some('N'))
                | (Some('='), Some('H'))
                | (Some('='), Some('P')) => {
                    r1_vec.push(true);
                    r2_vec.push(false);
                    r1_cigar = r1_cigarstring.next();
                    r2_cigar = r2_cigarstring.next();
                }
                (Some('D'), Some('M'))
                | (Some('N'), Some('M'))
                | (Some('H'), Some('M'))
                | (Some('P'), Some('M'))
                | (Some('D'), Some('S'))
                | (Some('N'), Some('S'))
                | (Some('H'), Some('S'))
                | (Some('P'), Some('S'))
                | (Some('D'), Some('X'))
                | (Some('N'), Some('X'))
                | (Some('H'), Some('X'))
                | (Some('P'), Some('X'))
                | (Some('D'), Some('='))
                | (Some('N'), Some('='))
                | (Some('H'), Some('='))
                | (Some('P'), Some('=')) => {
                    r1_vec.push(false);
                    r2_vec.push(true);
                    r1_cigar = r1_cigarstring.next();
                    r2_cigar = r2_cigarstring.next();
                }
                (None, None) | (None, Some(_)) | (Some(_), Some(_)) | (Some(_), None) => {
                    unreachable!()
                }
            };
        }
    }
    (r1_vec, r2_vec)
} */

fn cigar_has_softclips(rec: &bam::Record) -> bool {
    let cigar = rec.cigar_cached().unwrap();
    cigar.leading_softclips() > 0 || cigar.trailing_softclips() > 0
}

fn match_single_cigar(cigar: &Option<char>, first_vec: &mut Vec<bool>, second_vec: &mut Vec<bool>) {
    match cigar {
        Some('M') | Some('S') | Some('X') | Some('=') | Some('I') => {
            first_vec.push(true);
            second_vec.push(false);
        }
        Some(_) | None => {}
    };
}
pub enum RecordStorage {
    PairedRecords {
        r1_rec: IndexedRecord,
        r2_rec: Option<IndexedRecord>,
    },
    SingleRecord {
        rec: IndexedRecord,
    },
}

pub struct IndexedRecord {
    rec: bam::Record,
    rec_id: usize,
}

impl IndexedRecord {
    fn into_rec(self) -> bam::Record {
        self.rec
    }
}

impl Deref for IndexedRecord {
    type Target = bam::Record;
    fn deref(&self) -> &bam::Record {
        &self.rec
    }
}

pub enum CigarGroup {
    PairedRecords {
        r1_recs: Vec<bam::Record>,
        r2_recs: Vec<bam::Record>,
        r1_seqids: Vec<usize>,
        r2_seqids: Vec<usize>,
    },
    SingleRecords {
        recs: Vec<bam::Record>,
        seqids: Vec<usize>,
    },
}

#[derive(Hash, PartialEq, Eq, Clone)]
pub enum Cigar {
    Tuple {
        r1_cigar: Vec<u32>,
        r2_cigar: Vec<u32>,
    },
    Single {
        cigar: Vec<u32>,
    },
}
