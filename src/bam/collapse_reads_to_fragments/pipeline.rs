use super::calc_consensus::{CalcNonOverlappingConsensus, CalcOverlappingConsensus};
use super::unmark_record;
use anyhow::Result;
use bio::io::fastq;
use derive_new::new;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Read;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io;
use std::ops::Deref;
use std::ops::DerefMut;
use uuid::Uuid;

#[derive(new)]
pub struct CallConsensusRead<W: io::Write> {
    bam_reader: bam::Reader,
    fq1_writer: fastq::Writer<W>,
    fq2_writer: fastq::Writer<W>,
    fq_se_writer: fastq::Writer<W>,
    bam_skipped_writer: bam::Writer,
    annotate_record_ids: bool,
}

type Position = i64;
type GroupIDs = HashSet<GroupId>;
type RecordIDs = Vec<RecordId>;

#[derive(Hash, PartialEq, Eq, Debug)]
pub enum RecordId {
    Regular(Vec<u8>),
    Split(Vec<u8>),
}

#[derive(Hash, PartialEq, Eq, Clone, Debug)]
pub enum GroupId {
    Regular(u32),
    Split(u32),
}

#[derive(new, Debug)]
pub struct GroupEndIndex {
    #[new(default)]
    group_pos: HashMap<GroupId, Position>,
    #[new(default)]
    group_end_idx: BTreeMap<Position, GroupIDs>,
}

impl GroupEndIndex {
    ///Inserts a new group id at given position
    ///If position is already saved for the group id the group-end-index will be updated
    pub fn insert(&mut self, group_id: GroupId, end_pos: i64) -> Result<()> {
        let update_end_pos = match self.group_pos.get(&group_id) {
            Some(&current_end_pos) => match current_end_pos < end_pos {
                true => {
                    self.group_end_idx
                        .get_mut(&current_end_pos)
                        .map(|group_ids| group_ids.remove(&group_id));
                    true
                }
                false => false,
            },
            None => true,
        };
        if update_end_pos {
            self.group_pos.insert(group_id.clone(), end_pos);
            self.group_end_idx
                .entry(end_pos)
                .or_insert_with(HashSet::new)
                .insert(group_id);
        }
        Ok(())
    }

    pub fn cut_lower_group_ids(&mut self, current_pos: Option<i64>) -> Result<Vec<GroupId>> {
        let group_ids: Vec<GroupId> = self
            .group_end_idx
            .range(
                ..current_pos.unwrap_or(
                    self.group_end_idx
                        .iter()
                        .next_back()
                        .map_or(0, |(entry, _)| *entry)
                        + 1,
                ),
            )
            .flat_map(|(_, group_ids)| group_ids.clone())
            .collect();
        group_ids.iter().for_each(|group_id| {
            self.group_pos.remove(group_id);
        });
        match current_pos {
            Some(pos) => self.group_end_idx = self.group_end_idx.split_off(&pos),
            None => self.group_end_idx.clear(),
        }
        Ok(group_ids)
    }
}

impl<W: io::Write> CallConsensusRead<W> {
    pub fn call_consensus_reads(&mut self) -> Result<()> {
        let mut group_end_idx = GroupEndIndex::new();
        let mut duplicate_groups: HashMap<GroupId, RecordIDs> = HashMap::new();
        let mut record_storage: HashMap<RecordId, RecordStorage> = HashMap::new();
        let mut current_chrom = None;
        let mut read_ids: Option<HashMap<usize, Vec<u8>>> = if self.annotate_record_ids {
            Some(HashMap::new())
        } else {
            None
        };
        for (i, result) in self.bam_reader.records().enumerate() {
            let mut record = result?;
            if !record.is_unmapped() {
                let mut record_pos = None;
                match current_chrom == Some(record.tid()) {
                    true => record_pos = Some(record.pos()),
                    false => current_chrom = Some(record.tid()),
                }
                //Process completed duplicate groups
                calc_consensus_complete_groups(
                    &mut group_end_idx,
                    &mut duplicate_groups,
                    record_pos,
                    &mut record_storage,
                    &mut self.fq1_writer,
                    &mut self.fq2_writer,
                    &mut self.fq_se_writer,
                    &mut self.bam_skipped_writer,
                    &mut read_ids,
                )?;
            }
            if record.is_unmapped() || record.is_mate_unmapped() {
                unmark_record(&mut record)?;
                self.bam_skipped_writer.write(&record)?;
                continue;
            }
            if record.is_supplementary() {
                //TODO Supplementary Alignment
                continue;
            }
            record.cache_cigar();
            let duplicate_id_option = match record.aux(b"DI") {
                Ok(Aux::I8(duplicate_id)) => Some(duplicate_id as u32),
                Ok(Aux::I16(duplicate_id)) => Some(duplicate_id as u32),
                Ok(Aux::I32(duplicate_id)) => Some(duplicate_id as u32),
                Ok(Aux::U8(duplicate_id)) => Some(duplicate_id as u32),
                Ok(Aux::U16(duplicate_id)) => Some(duplicate_id as u32),
                Ok(Aux::U32(duplicate_id)) => Some(duplicate_id),
                Err(_) => None,
                _ => unreachable!("Invalid type for tag 'DI'"),
            };
            let record_name = record.qname();
            read_ids.as_mut().map(|x| x.insert(i, record_name.to_vec()));
            //Check if record has duplicate ID
            match duplicate_id_option {
                //Case: duplicate ID exists
                Some(duplicate_id) => {
                    let regular_id = RecordId::Regular(record_name.to_owned());
                    let record_end_pos = record.cigar_cached().unwrap().end_pos() - 1;
                    match record_storage.get_mut(&regular_id) {
                        //Case: Right record
                        Some(storage_entry) => {
                            //For right record save end position and duplicate group ID
                            let group_id_opt = match storage_entry {
                                RecordStorage::PairedRecords {
                                    ref mut r1_rec,
                                    ref mut r2_rec,
                                } => {
                                    let group_id = if cigar_has_softclips(r1_rec)
                                        || cigar_has_softclips(&record)
                                    {
                                        unmark_record(r1_rec)?;
                                        self.bam_skipped_writer.write(r1_rec)?;
                                        unmark_record(&mut record)?;
                                        self.bam_skipped_writer.write(&record)?;
                                        None
                                    } else {
                                        duplicate_groups
                                            .entry(GroupId::Regular(duplicate_id))
                                            .or_insert_with(Vec::new)
                                            .push(RecordId::Regular(record_name.to_owned()));
                                        r2_rec.get_or_insert(IndexedRecord {
                                            rec: record,
                                            rec_id: i,
                                        });
                                        Some(GroupId::Regular(duplicate_id))
                                    };
                                    group_id
                                }
                                // This arm is reached if a mate is mapped to another chromosome.
                                // In that case a new duplicate and record ID is required
                                RecordStorage::SingleRecord { rec } => {
                                    let group_id = if cigar_has_softclips(rec)
                                        || cigar_has_softclips(&record)
                                    {
                                        unmark_record(rec)?;
                                        self.bam_skipped_writer.write(rec)?;
                                        unmark_record(&mut record)?;
                                        self.bam_skipped_writer.write(&record)?;
                                        None
                                    } else {
                                        duplicate_groups
                                            .entry(GroupId::Split(duplicate_id))
                                            .or_insert_with(Vec::new)
                                            .push(RecordId::Split(record_name.to_owned()));
                                        record_storage.insert(
                                            RecordId::Split(record_name.to_owned()),
                                            RecordStorage::SingleRecord {
                                                rec: IndexedRecord {
                                                    rec: record,
                                                    rec_id: i,
                                                },
                                            },
                                        );
                                        Some(GroupId::Split(duplicate_id))
                                    };
                                    group_id
                                }
                            };
                            if let Some(group_id) = group_id_opt {
                                group_end_idx.insert(group_id, record_end_pos)?;
                            } else {
                                record_storage.remove(&regular_id);
                            };
                        }
                        //Case: Left record or record w/o mate
                        None => {
                            if !record.is_paired() {
                                //If right or single record save end position and duplicate group ID
                                if cigar_has_softclips(&record) {
                                    unmark_record(&mut record)?;
                                    self.bam_skipped_writer.write(&record)?;
                                } else {
                                    duplicate_groups
                                        .entry(GroupId::Regular(duplicate_id))
                                        .or_insert_with(Vec::new)
                                        .push(RecordId::Regular(record_name.to_owned()));

                                    group_end_idx
                                        .insert(GroupId::Regular(duplicate_id), record_end_pos)?;
                                    record_storage.insert(
                                        RecordId::Regular(record_name.to_owned()),
                                        RecordStorage::SingleRecord {
                                            rec: IndexedRecord {
                                                rec: record,
                                                rec_id: i,
                                            },
                                        },
                                    );
                                }
                            } else {
                                record_storage.insert(
                                    RecordId::Regular(record_name.to_owned()),
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
                //Record is written to bam file if it or its mate is unmapped
                //If record is right mate consensus is calculated
                //Else record is added to hashMap
                None => {
                    match record_storage.get_mut(&RecordId::Regular(record_name.to_owned())) {
                        //Case: Left record
                        None => {
                            if !record.is_paired() || record.tid() != record.mtid() {
                                unmark_record(&mut record)?;
                                self.bam_skipped_writer.write(&record)?;
                            } else {
                                record_storage.insert(
                                    RecordId::Regular(record_name.to_owned()),
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
                        //Case: Left record already stored
                        Some(_record_pair) => {
                            let (rec_id, mut l_rec) = match record_storage
                                .remove(&RecordId::Regular(record_name.to_owned()))
                                .unwrap()
                            {
                                RecordStorage::PairedRecords { r1_rec, .. } => {
                                    (r1_rec.rec_id, r1_rec.into_rec())
                                }
                                RecordStorage::SingleRecord { .. } => unreachable!(),
                            };
                            if cigar_has_softclips(&l_rec) || cigar_has_softclips(&record) {
                                unmark_record(&mut l_rec)?;
                                self.bam_skipped_writer.write(&l_rec)?;
                                unmark_record(&mut record)?;
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
                                                &mut read_ids,
                                            )
                                            .calc_consensus()
                                            .0,
                                        )?;
                                    }
                                    None => {
                                        unmark_record(&mut l_rec)?;
                                        self.bam_skipped_writer.write(&l_rec)?;
                                        unmark_record(&mut record)?;
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
            &mut read_ids,
        )?;
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
pub fn calc_consensus_complete_groups<'a, W: io::Write>(
    group_end_idx: &mut GroupEndIndex,
    duplicate_groups: &mut HashMap<GroupId, RecordIDs>,
    end_pos: Option<i64>,
    record_storage: &mut HashMap<RecordId, RecordStorage>,
    fq1_writer: &'a mut fastq::Writer<W>,
    fq2_writer: &'a mut fastq::Writer<W>,
    fq_se_writer: &'a mut fastq::Writer<W>,
    bam_skipped_writer: &'a mut bam::Writer,
    read_ids: &'a mut Option<HashMap<usize, Vec<u8>>>,
) -> Result<()> {
    let group_ids = group_end_idx.cut_lower_group_ids(end_pos)?;
    for group_id in group_ids {
        let cigar_groups =
            group_reads(duplicate_groups.remove(&group_id).unwrap(), record_storage)?;
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
                                    r1_recs,
                                    r2_recs,
                                    &r1_alignment,
                                    &r2_alignment,
                                    &seqids,
                                    uuid,
                                    read_ids,
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
                                    &CalcNonOverlappingConsensus::new(
                                        r1_recs, r1_seqids, uuid, read_ids,
                                    )
                                    .calc_consensus()
                                    .0,
                                )?;
                                fq2_writer.write_record(
                                    &CalcNonOverlappingConsensus::new(
                                        r2_recs, r2_seqids, uuid, read_ids,
                                    )
                                    .calc_consensus()
                                    .0,
                                )?;
                            } else {
                                let mut r1_rec = r1_recs[0].clone();
                                unmark_record(&mut r1_rec)?;
                                bam_skipped_writer.write(&r1_rec)?;
                                let mut r2_rec = r2_recs[0].clone();
                                unmark_record(&mut r2_rec)?;
                                bam_skipped_writer.write(&r2_rec)?;
                            }
                        }
                    };
                }
                CigarGroup::SingleRecords { recs, seqids } => match recs.len().cmp(&1) {
                    Ordering::Greater => {
                        let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                        fq_se_writer.write_record(
                            &CalcNonOverlappingConsensus::new(recs, seqids, uuid, read_ids)
                                .calc_consensus()
                                .0,
                        )?;
                    }
                    _ => {
                        let mut rec = recs[0].clone();
                        unmark_record(&mut rec)?;
                        bam_skipped_writer.write(&rec)?;
                    }
                },
            }
        }
    }
    Ok(())
}

// Final groups for consensus reads are created based on record and read orientations
fn group_reads(
    record_ids: Vec<RecordId>,
    record_storage: &mut HashMap<RecordId, RecordStorage>,
) -> Result<HashMap<(Cigar, String), CigarGroup>> {
    let mut final_groups: HashMap<(Cigar, String), CigarGroup> = HashMap::new();
    for rec_id in record_ids {
        let storage_entry = record_storage.remove(&rec_id).unwrap();
        storage_entry.add_to_group(&mut final_groups)?;
    }
    Ok(final_groups)
}
fn calc_read_alignments(
    r1_rec: &bam::Record,
    r2_rec: &bam::Record,
) -> Option<(Vec<bool>, Vec<bool>)> {
    let r1_start = r1_rec.pos();
    let r1_end = r1_rec.cigar_cached().unwrap().end_pos();
    let r2_start = r2_rec.pos();
    let r2_end = r1_rec.cigar_cached().unwrap().end_pos();
    if r1_rec.tid() != r2_rec.tid() {
        None
    } else if r1_start <= r2_start {
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
    let mut intersection_entry_passed = false;
    loop {
        if r2_cigar.is_none() {
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
        } else if r1_cigar.is_none() {
            match_single_cigar(&r2_cigar, &mut r2_vec, &mut r1_vec);
            r2_cigar = r2_cigarstring.next();
        } else if r1_cigar != r2_cigar {
            if !intersection_entry_passed && r1_cigar == Some('I') {
                r1_vec.push(true);
                r2_vec.push(false);
                r1_cigar = r1_cigarstring.next();
            } else {
                return None;
            }
        } else {
            intersection_entry_passed = true; // Can this me somehow only be called once?!
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

fn cigar_has_softclips(rec: &bam::Record) -> bool {
    for cigar_operation in rec.cigar_cached().unwrap().iter() {
        if let bam::record::Cigar::SoftClip(_) = cigar_operation {
            return true;
        }
    }
    false
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

impl RecordStorage {
    fn add_to_group(self, final_groups: &mut HashMap<(Cigar, String), CigarGroup>) -> Result<()> {
        let (r1_rec_entry, r1_rec_id, r2_rec_entry, r2_rec_id, cigar_tuple, orientation) =
            match self {
                RecordStorage::PairedRecords { r1_rec, r2_rec } => {
                    let r1_rec_id = r1_rec.rec_id;
                    let r1_rec_entry = r1_rec.into_rec();
                    let r2_rec_unwrapped = r2_rec.unwrap();
                    let r2_rec_id = r2_rec_unwrapped.rec_id;
                    let r2_rec_entry = r2_rec_unwrapped.into_rec();
                    let cigar_tuple = Cigar::Tuple {
                        r1_cigar: r1_rec_entry.raw_cigar().to_vec(),
                        r2_cigar: r2_rec_entry.raw_cigar().to_vec(),
                    };
                    let pair_orientation =
                        r1_rec_entry.read_pair_orientation().as_ref().to_string();
                    final_groups
                        .entry((cigar_tuple.clone(), pair_orientation.to_string()))
                        .or_insert_with(|| CigarGroup::PairedRecords {
                            r1_recs: Vec::new(),
                            r2_recs: Vec::new(),
                            r1_seqids: Vec::new(),
                            r2_seqids: Vec::new(),
                        });
                    (
                        r1_rec_entry,
                        r1_rec_id,
                        Some(r2_rec_entry),
                        Some(r2_rec_id),
                        cigar_tuple,
                        pair_orientation,
                    )
                }
                RecordStorage::SingleRecord { rec } => {
                    let rec_id = rec.rec_id;
                    let rec_entry = rec.into_rec();
                    let cigar_single = Cigar::Single {
                        cigar: rec_entry.raw_cigar().to_vec(),
                    };
                    let read_orientation = if rec_entry.is_reverse() {
                        "reverse"
                    } else {
                        "forward"
                    }
                    .to_string();
                    final_groups
                        .entry((cigar_single.clone(), read_orientation.clone()))
                        .or_insert_with(|| CigarGroup::SingleRecords {
                            recs: Vec::new(),
                            seqids: Vec::new(),
                        });
                    (
                        rec_entry,
                        rec_id,
                        None,
                        None,
                        cigar_single,
                        read_orientation,
                    )
                }
            };
        match final_groups.get_mut(&(cigar_tuple, orientation)) {
            Some(CigarGroup::PairedRecords {
                r1_recs,
                r2_recs,
                r1_seqids,
                r2_seqids,
            }) => {
                r1_recs.push(r1_rec_entry);
                r2_recs.push(r2_rec_entry.unwrap());
                r1_seqids.push(r1_rec_id);
                r2_seqids.push(r2_rec_id.unwrap());
            }
            Some(CigarGroup::SingleRecords { recs, seqids }) => {
                recs.push(r1_rec_entry);
                seqids.push(r1_rec_id);
            }
            None => unreachable!(),
        }

        Ok(())
    }
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

impl DerefMut for IndexedRecord {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.rec
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
