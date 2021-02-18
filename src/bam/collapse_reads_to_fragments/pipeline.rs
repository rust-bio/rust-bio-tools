use super::calc_consensus::{CalcNonOverlappingConsensus, CalcOverlappingConsensus};
use bio::io::fastq;
use derive_new::new;
use rust_htslib::bam;
use rust_htslib::bam::Read;
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
            } else {
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
                        Some(record_pair) => {
                            //For right record save end position and duplicate group ID
                            let record_end_pos = record.cigar_cached().unwrap().end_pos() - 1;
                            let group_id = match record_pair {
                                RecordStorage::PairedRecords { ref mut r_rec, .. } => {
                                    let group_id = duplicate_id.integer();
                                    r_rec.get_or_insert(IndexedRecord {
                                        rec: record,
                                        rec_id: i,
                                    });
                                    GroupID::Regular(group_id)
                                }
                                // This arm is reached if a mate is mapped to another chromosome.
                                // In that case a new duplicate and record ID is required
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
                            duplicate_groups
                                .entry(GroupID::Regular(duplicate_id.integer()))
                                .or_insert_with(Vec::new)
                                .push(RecordID::Regular(record_id.to_owned()));
                            if !record.is_paired()
                                || record.is_mate_unmapped()
                                || (record.tid() != record.mtid())
                            {
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
                                        l_rec: IndexedRecord {
                                            rec: record,
                                            rec_id: i,
                                        },
                                        r_rec: None,
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
                    if record.is_unmapped()
                        || record.is_mate_unmapped()
                        || (record.tid() != record.mtid())
                    {
                        self.bam_skipped_writer.write(&record)?;
                    } else {
                        match record_storage.get_mut(&RecordID::Regular(record_id.to_owned())) {
                            //Case: Left record
                            None => {
                                record_storage.insert(
                                    RecordID::Regular(record_id.to_owned()),
                                    RecordStorage::PairedRecords {
                                        l_rec: IndexedRecord {
                                            rec: record,
                                            rec_id: i,
                                        },
                                        r_rec: None,
                                    },
                                );
                            }
                            //Case: Left record already stored
                            Some(_record_pair) => {
                                let (rec_id, l_rec) = match record_storage
                                    .remove(&RecordID::Regular(record_id.to_owned()))
                                    .unwrap()
                                {
                                    RecordStorage::PairedRecords { l_rec, .. } => {
                                        (l_rec.rec_id, l_rec.into_rec())
                                    }
                                    RecordStorage::SingleRecord { .. } => unreachable!(),
                                };
                                if cigar_has_softclips(&l_rec) || cigar_has_softclips(&record) {
                                    self.bam_skipped_writer.write(&l_rec)?;
                                    self.bam_skipped_writer.write(&record)?;
                                } else {
                                    //TODO Alignment vectors need to include alignment of insertions and deletions
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
        let mut l_recs = Vec::new();
        let mut r_recs = Vec::new();
        let mut l_seqids = Vec::new();
        let mut r_seqids = Vec::new();
        for rec_id in duplicate_groups.remove(&group_id).unwrap() {
            match record_storage.remove(&rec_id).unwrap() {
                RecordStorage::PairedRecords { l_rec, r_rec } => {
                    let r_rec_unwraped = r_rec.unwrap();
                    let r_rec_id = r_rec_unwraped.rec_id;
                    let r_rec_enty = r_rec_unwraped.into_rec();
                    if cigar_has_softclips(&r_rec_enty) || cigar_has_softclips(&r_rec_enty) {
                        bam_skipped_writer.write(&l_rec.into_rec())?;
                        bam_skipped_writer.write(&r_rec_enty)?;
                    } else {
                        l_seqids.push(l_rec.rec_id);
                        l_recs.push(l_rec.into_rec());
                        r_seqids.push(r_rec_id);
                        r_recs.push(r_rec_enty);
                    }
                }
                RecordStorage::SingleRecord { rec } => {
                    let rec = rec.into_rec();
                    if cigar_has_softclips(&rec) {
                        bam_skipped_writer.write(&rec)?;
                    } else {
                        l_recs.push(rec);
                    }
                }
            };
        }

        if !r_recs.is_empty() {
            //TODO Alignment vectors need to include alignment of insertions and deletions
            let alignment_vectors = calc_read_alignments(&l_recs[0], &r_recs[0]);
            match alignment_vectors {
                Some((r1_alignment, r2_alignment)) => {
                    let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                    l_seqids.append(&mut r_seqids);
                    fq_se_writer.write_record(
                        &CalcOverlappingConsensus::new(
                            &l_recs,
                            &r_recs,
                            &r1_alignment,
                            &r2_alignment,
                            &l_seqids,
                            uuid,
                            verbose_read_names,
                        )
                        .calc_consensus()
                        .0,
                    )?;
                }
                None => {
                    let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                    fq1_writer.write_record(
                        &CalcNonOverlappingConsensus::new(&l_recs, &l_seqids, uuid)
                            .calc_consensus()
                            .0,
                    )?;
                    fq2_writer.write_record(
                        &CalcNonOverlappingConsensus::new(&r_recs, &r_seqids, uuid)
                            .calc_consensus()
                            .0,
                    )?;
                }
            };
        } else if !l_recs.is_empty() {
            //TODO Calculate alignment vectors
            let uuid = &Uuid::new_v4().to_hyphenated().to_string();
            fq_se_writer.write_record(
                &CalcNonOverlappingConsensus::new(&l_recs, &l_seqids, uuid)
                    .calc_consensus()
                    .0,
            )?;
        }
    }
    Ok(())
}

fn calc_read_alignments(
    r1_rec: &bam::Record,
    r2_rec: &bam::Record,
) -> Option<(Vec<bool>, Vec<bool>)> {
    let r1_start_softclips = r1_rec.cigar_cached().unwrap().leading_softclips();
    let r1_start = r1_rec.pos() - r1_start_softclips as i64;

    let r1_end_softclips = r1_rec.cigar_cached().unwrap().trailing_softclips();
    let r1_end = r1_rec.cigar_cached().unwrap().end_pos() + r1_end_softclips as i64;

    let r2_start_softclips = r2_rec.cigar_cached().unwrap().leading_softclips();
    let r2_start = r2_rec.pos() - r2_start_softclips as i64;

    let r2_end_softclips = r2_rec.cigar_cached().unwrap().trailing_softclips();
    let r2_end = r1_rec.cigar_cached().unwrap().end_pos() + r2_end_softclips as i64;

    if r1_start <= r2_start {
        //Check if reads overlap
        if r1_end >= r2_start {
            let offset = r2_start - r1_start;
            let (r1_vec, r2_vec) = calc_alignment_vectors(offset, r1_rec, r2_rec);
            Some((r1_vec, r2_vec))
        } else {
            //Reads do not overlap
            None
        }
    } else {
        //R2 starts before R1
        if r2_end >= r1_start {
            let offset = r1_start - r2_start;
            let (r2_vec, r1_vec) = calc_alignment_vectors(offset, r2_rec, r1_rec);
            Some((r1_vec, r2_vec))
        } else {
            None
        }
    }
}

fn calc_alignment_vectors(
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
}

fn cigar_has_softclips(rec: &bam::Record) -> bool {
    let cigar = rec.cigar_cached().unwrap();
    if cigar.leading_softclips() > 0 || cigar.trailing_softclips() > 0 {
        return true;
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
        l_rec: IndexedRecord,
        r_rec: Option<IndexedRecord>,
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
