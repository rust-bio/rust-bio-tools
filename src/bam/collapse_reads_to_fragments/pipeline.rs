use super::calc_consensus::{CalcNonOverlappingConsensus, CalcOverlappingConsensus};
use anyhow::Result;
use bio::io::fastq;
use derive_new::new;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Read;
use std::collections::{BTreeMap, HashMap, HashSet};
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
    pub fn call_consensus_reads(&mut self) -> Result<()> {
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
                        //TODO Handle intersecting reads mapped on different chromosomes
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
                                let overlap_opt = calc_overlap(&l_rec, &record);
                                //TODO overlap_opt ist part of skipping softclips
                                // Handle soft clips later
                                if let Some(overlap) = overlap_opt {
                                    if overlap > 0
                                        && is_valid_overlap(
                                            overlap as u32,
                                            l_rec.cigar_cached().unwrap().into_iter().rev(),
                                        )
                                        && is_valid_overlap(
                                            overlap as u32,
                                            record.cigar_cached().unwrap().into_iter(),
                                        )
                                    {
                                        let uuid = &Uuid::new_v4().to_hyphenated().to_string();

                                        self.fq_se_writer.write_record(
                                            &CalcOverlappingConsensus::new(
                                                &[l_rec],
                                                &[record],
                                                overlap as usize,
                                                &[rec_id, i],
                                                uuid,
                                                self.verbose_read_names,
                                            )
                                            .calc_consensus()
                                            .0,
                                        )?;
                                    } else {
                                        self.bam_skipped_writer.write(&l_rec)?;
                                        self.bam_skipped_writer.write(&record)?;
                                    }
                                } else {
                                    self.bam_skipped_writer.write(&l_rec)?;
                                    self.bam_skipped_writer.write(&record)?;
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
    verbose_read_names: bool,
) -> Result<()> {
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
                    l_seqids.push(l_rec.rec_id);
                    l_recs.push(l_rec.into_rec());
                    r_seqids.push(r_rec.as_ref().unwrap().rec_id);
                    r_recs.push(r_rec.unwrap().into_rec());
                }
                RecordStorage::SingleRecord { rec } => l_recs.push(rec.into_rec()),
            };
        }

        if !r_recs.is_empty() {
            let overlap_opt = calc_overlap(&l_recs[0], &r_recs[0]);
            if let Some(overlap) = overlap_opt {
                if overlap > 0
                    && is_valid_overlap(
                        overlap as u32,
                        l_recs[0].cigar_cached().unwrap().into_iter().rev(),
                    )
                    && is_valid_overlap(
                        overlap as u32,
                        r_recs[0].cigar_cached().unwrap().into_iter(),
                    )
                {
                    let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                    l_seqids.append(&mut r_seqids);
                    fq_se_writer.write_record(
                        &CalcOverlappingConsensus::new(
                            &l_recs,
                            &r_recs,
                            overlap as usize,
                            &l_seqids,
                            uuid,
                            verbose_read_names,
                        )
                        .calc_consensus()
                        .0,
                    )?;
                } else {
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
            } else {
                unreachable!()
            }
        } else {
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

fn is_valid_overlap<'a, I>(overlap: u32, cigar: I) -> bool
where
    I: Iterator<Item = &'a Cigar>,
{
    let mut i = 0;
    for c in cigar {
        match i < overlap {
            true => match c {
                Cigar::Ins(_) | Cigar::Del(_) => return false,
                Cigar::Match(l)
                | Cigar::RefSkip(l)
                | Cigar::SoftClip(l)
                | Cigar::Pad(l)
                | Cigar::Equal(l)
                | Cigar::Diff(l) => i += l,
                Cigar::HardClip(_) => {}
            },
            false => return true,
        }
    }
    true
}

fn calc_overlap(l_rec: &bam::Record, r_rec: &bam::Record) -> Option<i64> {
    let l_start_softclips = count_softclips(l_rec.cigar_cached().unwrap().into_iter());
    let l_start_pos = l_rec.pos() - l_start_softclips as i64;

    let l_end_softclips = count_softclips(l_rec.cigar_cached().unwrap().into_iter().rev());
    let l_end_pos = l_rec.cigar_cached().unwrap().end_pos() + l_end_softclips as i64;

    let r_start_softclips = count_softclips(r_rec.cigar_cached().unwrap().into_iter());
    let r_start_pos = r_rec.pos() - r_start_softclips as i64;

    let r_end_softclips = count_softclips(r_rec.cigar_cached().unwrap().into_iter().rev());
    let r_end_pos = l_rec.cigar_cached().unwrap().end_pos() + r_end_softclips as i64;

    //TODO Skipping soft clips here. Handle this correctly
    if (l_start_softclips > 0)
        | (l_end_softclips > 0)
        | (r_start_softclips > 0)
        | (r_end_softclips > 0)
    {
        return None;
    }
    //TODO if-closure is just a hotfix to ensure reads only overlap by end of r1 and start of r2
    // or are at exact same position
    // Fix this later by handling any other alignments
    if l_end_pos <= r_end_pos && l_start_pos <= r_start_pos {
        let left_overlap_pos = match l_start_pos >= r_start_pos {
            true => l_start_pos,
            false => r_start_pos,
        };
        let right_overlap_pos = match l_end_pos <= r_end_pos {
            true => l_end_pos,
            false => r_end_pos,
        };
        Some(right_overlap_pos - left_overlap_pos)
    } else {
        None
    }
}

//Gets an Iterator over Cigar-items and returns number of soft-clips at the beginning
fn count_softclips<'a, I>(cigar: I) -> i32
where
    I: Iterator<Item = &'a Cigar>,
{
    for c in cigar {
        match c {
            Cigar::HardClip(_) => {}
            Cigar::SoftClip(l) => return *l as i32,
            _ => return 0,
        }
    }
    unreachable!();
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
