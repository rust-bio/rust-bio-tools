use super::calc_consensus::{CalcNonOverlappingConsensus, CalcOverlappingConsensus};
use derive_new::new;
use indicatif;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Read;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::io::Write;
use std::ops::Deref;
use std::process::{Command, Stdio};
use uuid::Uuid;

#[derive(new)]
pub struct CallConsensusRead {
    bam_reader: bam::Reader,
    bam_writer: bam::Writer,
    seq_dist: usize,
    verbose_read_names: bool,
}

type Position = i32;
type GroupID = i64;
type GroupIDs = HashSet<GroupID>;
type RecordIDS = Vec<RecordID>;
type RecordID = Vec<u8>;

impl CallConsensusRead {
    pub fn call_consensus_reads(&mut self) -> Result<(), Box<dyn Error>> {
        let mut group_end_idx: BTreeMap<Position, GroupIDs> = BTreeMap::new();
        let mut duplicate_groups: HashMap<GroupID, RecordIDS> = HashMap::new();
        let mut record_storage: HashMap<RecordID, RecordStorage> = HashMap::new();

        for (i, result) in self.bam_reader.records().into_iter().enumerate() {
            let mut record = result?;
            if !record.is_unmapped() {
                //Process completed duplicate groups
                calc_consensus_complete_groups(
                    &mut group_end_idx,
                    &mut duplicate_groups,
                    Some(&record.pos()),
                    &mut record_storage,
                    &mut self.bam_writer,
                    self.seq_dist,
                    self.verbose_read_names,
                )?;
                group_end_idx = group_end_idx.split_off(&record.pos()); //Remove processed indexes
            } else {
                self.bam_writer.write(&record)?;
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
                    match record_storage.get_mut(record_id) {
                        //Case: Right record
                        Some(record_pair) => {
                            //For right record save end position and duplicate group ID
                            group_end_idx
                                .entry(record.cigar_cached().unwrap().end_pos()? - 1)
                                .or_insert_with(HashSet::new)
                                .insert(duplicate_id.integer());
                            match record_pair {
                                RecordStorage::PairedRecords { ref mut r_rec, .. } => r_rec
                                    .get_or_insert(IndexedRecord {
                                        rec: record,
                                        rec_id: i,
                                    }),
                                RecordStorage::SingleRecord { .. } => unreachable!(),
                            };
                        }
                        //Case: Left record or record w/o mate
                        None => {
                            duplicate_groups
                                .entry(duplicate_id.integer())
                                .or_insert_with(Vec::new)
                                .push(record_id.to_vec());
                            if !record.is_paired() || record.is_mate_unmapped() {
                                //If right or single record save end position and duplicate group ID
                                group_end_idx
                                    .entry(record.cigar_cached().unwrap().end_pos()? - 1)
                                    .or_insert_with(HashSet::new)
                                    .insert(duplicate_id.integer());
                                record_storage.insert(
                                    record_id.to_vec(),
                                    RecordStorage::SingleRecord {
                                        rec: IndexedRecord {
                                            rec: record,
                                            rec_id: i,
                                        },
                                    },
                                );
                            } else {
                                record_storage.insert(
                                    record_id.to_vec(),
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
                    if record.is_mate_unmapped() {
                        self.bam_writer.write(&record)?;
                    } else {
                        match record_storage.get_mut(record_id) {
                            //Case: Left record
                            None => {
                                record_storage.insert(
                                    record_id.to_vec(),
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
                                let (rec_id, l_rec) =
                                    match record_storage.remove(record_id).unwrap() {
                                        RecordStorage::PairedRecords { l_rec, .. } => {
                                            (l_rec.rec_id, l_rec.into_rec())
                                        }
                                        RecordStorage::SingleRecord { .. } => unreachable!(),
                                    };
                                let overlap = calc_overlap(&l_rec, &record)?;
                                if overlap > 0 {
                                    let uuid = &Uuid::new_v4().to_hyphenated().to_string();

                                    self.bam_writer.write(
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
                                    self.bam_writer.write(&l_rec)?;
                                    self.bam_writer.write(&record)?;
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
            &mut self.bam_writer,
            self.seq_dist,
            self.verbose_read_names,
        )?;
        Ok(())
    }
}

pub fn calc_consensus_complete_groups(
    group_end_idx: &mut BTreeMap<Position, GroupIDs>,
    duplicate_groups: &mut HashMap<GroupID, RecordIDS>,
    end_pos: Option<&i32>,
    record_storage: &mut HashMap<Vec<u8>, RecordStorage>,
    bam_writer: &mut bam::Writer,
    seq_dist: usize,
    verbose_read_names: bool,
) -> Result<(), Box<dyn Error>> {
    let spinner_style = indicatif::ProgressStyle::default_spinner()
        .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ")
        .template("{prefix:.bold.dim} {spinner} {wide_msg}");

    let group_ids: HashSet<i64> = group_end_idx
        .range(..end_pos.unwrap_or(&(group_end_idx.len() as i32)))
        .flat_map(|(_, group_ids)| group_ids.clone())
        .collect();

    // prepare spinner for user feedback
    let pb = indicatif::ProgressBar::new_spinner();
    pb.set_style(spinner_style.clone());
    pb.set_prefix(&format!(
        "Clustering duplicated records by sequence using starcode."
    ));
    for (i, group_id) in group_ids.into_iter().enumerate() {
        pb.inc(1);
        let mut read_id_storage = Vec::new();
        let mut seq_cluster = Command::new("starcode")
            .arg("--dist")
            .arg(format!("{}", seq_dist))
            .arg("--seq-id")
            .arg("-s")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()?;
        let rec_ids = duplicate_groups.remove(&group_id).unwrap();
        for rec_id in rec_ids {
            let seq = match record_storage.get(&rec_id).unwrap() {
                RecordStorage::PairedRecords { l_rec, r_rec } => [
                    &l_rec.seq().as_bytes()[..],
                    &r_rec.as_ref().unwrap().seq().as_bytes()[..],
                ]
                .concat(),
                RecordStorage::SingleRecord { rec } => rec.seq().as_bytes(),
            };
            seq_cluster.stdin.as_mut().unwrap().write(&seq)?;
            seq_cluster.stdin.as_mut().unwrap().write(b"\n")?;
            read_id_storage.push(rec_id);
        }
        seq_cluster.stdin.as_mut().unwrap().flush()?;
        drop(seq_cluster.stdin.take());
        pb.finish_with_message(&format!("Done. Analyzed {} group(s).", i));
        for record in csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(seq_cluster.stdout.as_mut().unwrap())
            .records()
        {
            let seqids = parse_cluster(record?)?;
            let mut l_recs = Vec::new();
            let mut r_recs = Vec::new();
            let mut l_seqids = Vec::new();
            let mut r_seqids = Vec::new();
            for seqid in seqids {
                match record_storage.remove(&read_id_storage[seqid - 1]).unwrap() {
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
                let overlap = calc_overlap(&l_recs[0], &r_recs[0])?;
                if overlap > 0 {
                    let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                    l_seqids.append(&mut r_seqids);
                    bam_writer.write(
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
                    bam_writer.write(
                        &CalcNonOverlappingConsensus::new(
                            &l_recs,
                            &l_seqids,
                            uuid,
                            verbose_read_names,
                        )
                        .calc_consensus()
                        .0,
                    )?;
                    let r_uuid = &Uuid::new_v4().to_hyphenated().to_string();
                    bam_writer.write(
                        &CalcNonOverlappingConsensus::new(
                            &r_recs,
                            &r_seqids,
                            r_uuid,
                            verbose_read_names,
                        )
                        .calc_consensus()
                        .0,
                    )?;
                }
            } else {
                let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                bam_writer.write(
                    &CalcNonOverlappingConsensus::new(&l_recs, &l_seqids, uuid, verbose_read_names)
                        .calc_consensus()
                        .0,
                )?;
            }
        }
    }
    Ok(())
}

fn calc_overlap(l_rec: &bam::Record, r_rec: &bam::Record) -> Result<i32, Box<dyn Error>> {
    let l_end_pos = l_rec.cigar_cached().unwrap().end_pos()?;
    let r_start_pos = r_rec.pos();
    let l_softclips = count_softclips(l_rec.cigar_cached().unwrap().into_iter().rev())?;
    let r_softclips = count_softclips(r_rec.cigar_cached().unwrap().into_iter())?;
    Ok((l_end_pos + l_softclips) - (r_start_pos - r_softclips))
}

//Gets an Iterator over Cigar-items and returns number of soft-clips at the beginning
fn count_softclips<'a, I>(cigar: I) -> Result<i32, Box<dyn Error>>
where
    I: Iterator<Item = &'a Cigar>,
{
    for c in cigar {
        match c {
            Cigar::HardClip(_) => {}
            Cigar::SoftClip(l) => return Ok(*l as i32),
            _ => return Ok(0),
        }
    }
    unreachable!();
}

/// Interpret a cluster returned by starcode
fn parse_cluster(record: csv::StringRecord) -> Result<Vec<usize>, Box<dyn Error>> {
    let seqids = &record[2];
    Ok(csv::ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(false)
        .from_reader(seqids.as_bytes())
        .deserialize()
        .next()
        .unwrap()?)
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
