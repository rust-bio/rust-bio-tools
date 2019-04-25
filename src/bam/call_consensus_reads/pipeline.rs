use super::calc_consensus::{CalcNonOverlappingConsensus, CalcOverlappingConsensus};
use indicatif;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::io::Write;
use std::process::{Command, Stdio};
use uuid::Uuid;

pub struct CallConsensusRead {
    bam_reader: bam::Reader,
    bam_writer: bam::Writer,
    seq_dist: usize,
}

type Position = i32;
type GroupID = i64;
type GroupIDs = HashSet<GroupID>;
type ReadIDs = Vec<ReadID>;
type ReadID = Vec<u8>;

impl CallConsensusRead {
    pub fn new(bam_reader: bam::Reader, bam_writer: bam::Writer, seq_dist: usize) -> Self {
        CallConsensusRead {
            bam_reader,
            bam_writer,
            seq_dist,
        }
    }
    pub fn call_consensus_reads(&mut self) -> Result<(), Box<dyn Error>> {
        let mut group_end_idx: BTreeMap<Position, GroupIDs> = BTreeMap::new();
        let mut duplicate_groups: HashMap<GroupID, ReadIDs> = HashMap::new();
        let mut record_storage: HashMap<ReadID, RecordStorage> = HashMap::new();

        for result in self.bam_reader.records() {
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
            let read_id = record.qname();
            //Check if record has duplicate ID
            match duplicate_id_option {
                //Case: duplicate ID exists
                Some(duplicate_id) => {
                    match record_storage.get_mut(read_id) {
                        //Case: Right record
                        Some(record_pair) => {
                            //For reverse read save end position and duplicate group ID
                            group_end_idx
                                .entry(record.cigar_cached().unwrap().end_pos()? - 1)
                                .or_insert_with(HashSet::new)
                                .insert(duplicate_id.integer());
                            match record_pair {
                                RecordStorage::PairedReads { ref mut r_rec, .. } => {
                                    r_rec.get_or_insert(record)
                                }
                                RecordStorage::SingleRead { .. } => unreachable!(),
                            };
                        }
                        //Case: Left record or record w/o mate
                        None => {
                            duplicate_groups
                                .entry(duplicate_id.integer())
                                .or_insert_with(Vec::new)
                                .push(read_id.to_vec());
                            if !record.is_paired() || record.is_mate_unmapped() {
                                //For reverse read save end position and duplicate group ID
                                group_end_idx
                                    .entry(record.cigar_cached().unwrap().end_pos()? - 1)
                                    .or_insert_with(HashSet::new)
                                    .insert(duplicate_id.integer());
                                record_storage.insert(
                                    read_id.to_vec(),
                                    RecordStorage::SingleRead { rec: record },
                                );
                            } else {
                                record_storage.insert(
                                    read_id.to_vec(),
                                    RecordStorage::PairedReads {
                                        f_rec: record,
                                        r_rec: None,
                                    },
                                );
                            }
                        }
                    }
                }
                //Duplicate ID not existing
                //Record is writen to bam file if read or its mate is unmapped
                //If Read is reverse and has mate consensus is calculated
                //Else record is added to hashMap
                None => {
                    if record.is_mate_unmapped() {
                        self.bam_writer.write(&record)?;
                    } else {
                        match record_storage.get_mut(read_id) {
                            //Case: Left record
                            None => {
                                record_storage.insert(
                                    read_id.to_vec(),
                                    RecordStorage::PairedReads {
                                        f_rec: record,
                                        r_rec: None,
                                    },
                                );
                            }
                            //Case: Left record already stored
                            Some(_record_pair) => {
                                let f_rec = match record_storage.remove(read_id).unwrap() {
                                    RecordStorage::PairedReads { f_rec, .. } => f_rec,
                                    RecordStorage::SingleRead { .. } => unreachable!(),
                                };
                                //TODO Consider softclipping in Overlap
                                let overlap =
                                    f_rec.cigar_cached().unwrap().end_pos()? - record.pos();
                                if overlap > 0 {
                                    let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                                    self.bam_writer.write(
                                        &CalcOverlappingConsensus::new(
                                            &[f_rec],
                                            &[record],
                                            overlap as usize,
                                            uuid,
                                        )
                                        .calc_consensus()
                                        .0,
                                    )?;
                                } else {
                                    self.bam_writer.write(&f_rec)?;
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
        )?;
        Ok(())
    }
}

pub fn calc_consensus_complete_groups(
    group_end_idx: &mut BTreeMap<Position, GroupIDs>,
    duplicate_groups: &mut HashMap<GroupID, ReadIDs>,
    end_pos: Option<&i32>,
    record_storage: &mut HashMap<Vec<u8>, RecordStorage>,
    bam_writer: &mut bam::Writer,
    seq_dist: usize,
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
        "Clustering duplicated reads by sequence using starcode."
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
        let read_ids = duplicate_groups.remove(&group_id).unwrap();
        for read_id in read_ids {
            let seq = match record_storage.get(&read_id).unwrap() {
                RecordStorage::PairedReads { f_rec, r_rec } => [
                    &f_rec.seq().as_bytes()[..],
                    &r_rec.clone().unwrap().seq().as_bytes()[..],
                ]
                .concat(),
                RecordStorage::SingleRead { rec } => rec.seq().as_bytes(),
            };
            seq_cluster.stdin.as_mut().unwrap().write(&seq)?;
            seq_cluster.stdin.as_mut().unwrap().write(b"\n")?;
            read_id_storage.push(read_id);
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
            let mut f_recs = Vec::new();
            let mut r_recs = Vec::new();
            for seqid in seqids {
                match record_storage.remove(&read_id_storage[seqid - 1]).unwrap() {
                    RecordStorage::PairedReads { f_rec, r_rec } => {
                        f_recs.push(f_rec);
                        r_recs.push(r_rec.unwrap());
                    }
                    RecordStorage::SingleRead { rec } => f_recs.push(rec)
                };
            }
            if !r_recs.is_empty() {
                let overlap = f_recs[0].cigar_cached().unwrap().end_pos()? - r_recs[0].pos(); //TODO
                if overlap > 0 {
                    let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                    bam_writer.write(
                        &CalcOverlappingConsensus::new(&f_recs, &r_recs, overlap as usize, uuid)
                            .calc_consensus()
                            .0,
                    )?;
                } else {
                    let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                    bam_writer.write(
                        &CalcNonOverlappingConsensus::new(&f_recs, uuid)
                            .calc_consensus()
                            .0,
                    )?;
                    let r_uuid = &Uuid::new_v4().to_hyphenated().to_string();
                    bam_writer.write(
                        &CalcNonOverlappingConsensus::new(&r_recs, r_uuid)
                            .calc_consensus()
                            .0,
                    )?;
                }
            } else {
                let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                bam_writer.write(
                    &CalcNonOverlappingConsensus::new(&f_recs, uuid)
                        .calc_consensus()
                        .0,
                )?;
            }
        }
    }
    Ok(())
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
    PairedReads {
        f_rec: bam::Record,
        r_rec: Option<bam::Record>,
    },
    SingleRead {
        rec: bam::Record,
    },
}
