use super::calc_consensus::{CalcNonOverlappingConsensus, CalcOverlappingConsensus};
use rocksdb::DB;
use rust_htslib::bam;
use rust_htslib::bam::header::Header;
use rust_htslib::bam::Read;
use serde_json;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::io::Write;
use std::mem;
use std::process::{Command, Stdio};
use std::str;
use tempfile::tempdir;
use uuid::Uuid;

pub struct CallConsensusRead<'a> {
    bam_reader: &'a mut bam::Reader,
    bam_out: &'a str,
    seq_dist: usize,
}

type Position = i32;
type GroupID = i32;
type GroupIDs = HashSet<GroupID>;
type ReadIDs = Vec<ReadID>;
type ReadID = Vec<u8>;

impl<'a> CallConsensusRead<'a> {
    pub fn new(bam_reader: &'a mut bam::Reader, bam_out: &'a str, seq_dist: usize) -> Self {
        CallConsensusRead {
            bam_reader,
            bam_out,
            seq_dist,
        }
    }
    pub fn call_consensus_reads(&'a mut self) -> Result<(), Box<dyn Error>> {
        let mut bam_writer = bam::Writer::from_path(
            self.bam_out,
            &Header::from_template(self.bam_reader.header()),
        )?;
        let mut group_end_idx: BTreeMap<Position, GroupIDs> = BTreeMap::new();
        let mut duplicate_groups: HashMap<GroupID, ReadIDs> = HashMap::new();
        let mut read_pairs: HashMap<ReadID, PairedReads> = HashMap::new();

        for result in self.bam_reader.records() {
            let record = result?;
            let duplicate_id_option = record.aux(b"DI");
            let read_id = record.qname();
            //Check if record has duplicate ID
            match duplicate_id_option {
                // If duplicate ID exists add record to HashMap
                Some(duplicate_id) => {
                    match read_pairs.get_mut(read_id) {
                        //Reverse Read
                        Some(read_pair) => {
                            //For reverse read save end position and duplicate group ID
                            match group_end_idx.get_mut(&(&record.cigar().end_pos()? - 1)) {
                                None => {
                                    let mut group_set = HashSet::new();
                                    group_set.insert(duplicate_id.integer());
                                    group_end_idx.insert(&record.cigar().end_pos()? - 1, group_set);
                                }
                                Some(group_set) => {
                                    group_set.insert(duplicate_id.integer());
                                }
                            }
                            read_pair.r_rec = Some(record);
                        }
                        //Forward Read
                        //Structure should be done
                        None => {
                            //Process completed duplicate groups
                            calc_consensus_complete_groups(
                                &mut group_end_idx,
                                &mut duplicate_groups,
                                &record.pos(),
                                &mut read_pairs,
                                &mut bam_writer,
                                self.seq_dist,
                            )?;
                            group_end_idx = group_end_idx.split_off(&record.pos()); //Remove processed indexes

                            match duplicate_groups.get_mut(&duplicate_id.integer()) {
                                None => {
                                    duplicate_groups
                                        .insert(duplicate_id.integer(), vec![read_id.to_vec()]);
                                }
                                Some(read_ids) => read_ids.push(read_id.to_vec()),
                            }
                            read_pairs.insert(
                                read_id.to_vec(),
                                PairedReads {
                                    f_rec: record,
                                    r_rec: None,
                                },
                            );
                        }
                    }
                }
                //If duplicate ID not existing add record to hashMap if not existing
                //else calc consensus and write to bam file
                None => match read_pairs.get_mut(read_id) {
                    None => {
                        calc_consensus_complete_groups(
                            &mut group_end_idx,
                            &mut duplicate_groups,
                            &record.pos(),
                            &mut read_pairs,
                            &mut bam_writer,
                            self.seq_dist,
                        )?;
                        group_end_idx = group_end_idx.split_off(&record.pos()); //Remove processed indexes
                        read_pairs.insert(
                            read_id.to_vec(),
                            PairedReads {
                                f_rec: record,
                                r_rec: None,
                            },
                        );
                    }
                    Some(_read_pair) => {
                        let f_rec = read_pairs.remove(read_id).unwrap().f_rec;
                        if (record.seq().len() + f_rec.seq().len()) > f_rec.insert_size() as usize {
                            bam_writer.write(&f_rec)?;
                            bam_writer.write(&record)?;
                        } else {
                            let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                            bam_writer.write(
                                &CalcOverlappingConsensus::new(&[f_rec], &[record], uuid)
                                    .calc_consensus()
                                    .0,
                            )?;
                        }
                    }
                },
            }
        }
        Ok(())
    }
}

pub fn calc_consensus_complete_groups(
    group_end_idx: &mut BTreeMap<Position, GroupIDs>,
    duplicate_groups: &mut HashMap<GroupID, ReadIDs>,
    end_pos: &i32,
    read_pairs: &mut HashMap<Vec<u8>, PairedReads>,
    bam_writer: &mut bam::Writer,
    seq_dist: usize,
) -> Result<(), Box<dyn Error>> {
    let group_ids: HashSet<i32> = group_end_idx
        .range(..end_pos)
        .flat_map(|(_, group_ids)| group_ids.clone())
        .collect();
    for group_id in group_ids {
        let mut read_storage = BAMStorage::new()?;
        let mut i = 0;

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
            let paired_record = read_pairs.remove(&read_id).unwrap();
            let f_rec = paired_record.f_rec;
            let r_rec = paired_record.r_rec.unwrap();
            read_storage.put(i, &f_rec, &r_rec);
            let seq = [&f_rec.seq().as_bytes()[..], &r_rec.seq().as_bytes()[..]].concat();
            seq_cluster.stdin.as_mut().unwrap().write(&seq)?;
            seq_cluster.stdin.as_mut().unwrap().write(b"\n")?;
            i += 1;
        }
        seq_cluster.stdin.as_mut().unwrap().flush()?;
        drop(seq_cluster.stdin.take());
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
                let (f_rec, r_rec) = read_storage.get(seqid - 1)?;
                f_recs.push(f_rec);
                r_recs.push(r_rec);
            }

            //Todo Starcode Clustering
            if f_recs[0].seq().len() + r_recs[0].seq().len() > f_recs[0].insert_size() as usize {
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
            } else {
                let uuid = &Uuid::new_v4().to_hyphenated().to_string();
                bam_writer.write(
                    &CalcOverlappingConsensus::new(&f_recs, &r_recs, uuid)
                        .calc_consensus()
                        .0,
                )?;
            }
        }
    }
    Ok(())
}

/// Used to store a mapping of read index to read sequence
#[derive(Debug)]
struct BAMStorage {
    db: DB,
    storage_dir: std::path::PathBuf,
}

impl BAMStorage {
    /// Create a new FASTQStorage using a Rocksdb database
    /// that maps read indices to read seqeunces.
    pub fn new() -> Result<Self, Box<Error>> {
        // Save storage_dir to prevent it from leaving scope and
        // in turn deleting the tempdir
        let storage_dir = tempdir()?.path().join("db");
        Ok(BAMStorage {
            db: DB::open_default(storage_dir.clone())?,
            storage_dir,
        })
    }

    fn as_key<'a>(i: u64) -> [u8; 8] {
        unsafe { mem::transmute::<u64, [u8; 8]>(i) }
    }

    /// Enter a (read index, read sequence) pair into the database.
    pub fn put(
        &mut self,
        i: usize,
        f_rec: &bam::Record,
        r_rec: &bam::Record,
    ) -> Result<(), Box<dyn Error>> {
        Ok(self.db.put(
            &Self::as_key(i as u64),
            serde_json::to_string(&(f_rec, r_rec))?.as_bytes(),
        )?)
    }

    /// Retrieve the read sequence of the read with index `i`.
    pub fn get(&self, i: usize) -> Result<(bam::Record, bam::Record), Box<dyn Error>> {
        Ok(serde_json::from_str(
            str::from_utf8(&self.db.get(&Self::as_key(i as u64))?.unwrap()).unwrap(),
        )?)
    }
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

pub struct PairedReads {
    f_rec: bam::Record,
    r_rec: Option<bam::Record>,
}
