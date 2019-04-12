use super::calc_consensus::{CalcNonOverlappingConsensus, CalcOverlappingConsensus};
use rust_htslib::bam;
use rust_htslib::bam::header::Header;
use rust_htslib::bam::Read;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use uuid::Uuid;

pub struct CallConsensusRead<'a> {
    bam_reader: &'a mut bam::Reader,
    bam_out: &'a str,
    seq_dist: usize,
}

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
        let mut group_end_idx: BTreeMap<i32, HashSet<i32>> = BTreeMap::new();
        let mut duplicate_groups: HashMap<i32, GroupData> = HashMap::new();
        let mut read_pairs: HashMap<Vec<u8>, PairedReads> = HashMap::new();

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
                            match duplicate_groups.get_mut(&duplicate_id.integer()) {
                                //TODO None case only occures if no forward Read exists
                                None => {}
                                Some(group_data) => {
                                    group_data.end_pos = Some(&record.cigar().end_pos()? - 1)
                                }
                            }
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
                            calc_consensus_complete_groups(&group_end_idx, &record.pos());
                            group_end_idx = group_end_idx.split_off(&record.pos()); //Remove processed indexes

                            match duplicate_groups.get_mut(&duplicate_id.integer()) {
                                None => {
                                    duplicate_groups.insert(
                                        duplicate_id.integer(),
                                        GroupData {
                                            read_ids: vec![read_id.to_vec()],
                                            end_pos: None,
                                        },
                                    );
                                }
                                Some(group_data) => group_data.read_ids.push(read_id.to_vec()),
                            }
                            read_pairs.insert(
                                read_id.to_vec(),
                                PairedReads {
                                    f_rec: record,
                                    r_rec: None,
                                },
                            );

                            //TODO If read start pos > group_end_pos calculate overlap for these groups
                        }
                    }
                }
                //If duplicate ID not existing add record to hashMap if not existing
                //else calc consensus and write to bam file
                None => match read_pairs.get_mut(read_id) {
                    None => {
                        calc_consensus_complete_groups(&group_end_idx, &record.pos());
                        group_end_idx = group_end_idx.split_off(&record.pos()); //Remove processed indexes
                        read_pairs.insert(
                            read_id.to_vec(),
                            PairedReads {
                                f_rec: record,
                                r_rec: None,
                            },
                        );
                    }
                    Some(read_pair) => {
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

pub fn calc_consensus_complete_groups(group_end_idx: &BTreeMap<i32, HashSet<i32>>, end_pos: &i32) {
    let end_idxs = group_end_idx.range(..end_pos).for_each(|(_, group_ids)| {
        group_ids.iter().map(|group_id| {});
    });
}

pub struct PairedReads {
    f_rec: bam::Record,
    r_rec: Option<bam::Record>,
}

pub struct GroupData {
    read_ids: Vec<Vec<u8>>,
    end_pos: Option<i32>,
}
