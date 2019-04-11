use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::{HashMap, HashSet};
use std::error::Error;

pub struct CallConsensusRead<'a> {
    bam_reader: &'a mut bam::Reader,
    bam_writer: &'a mut bam::Writer,
    seq_dist: usize,
}

impl<'a> CallConsensusRead<'a> {
    pub fn new(
        bam_reader: &'a mut bam::Reader,
        bam_writer: &'a mut bam::Writer,
        seq_dist: usize,
    ) -> Self {
        CallConsensusRead {
            bam_reader,
            bam_writer,
            seq_dist,
        }
    }
    pub fn call_consensus_reads(&'a mut self) -> Result<(), Box<dyn Error>> {
        let mut group_end_idx: HashMap<i32, HashSet<i32>> = HashMap::new();
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
                        Some(read_pair) => {
                            //TODO If reverse read add end pos to group_end_pos
                            read_pair.reverse_record(record);
                        }
                        None => {
                            //TODO Add read_id to duplicate_groups
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
                        let insert_size = record.cigar().end_pos()? - f_rec.pos();
                        //TODO Calc overlapping consensus
                        //TODO Write to bam writer
                    }
                },
            }
        }
        Ok(())
    }
}

pub struct PairedReads {
    f_rec: bam::Record,
    r_rec: Option<bam::Record>,
}

impl PairedReads {
    pub fn reverse_record(&mut self, record: bam::Record) {
        self.r_rec = Some(record)
    }
}

pub struct GroupData {
    read_ids: Vec<Vec<u8>>,
    r_rec: Option<usize>,
}