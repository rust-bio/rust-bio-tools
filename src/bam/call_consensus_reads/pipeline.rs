use rust_htslib::bam;
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
        // Loop over group of reads with same duplicate ID from bam_reader (Filter duplicate IDs?)
        // spawn starcode cluster
        // concat forward and reverse reads
        // add read to cluster
        // perform clustering
        // calc consensus for reads in each subcluster
        // write consensus to output
        // go on with next group of duplicates
        Ok(())
    }
}
