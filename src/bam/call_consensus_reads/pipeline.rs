use bio::io::fastq;
use rust_htslib::bam;
use std::error::Error;
use std::io;

pub struct CallConsensusRead<'a, W: io::Write> {
    bam_reader: &'a mut bam::Reader,
    fq1_writer: &'a mut fastq::Writer<W>,
    fq2_writer: &'a mut fastq::Writer<W>,
    fq3_writer: &'a mut fastq::Writer<W>,
    seq_dist: usize,
}

impl<'a, W: io::Write> CallConsensusRead<'a, W> {
    pub fn new(
        bam_reader: &'a mut bam::Reader,
        fq1_writer: &'a mut fastq::Writer<W>,
        fq2_writer: &'a mut fastq::Writer<W>,
        fq3_writer: &'a mut fastq::Writer<W>,
        seq_dist: usize,
    ) -> Self {
        CallConsensusRead {
            bam_reader,
            fq1_writer,
            fq2_writer,
            fq3_writer,
            seq_dist,
        }
    }
    pub fn call_consensus_reads(&'a mut self) -> Result<(), Box<dyn Error>> {
        Ok(())
    }
}
