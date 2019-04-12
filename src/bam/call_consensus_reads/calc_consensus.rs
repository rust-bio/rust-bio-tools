use bio::stats::probs::{LogProb, PHREDProb};
use rust_htslib::bam;

pub struct CalcOverlappingConsensus<'a> {
    recs1: &'a [bam::Record],
    recs2: &'a [bam::Record],
    uuid: &'a str,
}

impl<'a> CalcOverlappingConsensus<'a> {
    pub fn new(recs1: &'a [bam::Record], recs2: &'a [bam::Record], uuid: &'a str) -> Self {
        CalcOverlappingConsensus { recs1, recs2, uuid }
    }
    pub fn calc_consensus(&self) -> (bam::Record, LogProb) {
        (bam::Record::new(), LogProb::ln_one()) //Placeholder
    }
}

pub struct CalcNonOverlappingConsensus<'a> {
    recs: &'a [bam::Record],
    uuid: &'a str,
}

impl<'a> CalcNonOverlappingConsensus<'a> {
    pub fn new(recs: &'a [bam::Record], uuid: &'a str) -> Self {
        CalcNonOverlappingConsensus { recs, uuid }
    }
    pub fn calc_consensus(&self) -> (bam::Record, LogProb) {
        (bam::Record::new(), LogProb::ln_one()) //Placeholder
    }
}
