mod calc_consensus;
mod pipeline;

use pipeline::CallConsensusRead;
use rust_htslib::bam;
use std::error::Error;

pub fn call_consensus_reads_from_paths(
    bam_in: &str,
    bam_out: &str,
    seq_dist: usize,
) -> Result<(), Box<dyn Error>> {
    eprintln!("Reading input files:\n    {}", bam_in);
    eprintln!("Writing output to:\n    {}", bam_out);
    CallConsensusRead::new(&mut bam::Reader::from_path(bam_in)?, &bam_out, seq_dist)
        .call_consensus_reads()
}
