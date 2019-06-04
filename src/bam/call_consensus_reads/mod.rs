mod calc_consensus;
mod pipeline;

use log::info;
use pipeline::CallConsensusRead;
use rust_htslib::bam;
use rust_htslib::bam::{Header, Read};
use std::error::Error;

pub fn call_consensus_reads_from_paths(
    bam_in: &str,
    bam_out: &str,
    seq_dist: usize,
    verbose_read_names: bool,
) -> Result<(), Box<dyn Error>> {
    info!("Reading input files:\n    {}", bam_in);
    info!("Writing output to:\n    {}", bam_out);
    let bam_reader = bam::Reader::from_path(bam_in)?;
    let bam_writer = bam::Writer::from_path(bam_out, &Header::from_template(bam_reader.header()))?;
    CallConsensusRead::new(bam_reader, bam_writer, seq_dist, verbose_read_names)
        .call_consensus_reads()
}
