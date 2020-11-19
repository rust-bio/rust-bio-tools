mod calc_consensus;
mod pipeline;

use bio::io::fastq;
use log::info;
use pipeline::CallConsensusRead;
use rust_htslib::bam;
use rust_htslib::bam::{Format, Header, Read};
use std::error::Error;

pub fn call_consensus_reads_from_paths(
    bam_in: &str,
    fq1: &str,
    fq2: &str,
    fq_se: &str,
    bam_skipped_out: &str,
    verbose_read_names: bool,
) -> Result<(), Box<dyn Error>> {
    info!("Reading input files:\n    {}", bam_in);
    info!("Writing forward consensus reads to:\n    {}", fq1);
    info!("Writing reverse consensus reads to:\n    {}", fq2);
    info!("Writing single end consensus reads to:\n    {}", fq_se);
    info!("Writing skipped reads to:\n    {}", bam_skipped_out);
    let bam_reader = bam::Reader::from_path(bam_in)?;
    let fq1_writer = fastq::Writer::to_file(fq1)?;
    let fq2_writer = fastq::Writer::to_file(fq2)?;
    let fq_se_writer = fastq::Writer::to_file(fq_se)?;
    let bam_skipped_writer = bam::Writer::from_path(
        bam_skipped_out,
        &Header::from_template(bam_reader.header()),
        Format::BAM,
    )?;
    CallConsensusRead::new(
        bam_reader,
        fq1_writer,
        fq2_writer,
        fq_se_writer,
        bam_skipped_writer,
        verbose_read_names,
    )
    .call_consensus_reads()
}
