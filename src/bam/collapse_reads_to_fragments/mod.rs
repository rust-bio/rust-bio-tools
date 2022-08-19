mod calc_consensus;
mod pipeline;

use anyhow::Result;
use bio::io::fastq;
use log::info;
use pipeline::CallConsensusRead;
use rust_htslib::bam;
use rust_htslib::bam::{Format, Header, Read};
use std::path::Path;

pub fn call_consensus_reads_from_paths<P: AsRef<Path>>(
    bam_in: P,
    fq1: P,
    fq2: P,
    fq_se: P,
    bam_skipped_out: P,
    verbose_read_names: bool,
) -> Result<()> {
    info!("Reading input files:\n    {}", bam_in.as_ref().display());
    info!(
        "Writing forward consensus reads to:\n    {}",
        fq1.as_ref().display()
    );
    info!(
        "Writing reverse consensus reads to:\n    {}",
        fq2.as_ref().display()
    );
    info!(
        "Writing single end consensus reads to:\n    {}",
        fq_se.as_ref().display()
    );
    info!(
        "Writing skipped reads to:\n    {}",
        bam_skipped_out.as_ref().display()
    );
    let bam_reader = bam::Reader::from_path(bam_in)?;
    let fq1_writer = fastq::Writer::to_file(fq1)?;
    let fq2_writer = fastq::Writer::to_file(fq2)?;
    let fq_se_writer = fastq::Writer::to_file(fq_se)?;
    let bam_skipped_writer = bam::Writer::from_path(
        bam_skipped_out,
        &Header::from_template(bam_reader.header()),
        Format::Bam,
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

pub fn unmark_record(record: &mut bam::record::Record) -> Result<()> {
    record.unset_duplicate();
    let _ = record.remove_aux(b"PG");
    let _ = record.remove_aux(b"DI");
    let _ = record.remove_aux(b"DS");
    Ok(())
}
