mod pipeline;

use bio::io::fastq;
use flate2::write::GzEncoder;
use flate2::Compression;
use pipeline::CallConsensusRead;
use rust_htslib::bam;
use std::error::Error;
use std::fs;

pub fn call_consensus_reads_from_paths(
    bam: &str,
    fq1_out: &str,
    fq2_out: &str,
    fq3_out: &str,
    seq_dist: usize,
) -> Result<(), Box<dyn Error>> {
    eprintln!("Reading input files:\n    {}", bam);
    eprintln!(
        "Writing output to:\n    {}\n    {}\n    {}",
        fq1_out, fq2_out, fq3_out
    );
    match (fq1_out.ends_with(".gz"), fq2_out.ends_with(".gz"), fq3_out.ends_with(".gz")) {
        (false, false, false) => CallConsensusRead::new(
            &mut bam::Reader::from_path(bam)?,
            &mut fastq::Writer::to_file(fq1_out)?,
            &mut fastq::Writer::to_file(fq2_out)?,
            &mut fastq::Writer::to_file(fq3_out)?,
            seq_dist,
        ).call_consensus_reads(),
        (true, true, true) => CallConsensusRead::new(
            &mut bam::Reader::from_path(bam)?,
            &mut fastq::Writer::new(GzEncoder::new(fs::File::create(fq1_out)?, Compression::default())),
            &mut fastq::Writer::new(GzEncoder::new(fs::File::create(fq2_out)?, Compression::default())),
            &mut fastq::Writer::new(GzEncoder::new(fs::File::create(fq3_out)?, Compression::default())),
            seq_dist,
        ).call_consensus_reads(),
        _ => panic!("Invalid combination of files. Each pair of files (output) need to be both gzipped or both not zipped.")
    }
}
