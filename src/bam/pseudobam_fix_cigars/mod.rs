use log::info;
use std::error::Error;
use std::str;

use rust_htslib::bam;
use rust_htslib::bam::{Read, Record};
use rust_htslib::bam::record::{Cigar, CigarString};

use bio::io::fasta;
use bio::alignment::{Alignment, pairwise};
use itertools::Itertools;
use std::collections::BTreeMap;

pub fn pseudobam_fix_cigars(
    fasta_in: &str,
    bam_in: &str,
    bam_out: &str,
) -> Result<(), Box<dyn Error>> {
    info!("Reading reference from:\n  {}", bam_in);
    let mut fasta_reader = fasta::IndexedReader::from_file(&fasta_in)?;

    info!("Reading pseudobam input file:\n  {}", bam_in);
    let bam_reader = bam::IndexedReader::from_path(&bam_in)?;

    info!("Writing output bam with accurate correct CIGAR strings to:\n  {}", bam_out);
    let mut bam_writer = bam::Writer::from_path(
        bam_out,
        &bam::Header::from_template(&bam_reader.header()),
        bam::Format::BAM,
    )?;

    let mut bam_buffer = bam::buffer::RecordBuffer::new(bam_reader, false);

    // values as used by bwa, see: http://bio-bwa.sourceforge.net/bwa.shtml
    let scoring = pairwise::Scoring {
        gap_open: -6,
        gap_extend: -1,
        match_fn: |a: u8, b: u8| if a == b {1i32} else {-4i32},
        match_scores: Some((1, -4)),
        xclip_prefix: -5,
        xclip_suffix: -5,
        yclip_prefix: 0,
        yclip_suffix: 0
    };
    let mut aligner = pairwise::Aligner::with_scoring(scoring);
    let check_upstream= 100;

    let mut reference_seq = Vec::new();
    for reference_name in bam_writer.header().target_names() {
        if let Ok(_) = fasta_reader.fetch_all(str::from_utf8(reference_name)? ) {
            let mut first_in_pairs: BTreeMap<&[u8], Record> = BTreeMap::new();
            fasta_reader.read(&mut reference_seq)?;
            if let Ok((_, _)) = bam_buffer.fetch(reference_name, 0, reference_seq.len() as u32) {
                for record in bam_buffer.iter_mut().collect_vec() {
                    let mut new_record = record.clone();
                    let start_align: usize = if ( record.pos() as usize <= check_upstream ) {
                        0
                    } else {
                        record.pos() as usize - check_upstream
                    };
                    let alignment = aligner.semiglobal( &record.seq().as_bytes().as_slice(), &reference_seq[start_align..]);
                    // new function cigar_string_from_alignment_operations()
                    let cigar_string = CigarString::from_alignment(&alignment, false);
                    new_record.set(record.qname(), Some(&cigar_string), record.seq().encoded, record.qual());
                    new_record.set_pos( (aligner.ystart + start_align) as i32);
                    if new_record.is_paired() {
                        if new_record.is_first_in_template() {
                            // need to wait for the second in the template to
                            // include its position and insert size in the metadata
                            first_in_pairs.insert(new_record.qname(), new_record);
                        } else if new_record.is_last_in_template() {
                            let mut first_in_pair = first_in_pairs.remove(new_record.qname())?;
                            let insert_size = new_record.pos() - first_in_pair.pos() + ( (alignment.yend + start_align) as i32 - new_record.pos() );
                            first_in_pair.set_mpos(new_record.pos());
                            new_record.set_mpos(first_in_pair.pos());
                            first_in_pair.set_insert_size(insert_size);
                            new_record.set_insert_size(insert_size);
                            bam_writer.write(&first_in_pair)?;
                            bam_writer.write(&new_record)?;
                        } else {
                            panic!("Read is paired but neither first nor last in template.\n\
                                    This should not happen in a pseudobam where only one mapping\n\
                                    per transcript isoform is expected for each read in a pair.\n");
                        }
                    } else {
                        bam_writer.write(&new_record)?;
                    }
                }
            }
        }
    }

    Ok(())
}
