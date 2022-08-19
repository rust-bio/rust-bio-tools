use anyhow::Result;
use bio::io::fasta;
use rand::prelude::{SliceRandom, ThreadRng};
use rand::seq::IteratorRandom;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::ops::Range;
use std::path::Path;
use uuid::Uuid;

pub fn anonymize_reads<P: AsRef<Path> + std::fmt::Debug>(
    bam: P,
    input_ref: P,
    output_bam: P,
    output_ref: P,
    chr: String,
    interval: Range<u64>,
    keep_only_pairs: bool,
) -> Result<()> {
    let start = interval.start;
    let end = interval.end;
    let mut fasta_reader = fasta::IndexedReader::from_file(&input_ref)?;
    fasta_reader.fetch(&chr, start, end)?;
    let mut reference = Vec::new();
    reference.resize((end - start) as usize, 0);
    fasta_reader.read(&mut reference)?;
    let mut rng = rand::thread_rng();
    let alphabet = [b'A', b'C', b'G', b'T'];

    //Build artificial reference
    let mut artificial_reference = Vec::new();
    add_random_bases(end - start, &mut artificial_reference, &mut rng, &alphabet)?;
    let mut altered_bases = init_altered_bases(&reference, &artificial_reference)?;
    let mut fa_writer = fasta::Writer::to_file(output_ref)?;
    let ref_id = Uuid::new_v4().to_hyphenated().to_string();
    fa_writer.write(&ref_id, None, &artificial_reference)?;

    let mut bam_reader = bam::IndexedReader::from_path(bam)?;
    bam_reader.fetch((chr.as_bytes(), start, end + 1))?;

    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", &ref_id)
            .push_tag(b"LN", &(end - start)),
    );
    let mut bam_writer = bam::Writer::from_path(output_bam, &header, bam::Format::Bam)?;
    let mate_in_range = |record: &bam::Record| -> bool {
        (record.mtid() == record.tid())
            && (record.mpos() >= (start as i64))
            && (record.mpos() < (end as i64))
    };
    for result in bam_reader.records() {
        let mut record = result?;
        if (record.pos() >= start as i64)
            && (record.cigar().end_pos() < end as i64)
            && (!keep_only_pairs || mate_in_range(&record))
        {
            record.cache_cigar();
            //Check if mate record end within region
            let artificial_seq = if record.is_unmapped() || record.seq_len() == 0 {
                let mut seq = Vec::new();
                add_random_bases(record.seq_len() as u64, &mut seq, &mut rng, &alphabet)?;
                seq
            } else {
                build_sequence(
                    &mut altered_bases,
                    &record,
                    start as usize,
                    &mut rng,
                    &alphabet,
                )?
            };
            let artificial_record = build_record(&record, &artificial_seq, start as i64)?;
            bam_writer.write(&artificial_record)?;
        }
    }
    Ok(())
}

fn init_altered_bases(
    original_ref: &[u8],
    artificial_reference: &[u8],
) -> Result<HashMap<usize, HashMap<u8, u8>>> {
    let mut altered_bases = HashMap::new();
    for (i, (artifical_base, original_base)) in artificial_reference
        .iter()
        .zip(original_ref.iter())
        .enumerate()
    {
        altered_bases
            .entry(i)
            .or_insert_with(HashMap::new)
            .insert(*original_base, *artifical_base);
    }
    Ok(altered_bases)
}

fn build_record(record: &bam::Record, artificial_seq: &[u8], offset: i64) -> Result<bam::Record> {
    let mut artificial_record = bam::record::Record::new();
    artificial_record.set(
        record.qname(),
        Some(&record.cigar()),
        artificial_seq,
        record.qual(),
    );
    set_mandatory_fields(&mut artificial_record, record, offset)?;
    for aux_result in record.aux_iter() {
        let (tag, aux_field) = aux_result?;
        artificial_record.push_aux(tag, aux_field)?;
    }
    Ok(artificial_record)
}

fn build_sequence(
    altered_bases: &mut HashMap<usize, HashMap<u8, u8>>,
    record: &bam::Record,
    offset: usize,
    rng: &mut ThreadRng,
    alphabet: &[u8],
) -> Result<Vec<u8>> {
    let mut artificial_seq = Vec::new();
    let record_seq = record.seq().as_bytes();
    let mut record_pos = 0;
    let mut ref_pos = record.pos() as usize - offset;
    //Create random seq for leading softclips
    for cigar in record.cigar_cached().unwrap().iter() {
        match cigar.char() {
            'S' => {
                add_random_bases(cigar.len() as u64, &mut artificial_seq, rng, alphabet)?;
                record_pos += cigar.len() as usize;
            }
            'M' | 'X' | '=' => {
                (0..cigar.len()).for_each(|_| {
                    let base_mappings = altered_bases.get(&ref_pos).unwrap().clone();
                    let altered_base = *altered_bases
                        .get_mut(&ref_pos)
                        .unwrap()
                        .entry(*record_seq.get(record_pos).unwrap())
                        .or_insert_with(|| {
                            *alphabet
                                .iter()
                                .filter(|&x| !base_mappings.values().any(|y| x == y))
                                .choose(rng)
                                .unwrap()
                        });
                    artificial_seq.push(altered_base);
                    ref_pos += 1;
                    record_pos += 1;
                });
                // Add reference bases except for mismatches
            }
            'I' => {
                add_random_bases(cigar.len() as u64, &mut artificial_seq, rng, alphabet)?;
                record_pos += cigar.len() as usize;
            }
            'D' | 'N' => {
                ref_pos += cigar.len() as usize;
            }
            _ => {}
        }
    }

    Ok(artificial_seq)
}

fn set_mandatory_fields(
    target_rec: &mut bam::Record,
    source_rec: &bam::Record,
    offset: i64,
) -> Result<()> {
    target_rec.set_pos(source_rec.pos() - offset);
    target_rec.set_tid(0);
    let (mtid, mpos) = if source_rec.mtid() == -1 {
        (-1, -1)
    } else if source_rec.mtid() == source_rec.tid() {
        (0, source_rec.mpos() - offset)
    } else {
        (1, source_rec.mpos())
    };
    target_rec.set_mtid(mtid);
    target_rec.set_mpos(mpos);
    target_rec.set_flags(source_rec.flags());
    target_rec.set_insert_size(source_rec.insert_size());
    target_rec.set_mapq(source_rec.mapq());
    Ok(())
}

fn add_random_bases(
    length: u64,
    seq: &mut Vec<u8>,
    rng: &mut ThreadRng,
    alphabet: &[u8],
) -> Result<()> {
    (0..length).for_each(|_| seq.push(*alphabet.choose(rng).unwrap()));
    Ok(())
}
