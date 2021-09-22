use anyhow::Result;
use bio::io::fasta;
use rand::prelude::ThreadRng;
use rand::Rng;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::path::Path;
use uuid::Uuid;

pub fn simulate_reads<P: AsRef<Path>>(
    bam: P,
    input_ref: P,
    output_bam: P,
    output_ref: P,
    chr: String,
    start: u64,
    end: u64,
) -> Result<()> {
    let mut fasta_reader = fasta::IndexedReader::from_file(&input_ref)?;
    fasta_reader.fetch(&chr, start - 1, end - 1)?;
    let mut reference = Vec::new();
    fasta_reader.read(&mut reference)?;
    let mut rng = rand::thread_rng();
    let alphabet = [b'A', b'C', b'G', b'T'];

    //Build artificial reference
    let mut artificial_reference = Vec::new();
    add_random_bases(end - start, &mut artificial_reference, &mut rng, &alphabet)?;
    let mut fa_writer = fasta::Writer::to_file(output_ref)?;
    let ref_id = Uuid::new_v4().to_hyphenated().to_string();
    fa_writer.write(&ref_id, None, &artificial_reference)?;

    let mut bam_reader = bam::IndexedReader::from_path(bam)?;
    bam_reader.fetch((chr.as_bytes(), start - 1, end))?;

    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", &ref_id)
            .push_tag(b"LN", &(end - start)),
    );
    let mut bam_writer = bam::Writer::from_path(output_bam, &header, bam::Format::Bam)?;
    for result in bam_reader.records() {
        let mut record = result?;
        if (record.pos() >= (start - 1) as i64)
            && (record.cigar().end_pos() < (end - 1) as i64)
            && (record.mtid() == record.tid())
            && (record.mpos() + 1 >= (start as i64))
            && (record.mpos() + 1 < (end as i64))
        //TODO Check if mate is fully in interval
        {
            record.cache_cigar();
            //Check if mate record end within region
            let artificial_seq = if record.is_unmapped() {
                let mut seq = Vec::new();
                add_random_bases(record.seq_len() as u64, &mut seq, &mut rng, &alphabet)?;
                seq
            } else {
                build_sequence(
                    &reference,
                    &artificial_reference,
                    &record,
                    (start - 1) as usize,
                    &mut rng,
                    &alphabet,
                )?
            };
            let artificial_record = build_record(&record, &artificial_seq, (start - 1) as i64)?;
            bam_writer.write(&artificial_record)?;
        }
    }
    Ok(())
}

fn build_record(record: &bam::Record, artificial_seq: &[u8], offset: i64) -> Result<bam::Record> {
    let mut artificial_record = bam::record::Record::new();
    if let Ok(mate_cigar) = record.aux(b"MC") {
        artificial_record.push_aux(b"MC", mate_cigar)?;
    }
    artificial_record.set(
        record.qname(),
        Some(&record.cigar()),
        artificial_seq,
        record.qual(),
    );
    artificial_record.set_pos(record.pos() - offset);
    artificial_record.set_tid(0);
    artificial_record.set_mtid(0);
    artificial_record.set_mpos(record.mpos() - offset);
    artificial_record.set_flags(record.flags());
    artificial_record.set_insert_size(record.insert_size());
    artificial_record.set_mapq(record.mapq());
    Ok(artificial_record)
}

fn build_sequence(
    reference: &[u8],
    artificial_reference: &[u8],
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
    let mut ref_base = artificial_reference.get(ref_pos).unwrap();
    for cigar in record.cigar_cached().unwrap().iter() {
        match cigar.char() {
            'S' => {
                add_random_bases(cigar.len() as u64, &mut artificial_seq, rng, alphabet)?;
                record_pos += cigar.len() as usize;
            }
            'M' => {
                (0..cigar.len()).for_each(|_| {
                    ref_base = artificial_reference.get(ref_pos).unwrap();
                    if record_seq.get(record_pos).unwrap() == reference.get(ref_pos).unwrap() {
                        artificial_seq.push(*ref_base);
                    } else {
                        let mut reduced_alphabet = alphabet.to_vec();
                        reduced_alphabet.retain(|x| x != ref_base);
                        add_random_bases(1, &mut artificial_seq, rng, &reduced_alphabet).unwrap();
                    }
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
            'X' => {
                (0..cigar.len()).for_each(|_| {
                    ref_base = artificial_reference.get(ref_pos).unwrap();
                    let mut reduced_alphabet = alphabet.to_vec();
                    reduced_alphabet.retain(|x| x != ref_base);
                    add_random_bases(1, &mut artificial_seq, rng, &reduced_alphabet).unwrap();
                    ref_pos += 1;
                });
                record_pos += cigar.len() as usize;
            }
            '=' => {
                (0..cigar.len()).for_each(|_| {
                    let ref_base = artificial_reference.get(ref_pos).unwrap();
                    artificial_seq.push(*ref_base);
                    ref_pos += 1;
                });
                record_pos += cigar.len() as usize;
            }
            _ => {}
        }
    }

    Ok(artificial_seq)
}

fn add_random_bases(
    length: u64,
    seq: &mut Vec<u8>,
    rng: &mut ThreadRng,
    alphabet: &[u8],
) -> Result<()> {
    let alphabet_size = alphabet.len();
    (0..length).for_each(|_| {
        let idx = rng.gen_range(0, alphabet_size);
        seq.push(alphabet[idx])
    });
    Ok(())
}
