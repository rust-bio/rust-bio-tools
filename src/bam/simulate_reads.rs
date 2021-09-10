use anyhow::Result;
use bio::io::fasta;
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
    fasta_reader.fetch(&chr, start, end)?;
    let mut reference = Vec::new();
    fasta_reader.read(&mut reference)?;

    //Build artificial reference
    let mut artificial_reference = Vec::new();
    add_random_bases(end - start, &mut artificial_reference)?;
    let mut fa_writer = fasta::Writer::to_file(output_ref)?;
    let ref_id = Uuid::new_v4().to_hyphenated().to_string();
    fa_writer.write(&ref_id, None, &artificial_reference)?;

    let mut bam_reader = bam::IndexedReader::from_path(bam)?;
    bam_reader.fetch((chr.as_bytes(), start, end))?;

    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", &ref_id)
            .push_tag(b"LN", &(end - start)),
    );
    let mut bam_writer = bam::Writer::from_path(output_bam, &header, bam::Format::BAM)?;
    for result in bam_reader.records() {
        let mut record = result?;
        if (record.cigar().end_pos() < (end as i64))
            && (record.mtid() == record.tid())
            && (record.mpos() >= (start as i64))
            && (record.mpos() < (end as i64))
        {
            record.cache_cigar();
            //Check if mate record end within region
            let artifical_seq =
                build_sequence(&reference, &artificial_reference, &record, start as usize)?;
            let mut artifical_record = bam::record::Record::new();
            artifical_record.set(
                record.qname(),
                Some(&record.cigar()),
                &artifical_seq,
                record.qual(),
            );
            artifical_record.set_pos(record.pos() - start as i64);
            artifical_record.set_tid(0);
            artifical_record.set_mtid(0);
            artifical_record.set_mpos(record.mpos() - start as i64);
            bam_writer.write(&artifical_record)?;
        }
    }
    Ok(())
}

fn build_sequence(
    reference: &[u8],
    artificial_reference: &[u8],
    record: &bam::Record,
    offset: usize,
) -> Result<Vec<u8>> {
    let mut artificial_seq = Vec::new();
    let record_seq = record.seq();
    let mut record_pos = 0;
    let mut ref_pos = record.pos() as usize - offset;
    //Create random seq for leading softclips
    for cigar in record.cigar_cached().unwrap().iter() {
        match cigar.char() {
            'S' => {
                add_random_bases(cigar.len() as u64, &mut artificial_seq)?;
                //Add random Sequence of length cigar.len()
            }
            'M' => {
                (0..cigar.len()).for_each(|_| {
                    let ref_base = artificial_reference.get(ref_pos).unwrap();
                    if record_seq.encoded_base(record_pos) == *reference.get(ref_pos).unwrap() {
                        artificial_seq.push(*artificial_reference.get(ref_pos).unwrap());
                    } else {
                        add_random_base(&mut artificial_seq, ref_base);
                    }
                    ref_pos += 1;
                    record_pos += 1;
                });
                // Add reference bases except for mismatches
            }
            'I' => {
                add_random_bases(cigar.len() as u64, &mut artificial_seq)?;
                record_pos += cigar.len() as usize;
            }
            'D' | 'N' => {
                ref_pos += cigar.len() as usize;
            }
            'X' => {
                (0..cigar.len()).for_each(|_| {
                    let ref_base = artificial_reference.get(ref_pos).unwrap();
                    add_random_base(&mut artificial_seq, ref_base);
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

fn add_random_base(seq: &mut Vec<u8>, excluded_base: &u8) {
    let mut rng = rand::thread_rng();
    let mut alphabet = vec![b'A', b'C', b'G', b'T'];
    alphabet.retain(|x| x != excluded_base);
    let idx = rng.gen_range(0, 3);
    seq.push(alphabet[idx]);
}

fn add_random_bases(length: u64, seq: &mut Vec<u8>) -> Result<()> {
    let alphabet = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::thread_rng();
    (0..length).for_each(|_| {
        let idx = rng.gen_range(0, 4);
        seq.push(alphabet[idx])
    });
    Ok(())
}
