use std::io;
use std::cmp;
use std::error::Error;

use csv;

use rust_htslib::bam;
use rust_htslib::bam::Read;


#[derive(RustcDecodable, Debug)]
struct PosRecord {
    chrom: String,
    pos: u32
}


pub fn depth(
    bam_path: &str,
    max_read_length: u32,
    include_flags: u16,
    exclude_flags: u16,
    min_mapq: u8) -> Result<(), Box<Error>>
{
    let mut bam_reader = bam::IndexedReader::from_path(&bam_path)?;
    let mut pos_reader = csv::Reader::from_reader(io::stdin()).has_headers(false).delimiter(b'\t');
    let mut csv_writer = csv::Writer::from_buffer(io::BufWriter::new(io::stdout())).delimiter(b'\t');
    let mut last_tid = 0;
    let mut last_pos = 0;
    let mut exceeded = false;
    let mut pileup_iter = bam_reader.pileup();

    for (i, record) in pos_reader.decode().enumerate() {
        let record: PosRecord = try!(record);

        // jump to correct position if necessary
        let tid = bam_reader.header.tid(record.chrom.as_bytes()).unwrap();
        let start = cmp::max(record.pos as i32 - max_read_length as i32 - 1, 0) as u32;
        if exceeded || tid != last_tid || start > last_pos || record.pos - 1 <= last_pos {
            let n = bam_reader.header.target_len(tid).unwrap();
            try!(bam_reader.fetch(tid, start, n));
            pileup_iter = bam_reader.pileup();
        }

        // iterate over pileups
        let mut covered = false;
        exceeded = false;
        while let Some(pileup) = pileup_iter.next() {
            let pileup = try!(pileup);
            covered = pileup.pos() == record.pos - 1;
            last_tid = pileup.tid();
            last_pos = pileup.pos();

            if covered {
                let depth = pileup.alignments().filter(|alignment| {
                    let record = alignment.record();
                    let flags = record.flags();
                    (!flags) & include_flags == 0 &&
                    flags & exclude_flags == 0 &&
                    record.mapq() >= min_mapq
                }).count();

                try!(csv_writer.encode((&record.chrom, record.pos, depth)));
                break;
            } else if pileup.pos() > record.pos {
                exceeded = true;
                break;
            }
        }
        if !covered {
            try!(csv_writer.encode((&record.chrom, record.pos, 0)));
        }

        if (i + 1) % 100 == 0 {
            info!("{} records written.", i + 1);
        }
    }
    Ok(())
}
