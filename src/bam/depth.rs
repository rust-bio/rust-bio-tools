use std::io;
use std::cmp;

use csv;

use rust_htslib::bam;
use rust_htslib::bam::Read;


#[derive(RustcDecodable, Debug)]
struct PosRecord {
    chrom: String,
    pos: u32
}


pub fn depth(bam_path: &str, max_read_length: u32) {
    let mut bam_reader = bam::IndexedReader::new(&bam_path).unwrap();
    let mut pos_reader = csv::Reader::from_reader(io::stdin()).has_headers(false).delimiter(b'\t');
    let mut csv_writer = csv::Writer::from_buffer(io::BufWriter::new(io::stdout())).delimiter(b'\t');

    for record in pos_reader.decode() {
        let record: PosRecord = record.unwrap();
        let tid = bam_reader.header.tid(record.chrom.as_bytes()).unwrap();
        bam_reader.seek(tid, cmp::max(record.pos as i32 - max_read_length as i32, 0) as u32, record.pos).unwrap();
        let mut covered = false;
        for pileup in bam_reader.pileup() {
            let pileup = pileup.unwrap();
            covered = pileup.pos() == record.pos;
            if covered {
                csv_writer.encode((&record.chrom, record.pos, pileup.depth())).unwrap();
                break;
            }
        }
        if !covered {
            csv_writer.encode((&record.chrom, record.pos, 0)).unwrap();
        }
    }
}
