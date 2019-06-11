use snafu::Snafu;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Snafu)]
#[snafu(visibility(pub))]
pub enum Error {
    #[snafu(display("Could not open input file from path {}: {:?}", filepath, source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    BamIndexedReaderError {
        filepath: String,
        source: rust_htslib::bam::IndexedReaderPathError,
    },

    #[snafu(display("Could not open input file {}: {:?}", filename, source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    BamReaderError {
        filename: String,
        source: rust_htslib::bam::ReaderPathError,
    },

    #[snafu(display("Could not open input file {}: {:?}", filename, source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    FastqReaderError {
        filename: String,
        source: std::io::Error,
    },

    #[snafu(display("Could not open output file {}: {:?}", filename, source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    BamWriterError {
        filename: String,
        source: rust_htslib::bam::WriterPathError,
    },

    #[snafu(display("Could not open output file {}: {:?}", filename, source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    FastqWriterError {
        filename: String,
        source: std::io::Error,
    },

    #[snafu(display("Could not write record {:?}: {:?}", record, source))]
    BamWriteError {
        record: Option<rust_htslib::bam::record::Record>,
        source: rust_htslib::bam::WriteError,
    },

    #[snafu(display("Could not write record {:?}: {:?}", record, source))]
    FastqWriteError {
        record: Option<bio::io::fastq::Record>,
        source: std::io::Error,
    },

    #[snafu(display("Could not read record: {:?}", source))]
    BamReadError { source: rust_htslib::bam::ReadError },

    #[snafu(display("Could not read csv record: {:?}", source))]
    CsvReadError { source: csv::Error },

    #[snafu(display("Could not write to csv record: {:?}", source))]
    CsvWriteError { source: csv::Error },

    #[snafu(display("Could not read record: {:?}", source))]
    FastqReadError { source: std::io::Error },

    // TODO rewrite this to make the forward reverse part of display a function
    // format!("Given FASTQ files have unequal lengths. Reverse file returned record {} as empty, forward record is not: id:'{}' seq:'{:?}'.", i, f_rec.id(), str::from_utf8(f_rec.seq()));
    #[snafu(display("Given FASTQ files have unequal lengths. Forward file returned record {} as empty, reverse record is not: id:'{}' seq:'{:?}'.", nr, name, seq))]
    OrphanPairedEndReadError {
        nr: usize,
        name: String,
        seq: String,
        forward_orphan: bool,
    },

    #[snafu(display("Failed to run UMI clustering with starcode: {:?}", source))]
    StarcodeCallError { source: std::io::Error },

    #[snafu(display("Failed to write {} to starcode via stdin: {:?}", payload, source))]
    StarcodeWriteError {
        payload: String,
        source: std::io::Error,
    },

    #[snafu(display(
        "Error parsing the following Starcode cluster in csv format {:?}: {}",
        record,
        source
    ))]
    StarcodeClusterParseError {
        record: csv::StringRecord,
        source: csv::Error,
    },

    #[snafu(display("Error parsing csv record: {}", source))]
    CsvRecordParseError { source: csv::Error },

    #[snafu(display("Failed to flush pipeline: {:?}", source))]
    WriteFlushError { source: std::io::Error },

    #[snafu(display("Pipeline Error with {}: {}", params, source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    PipelineError {
        params: String,
        source: Box<std::error::Error>,
    },

    #[snafu(display("Failed to create a temporary directory for RocksDB: {:?}", source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    TempdirCreationError { source: std::io::Error },

    #[snafu(display("Failed to create a RocksDB at location {}: {:?}", filename, source))]
    FastqStorageCreationError {
        filename: String,
        source: rocksdb::Error,
    },

    #[snafu(display("Failed to put read nr {} into the rocksDB.: {:?}", read_nr, source))]
    FastqStoragePutError {
        read_nr: usize,
        source: rocksdb::Error,
    },

    #[snafu(display("Failed to get read nr {} from the rocksDB.: {:?}", read_nr, source))]
    FastqStorageGetError {
        read_nr: usize,
        source: rocksdb::Error,
    },

    // TODO Would we like to also return the read name here?
    // If so, forward, reverse, or both?
    #[snafu(display("Serde failed to encode read pair nr {}.: {:?}", read_nr, source))]
    ReadSerializationError {
        read_nr: usize,
        source: serde_json::error::Error,
    },

    // TODO Would we like to also return the read name here?
    // If so, forward, reverse, or both?
    #[snafu(display("Serde failed to decode read pair nr {}.: {:?}", read_nr, source))]
    ReadDeserializationError {
        read_nr: usize,
        source: serde_json::error::Error,
    },

    #[snafu(display(
        "Record {:?} contains unsupported Cigar operation.: {:?}",
        record,
        source
    ))]
    BamCigarError {
        record: rust_htslib::bam::record::Record,
        source: rust_htslib::bam::record::CigarError,
    },

    #[snafu(display("Could not fetch record {}: {:?}", tid, source))]
    BamFetchError {
        tid: u32,
        source: rust_htslib::bam::FetchError,
    },

    #[snafu(display("Could not read pileup from {}: {:?}", path, source))]
    BamPileupError {
        path: String,
        source: rust_htslib::bam::pileup::PileupError,
    },

    #[snafu(display("Failed to read BCF file from stdin: {:?}", source))]
    BCFReaderStdinError { source: rust_htslib::bcf::BCFError },

    #[snafu(display(
        "Failed to open BCF writer for stdout with header {:?}: {:?}",
        header,
        source
    ))]
    BCFWriterStdoutError {
        header: String,
        source: rust_htslib::bcf::BCFError,
    },

    #[snafu(display(
        "Failed to interpret data as integer with format {:?} in BCF record: {:?}",
        fd,
        source
    ))]
    BCFFormatReadError {
        fd: String,
        source: rust_htslib::bcf::record::FormatReadError,
    },

    #[snafu(display("Failed to format {:?} as {:?}: {:?}", bafs, fd, source))]
    BCFTagWriteError {
        bafs: Vec<f32>,
        fd: String,
        source: rust_htslib::bcf::record::TagWriteError,
    },

    #[snafu(display("Failed write record {} to BCF file: {:?}", record, source))]
    BCFWriteError {
        record: String,
        source: rust_htslib::bcf::WriteError,
    },

    #[snafu(display("Failed to read BCF file with header {}: {:?}", header, source))]
    BCFReadError {
        header: String,
        source: rust_htslib::bcf::ReadError,
    },
}
