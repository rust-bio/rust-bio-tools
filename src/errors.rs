use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Error, Debug)]
pub enum Error {
    #[error("Could not open input file; {}", source)]
    BamIndexedReaderError { source: rust_htslib::errors::Error },

    #[error("Could not open input file; {}", source)]
    BamReaderError { source: rust_htslib::errors::Error },

    #[error("Could not open input file {}; {}", filename, source)]
    FastqReaderError {
        filename: String,
        source: std::io::Error,
    },

    #[error("Could not open output file; {}", source)]
    BamWriterError { source: rust_htslib::errors::Error },

    #[error("Could not open output file {}; {}", filename, source)]
    FastqWriterError {
        filename: String,
        source: std::io::Error,
    },

    #[error("Could not write record; {}", source)]
    BamWriteError { source: rust_htslib::errors::Error },

    #[error("Could not write FASTQ record; {}", source)]
    FastqWriteError { source: std::io::Error },

    #[error("Could not read record with index {}; {}", record_idx, source)]
    BamReadError {
        record_idx: usize,
        source: rust_htslib::errors::Error,
    },

    #[error("Could not read csv record; {}", source)]
    CsvReadError { source: csv::Error },

    #[error("Could not write to csv record; {}", source)]
    CsvWriteError { source: csv::Error },

    #[error("Could not read record with index {}; {}", record_idx, source)]
    FastqReadError {
        record_idx: usize,
        source: std::io::Error,
    },

    // TODO rewrite this to make the forward reverse part of display a function
    // format!("Given FASTQ files have unequal lengths. Reverse file returned record {} as empty, forward record is not: id:'{}' seq:'{:?}'.", i, f_rec.id(), str::from_utf8(f_rec.seq()));
    #[error("Given FASTQ files have unequal lengths. Forward file returned record {} as empty, reverse record is not: id:'{}' seq:'{:?}'.", nr, name, seq)]
    OrphanPairedEndReadError {
        nr: usize,
        name: String,
        seq: String,
        forward_orphan: bool,
    },

    #[error("Failed to run UMI clustering with starcode: {:?}", source)]
    StarcodeCallError { source: std::io::Error },

    #[error("Failed to write {} to starcode via stdin: {:?}", payload, source)]
    StarcodeWriteError {
        payload: String,
        source: std::io::Error,
    },

    #[error(
        "Error parsing the following Starcode cluster in csv format {:?}: {}",
        record,
        source
    )]
    StarcodeClusterParseError {
        record: csv::StringRecord,
        source: csv::Error,
    },

    #[error("Error parsing csv record: {}", source)]
    CsvRecordParseError { source: csv::Error },

    #[error("Failed to flush pipeline: {:?}", source)]
    WriteFlushError { source: std::io::Error },

    #[error("Pipeline Error with {}: {}", params, source)]
    PipelineError {
        params: String,
        source: Box<dyn std::error::Error>,
    },

    #[error("Failed to create a temporary directory for RocksDB: {:?}", source)]
    TempdirCreationError { source: std::io::Error },

    #[error("Failed to create a RocksDB at location {}: {:?}", filename, source)]
    FastqStorageCreationError {
        filename: String,
        source: rocksdb::Error,
    },

    #[error("Failed to put read nr {} into the rocksDB; {:?}", read_nr, source)]
    FastqStoragePutError {
        read_nr: usize,
        source: rocksdb::Error,
    },

    #[error("Failed to get read nr {} from the rocksDB; {:?}", read_nr, source)]
    FastqStorageGetError {
        read_nr: usize,
        source: rocksdb::Error,
    },

    // TODO Would we like to also return the read name here?
    // If so, forward, reverse, or both?
    #[error("Serde failed to encode read pair nr {}; {:?}", read_nr, source)]
    ReadSerializationError {
        read_nr: usize,
        source: serde_json::error::Error,
    },

    // TODO Would we like to also return the read name here?
    // If so, forward, reverse, or both?
    #[error("Serde failed to decode read pair nr {}; {:?}", read_nr, source)]
    ReadDeserializationError {
        read_nr: usize,
        source: serde_json::error::Error,
    },

    #[error("Cigar {} contains unsupported Cigar operation; {}", cigar, source)]
    BamCigarError {
        cigar: String,
        source: rust_htslib::errors::Error,
    },

    #[error("Could not fetch record {}; {}", tid, source)]
    BamFetchError {
        tid: u32,
        source: rust_htslib::errors::Error,
    },

    #[error("Could not read pileup from {}; {}", path, source)]
    BamPileupError {
        path: String,
        source: rust_htslib::errors::Error,
    },

    #[error("Failed to read BCF file from stdin; {}", source)]
    BCFReaderStdinError { source: rust_htslib::errors::Error },

    #[error(
        "Failed to open BCF writer for stdout with header {:?}; {}",
        header,
        source
    )]
    BCFWriterStdoutError {
        header: String,
        source: rust_htslib::errors::Error,
    },

    #[error("Failed to format data in BCF record; {}", source)]
    BCFFormatReadError { source: rust_htslib::errors::Error },

    #[error("Failed to format {:?} as {:?}; {}", data, fd, source)]
    BCFTagWriteError {
        data: String,
        fd: String,
        source: rust_htslib::errors::Error,
    },

    #[error("Failed write record to BCF file; {}", source)]
    BCFWriteError { source: rust_htslib::errors::Error },

    #[error("Failed to read BCF file with header {}; {}", header, source)]
    BCFReadError {
        header: String,
        source: rust_htslib::errors::Error,
    },

    #[error(
        "Failed to find record id {} in BCF with header {}; {}",
        rid,
        header,
        source
    )]
    BCFReadIdError {
        rid: u32,
        header: String,
        source: rust_htslib::errors::Error,
    },

    #[error("Failed to create new BCF reader from {}; {}", path, source)]
    BCFReaderFromPathError {
        path: String,
        source: rust_htslib::errors::Error,
    },

    //TODO This should be more clear
    #[error("Missing tag {}", tag)]
    MissingTagError { tag: String },

    //TODO This should be more clear
    #[error("Unsupported variant")]
    UnsupportedVariantError,

    #[error("{:?}", source)]
    StdStrUtf8Error { source: std::str::Utf8Error },

    #[error("Failed to create new BCF reader from {}; {}", path, source)]
    FileOpenError {
        path: String,
        source: std::io::Error,
    },

    #[error("Failed to parse tag length: {}. Currently, only R, A and 1 are supported multiplicities of tags", tag_length)]
    UnsupportedTagLengthError { tag_length: String },

    #[error("Failed to write to BCF writer; {}", source)]
    BCFStdIOWriteError { source: std::io::Error },

    #[error(
        "Failed to read info from BCF record with index {}; {}",
        record_idx,
        source
    )]
    BCFInfoReadError {
        record_idx: usize,
        source: rust_htslib::errors::Error,
    },

    #[error(
        "Invalid combination of files. Each pair of files \
         (input and output) need to be both gzipped or \
         both not zipped."
    )]
    InvalidFileCombinationError,

    #[error("There is no flag type for format")]
    FlagTypeError,
}
