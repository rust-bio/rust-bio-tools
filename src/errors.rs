use snafu::Snafu;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Snafu)]
#[snafu(visibility(pub))]
pub enum Error {
    #[snafu(display("Could not open input file {}: {:?}", filename, source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    FastqReaderError {
        filename: String,
        source: std::io::Error,
    },
    #[snafu(display("Could not open output file {}: {:?}", filename, source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    FastqWriterError {
        filename: String,
        source: std::io::Error,
    },
    #[snafu(display("Could not write record {:?}: {:?}", record, source))]
    FastqWriteError {
        record: Option<bio::io::fastq::Record>,
        source: std::io::Error,
    },
    #[snafu(display("Pipeline Error with {}: {}", params, source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    PipelineError {
        params: String,
        source: Box<std::error::Error>,
    },

    #[snafu(display("Error parsing the following Starcode cluster in csv format {:?}: {}", record, source))]
    StarcodeClusterParseError {
        record: csv::StringRecord,
        source: csv::Error,
    },

    #[snafu(display("Failed to create a temporary directory for RocksDB: {:?}", source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    TempdirCreationError{
        source: std::io::Error,
    },

    #[snafu(display("Failed to create a RocksDB at location {}: {:?}", filename, source))]
    FastqStorageCreationError{
        filename: String,
        source: rocksdb::Error,
    },

    #[snafu(display("Failed to put read nr {} into the rocksDB.: {:?}", read_nr, source))]
    FastqStoragePutError{
        read_nr: usize,
        source: rocksdb::Error,
    },

    #[snafu(display("Failed to get read nr {} from the rocksDB.: {:?}", read_nr, source))]
    FastqStorageGetError{
        read_nr: usize,
        source: rocksdb::Error,
    },

    // TODO Would we like to also return the read name here?
    // If so, forward, reverse, or both?
    #[snafu(display("Serde failed to encode read pair nr {}.: {:?}", read_nr, source))]
    ReadSerializationError{
        read_nr: usize,
        source: serde_json::error::Error,
    },

    // TODO Would we like to also return the read name here?
    // If so, forward, reverse, or both?
    #[snafu(display("Serde failed to decode read pair nr {}.: {:?}", read_nr, source))]
    ReadDeserializationError{
        read_nr: usize,
        source: serde_json::error::Error,
    },
}
