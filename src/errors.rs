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
    }

    
}
