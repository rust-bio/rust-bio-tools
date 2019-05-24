use snafu::Snafu;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Snafu)]
#[snafu(visibility(pub))]
pub enum Error {
    #[snafu(display("Could not open input file {}: {:?}", filename, source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    ReaderError {
        filename: String,
        source: std::io::Error,
    },
    #[snafu(display("Could not open output file {}: {:?}", filename, source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    WriterError {
        filename: String,
        source: std::io::Error,
    },
    #[snafu(display("Pipeline Error with {}: {}", params, source))]
    #[snafu(source(from((dyn std::error::Error + 'static), Box::new)))]
    PipelineError {
        params: String,
        source: Box<std::error::Error>,
    },
}
