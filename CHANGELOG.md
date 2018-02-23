# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.2.0] - 2018-02-23
## Added
- A tool to selectively remove records from a fastq file. 
- A tool to calculate b-allele frequencies (BAF) from a vcf file.

## [0.1.3] - 2017-12-05
## Changed
- Adapted to changes in Rust-Htslib.

## [0.1.2] - 2016-11-02
### Changed
- Improved robustness for exotic variant calls in vcf-match and vcf-to-txt.

## [0.1.1] - 2016-08-31
### Changed
- Hotfix of a bug in vcf-match.


## [0.1.0] - 2016-08-31
### Added
- Tool for matching VCF/BCF files.
- Tool for converting VCF/BCF files to tabular TXT.
- Tool for calculating BAM depth at given loci.
- Tool for splitting FASTQ into n chunks in a single pass using a round-robin procedure.
