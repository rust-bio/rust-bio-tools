# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).


## [0.3.0] - 2019-01-??
### Added
- A tool to merge FASTQ reads using unique molecular identifiers (UMIs).


## [0.2.7] - 2018-12-28
### Changed
- Handle missing values in vcf-to-txt.
- Bug fixes.

## [0.2.6] - 2018-07-04
### Added
- Support for Number="R" tag in VCF files.

## [0.2.5] - 2018-04-06
### Changed
- Fixed various small bugs.

## [0.2.0] - 2018-02-23
### Added
- A tool to selectively remove records from a fastq file. 
- A tool to calculate b-allele frequencies (BAF) from a vcf file.

## [0.1.3] - 2017-12-05
### Changed
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
