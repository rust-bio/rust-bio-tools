# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.5.0] - 2019-05-08
### Changed
- Fixed a bug in `call-consensus-reads` where read sequences were kept in memory.
  This significantly reduces memory footprint.
- `call-consensus-reads` now writes shorter name lines for consensus reads.
  Name lines of consensus reads now only contain a random uuid and the number of
  reads they were created from. To get the old behavior, i.e. a list of reads
  that were merged into this consensus read, the `--verbose-read-names` can be used.


## [0.4.1] - 2019-05-02
### Changed
- Fixed a bug in `call-consensus-reads` where BGZF-compressed files terminated early (after the first block).
- Extended documentation of CLI.
- Updated used versions of rust-htslib and rocksdb.

## [0.4.0] - 2019-04-10
### Changed
- Consensus reads generated from UMI-tagged reads with `call-consensus-reads` no longer contain the UMI.

## [0.3.0] - 2019-01-24
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
