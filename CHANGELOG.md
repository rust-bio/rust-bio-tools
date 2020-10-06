# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.13.0] - 2020-10-06
### Changed
- cli syntax for `rbt vcf-report` has been changed in order to allow BCF/VCF files holding variants from multiple samples

## [0.12.2] - 2020-09-30
### Added
- New flags `--info`, `--js` and `--format` for `rbt vcf-report`
- `--js` allows definition of a custom javascript template for creating vcf-reports
- `--info` and `--format` allow specifying INFO- and FORMAT-fields that can be used by the custom javascript template

## [0.12.1] - 2020-09-15
### Added
- `rbt vcf-report` allows to create an interactive HTML report allowing advanced analysis of variants from BCF/VCF and BAM files.

### Removed
- `rbt oncoprint` is now included in `rbt vcf-report` and therefore has been removed.

## [0.11.0] - 2020-08-05
### Added
- Flag --sources that allows to constrain data sources for dgidb annotation of VCF/BCFs.

## [0.10.1] - 2020-05-13
### Changed
- Only print recurrent genes in oncoprint, in order to save time and space.

## [0.10.0] - 2020-04-22
### Added
- `rbt oncoprint`, which allows to create an interactive HTML oncoprint from BCF/VCF files.

### Changed
- `rbt call-consensus-reads` has been renamed to `rbt collapse-reads-to-fragments`, which better highlights what is actually done.

## [0.9.2] - 2020-03-20
### Changed
- BCF output for vcf-annotate-dgidb

## [0.9.0] - 2020-01-12
### Added
- Subcommand vcf-fix-iupac-alleles for converting IUPAC codes in VCF/BCF files into Ns.

## [0.8.3] - 2019-11-08
### Changed
- Fixed a bug in `annotate-dgidb` where empty vcf input-files resulted in an error.

## [0.8.0] - 2019-10-15
### Added
- Annotation of VCF/BCF files with drug interactions from DGIdb.

## [0.6.0] - 2019-08-28
### Added
- Consensus reads can now be calculated from BAM files, based on marked PCR/optical duplicates from picard tools.

## [0.5.0] - 2019-05-08
### Changed
- Fixed a bug in `call-consensus-reads` where read sequences were kept in memory.
  This significantly reduces memory footprint.
- `call-consensus-reads` now writes shorter name lines for consensus reads.
  Name lines of consensus reads now only contain a random uuid and the number of
  reads they were created from. To get the old behavior, i.e. a list of reads
  that were merged into this consensus read,http://semver.org/ the `--verbose-read-names` can be used.


## [0.4.1] - 2019-05-02
### Changed
- Fixed a bug in `call-consensus-reads` where BGZF-compressed files terminated early (after the first block).
- Extended documentation of CLI.
- Updated used versions of rust-htslib and rocksdb.
- Renamed the `--reverse-umi` parameter of `call-consensus-reads` to `--umi-on-reverse`.

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
