# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.21.0] - 2021-04-30
### Changed
- Small fixes for `rbt vcf-report`, e.g. for handling unexpected multiple canonical transcripts.

## [0.20.6] - 2021-XX-XX
### Changed
- Bugfix for `rbt vcf-report` that eliminates the potential risk of receiving plots without any reads overlaying the variant when using `--max-read-depth`.

## [0.20.5] - 2021-04-19
### Changed
- Bugfix for `rbt vcf-report` that stops displaying undefined values in the table-report.
### Added
- New parameter `--pin-until` for `rbt csv-report`. `--pin-until` pins the table until the given column such that scrolling to the right does not hide the given column and those before.

## [0.20.4] - 2021-04-15
### Changed
- Fixed a JS bug in VCF report leading to an error with empty annotation fields.

## [0.20.3] - 2021-03-31
### Changed
- Cosmetic changed and bug fixes for `rbt csv-report`.

## [0.20.2] - 2021-03-19
### Changed
- Cleaned up histogram presentation for pathologic columns.

## [0.20.1] - 2021-03-17
### Changed
- Added button to highlight numeric values in column histogram (`rbt csv-report`).

## [0.20.0] - 2021-03-10
### Added
- New parameter `--formatter` for `rbt csv-report`.
- `--formatter` allows the specification of custom formatting functions for one or multiple columns of the given csv file.

## [0.19.7] - 2021-03-09
### Changed
- Performance improvement for `rbt csv-report`.

## [0.19.6] - 2021-03-09
### Changed
- Fixed a JSON syntax error in `rbt vcf-report`.

## [0.19.5] - 2021-03-05
### Changed
- Fixed internal file naming in order to avoid too long filenames.

## [0.19.4] - 2021-03-04
### Changed
- Fixed further bugs in the csv-report.

## [0.19.3] - 2021-03-02
### Changed
- Improved fasta sequence length retrieval in `rbt vcf-report` (@fxwiegand).
- Better error messages (@fxwiegand).
- Add check for empty values in `rbt csv-report` (@fxwiegand).
- Removed unused dependencies (@fxwiegand).

## [0.19.2] - 2021-02-12
### Changed
- Various small bug fixes in `rbt vcf-report`.

## [0.19.1] - 2021-02-08
### Changed
- Fix for javascript import order problem in `rbt vcf-report`.

## [0.19.0] - 2021-02-05
### Added
- New subcommand csv-report that allows to generate an interactive HTML report from a CSV/TSV table.
### Changed
- Some polishing for vcf-report.

## [0.18.1] - 2021-01-27
### Changed
- If CANONICAL field is not present in ANN of bcf, rbt vcf-report now assumes that transcript is not the canonical one.

## [0.18.0] - 2021-01-25
### Added
- Parameter `--plot-info` for `rbt vcf-report` that allows to plot arbitrary info fields next to the oncoprint matrix.
- Ability to display intergenic variants at the primary stage of `rbt vcf-report`.

## [0.17.0] - 2021-01-21
### Changed
- `rbt collapse-reads-to-fragments` now always writes FASTQ files since consensus reads need to be remapped anyways (MAPQ might change).
- `rbt vcf-report` can now be parallelized.
- Various fixes and improvements for `rbt vcf-report`
- Various fixes for `rbt collapse-reads-to-fragments`

## [0.16.0] - 2020-11-17
### Changed
- `rbt collapse-reads-to-fragments bam` writes skipped reads to separate bam file now
- `rbt collapse-reads-to-fragments bam` does not perform starcode clustering anymore
- Fixed bug in `rbt collapse-reads-to-fragments bam` failing when read mates map to different chromosomes

## [0.15.1] - 2020-11-12
### Changed
- Fixed bug in vcf-split that led to unequally filled splitted VCF/BCFs.

## [0.15.0] - 2020-11-11
### Changed
- Allow INFO field prefixes in vcf-report.
- Require annotation with --hgvsg in vcf-report.
- Fixes for vcf-report layout.
### Added
- New subcommand vcf-split for splitting VCF/BCF files into N equal chunks, including BND support.

## [0.14.0] - 2020-10-26
### Added
- New flag `--tsv` for `rbt vcf-report`
- `--tsv` allows adding a custom tsv file that will be visualized in the vcf-report

## [0.13.0] - 2020-10-06
### Changed
- CLI syntax for `rbt vcf-report` has been changed in order to allow BCF/VCF files holding variants from multiple samples.

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
- Flag --sources that allows to constrain data sources for DGIdb annotation of VCF/BCFs.

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
- BCF output for `vcf-annotate-dgidb`

## [0.9.0] - 2020-01-12
### Added
- Subcommand `rbt vcf-fix-iupac-alleles` for converting IUPAC codes in VCF/BCF files into Ns.

## [0.8.3] - 2019-11-08
### Changed
- Fixed a bug in `annotate-dgidb` where empty VCF input-files resulted in an error.

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
  that were merged into this consensus read, the `--verbose-read-names` flag can be used.


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
- Handle missing values in `vcf-to-txt`.
- Bug fixes.

## [0.2.6] - 2018-07-04
### Added
- Support for Number="R" tag in VCF files.

## [0.2.5] - 2018-04-06
### Changed
- Fixed various small bugs.

## [0.2.0] - 2018-02-23
### Added
- A tool to selectively remove records from a FASTQ file.
- A tool to calculate b-allele frequencies (BAF) from a VCF file.

## [0.1.3] - 2017-12-05
### Changed
- Adapted to changes in Rust-Htslib.

## [0.1.2] - 2016-11-02
### Changed
- Improved robustness for exotic variant calls in `vcf-match` and `vcf-to-txt`.

## [0.1.1] - 2016-08-31
### Changed
- Hotfix of a bug in `vcf-match`.


## [0.1.0] - 2016-08-31
### Added
- Tool for matching VCF/BCF files.
- Tool for converting VCF/BCF files to tabular TXT.
- Tool for calculating BAM depth at given loci.
- Tool for splitting FASTQ into n chunks in a single pass using a round-robin procedure.
