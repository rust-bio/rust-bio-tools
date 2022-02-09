# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

### [0.38.3](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.38.2...v0.38.3) (2022-02-09)


### Bug Fixes

* Fix vcf-report allele frequency heatmap width ([#236](https://www.github.com/rust-bio/rust-bio-tools/issues/236)) ([6f119d3](https://www.github.com/rust-bio/rust-bio-tools/commit/6f119d3bab34594029a846172e29623fb9ed71aa))
* Fix vcf-report read plot height when using --max-read-depth ([#235](https://www.github.com/rust-bio/rust-bio-tools/issues/235)) ([4350a1a](https://www.github.com/rust-bio/rust-bio-tools/commit/4350a1a3887d7f22e16de8ed8944981d138982d6))

### [0.38.2](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.38.1...v0.38.2) (2022-01-14)


### Bug Fixes

* add popper in csv report ([#233](https://www.github.com/rust-bio/rust-bio-tools/issues/233)) ([a095755](https://www.github.com/rust-bio/rust-bio-tools/commit/a0957551e39928a510c7abd408378c274a3b5c33))

### [0.38.1](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.38.0...v0.38.1) (2022-01-06)


### Bug Fixes

* Report 1-based position in vcf-report log messages ([#231](https://www.github.com/rust-bio/rust-bio-tools/issues/231)) ([bcf6260](https://www.github.com/rust-bio/rust-bio-tools/commit/bcf62609b52578a30a21a42971e0001701941681))

## [0.38.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.37.0...v0.38.0) (2021-12-25)


### Features

* Add aux field to tooltip ([#226](https://www.github.com/rust-bio/rust-bio-tools/issues/226)) ([cfc8fe5](https://www.github.com/rust-bio/rust-bio-tools/commit/cfc8fe572fbc2e498d3d3cc0fe0de9e320b773d1))


### Bug Fixes

* Apply latest clippy lints ([#227](https://www.github.com/rust-bio/rust-bio-tools/issues/227)) ([d1a60cc](https://www.github.com/rust-bio/rust-bio-tools/commit/d1a60cc8876eb8ea2b6cffb21bdcee44a00d2dc4))

## [0.37.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.36.0...v0.37.0) (2021-12-19)


### Features

* anonymize bam - adapt mandatory and aux fields ([#219](https://www.github.com/rust-bio/rust-bio-tools/issues/219)) ([e6ffdea](https://www.github.com/rust-bio/rust-bio-tools/commit/e6ffdea140704e305d8b30213c57d76bf7841fc2))

## [0.36.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.35.0...v0.36.0) (2021-12-15)


### Features

* Use Feature before HGVSg for level 1 naming ([#224](https://www.github.com/rust-bio/rust-bio-tools/issues/224)) ([4377e42](https://www.github.com/rust-bio/rust-bio-tools/commit/4377e428e9bca2b806e666a8234da35e027fd523))
* Warn instead of panic for missing HGVSg values in vcf-report ([#220](https://www.github.com/rust-bio/rust-bio-tools/issues/220)) ([5f16505](https://www.github.com/rust-bio/rust-bio-tools/commit/5f1650506fb502548e07d2756ff99557bb562580))


### Bug Fixes

* Fix possible broken link to third stage of vcf-report ([#221](https://www.github.com/rust-bio/rust-bio-tools/issues/221)) ([66299d2](https://www.github.com/rust-bio/rust-bio-tools/commit/66299d2e0e5ac96f67db8fcf4fef7834c2a2de05))

## [0.35.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.34.0...v0.35.0) (2021-12-14)


### Features

* Add --annotation-field option to vcf-report ([#211](https://www.github.com/rust-bio/rust-bio-tools/issues/211)) ([29719d5](https://www.github.com/rust-bio/rust-bio-tools/commit/29719d51f0169de39570c6715af5b45336691ebb))


### Bug Fixes

* add missing js-file to vcf-report ([#212](https://www.github.com/rust-bio/rust-bio-tools/issues/212)) ([55b3166](https://www.github.com/rust-bio/rust-bio-tools/commit/55b3166516b6f369577df74f8cc020fbe9629721))
* release version ([5de5480](https://www.github.com/rust-bio/rust-bio-tools/commit/5de548073816dc46ec68ad352c3b91a738971e20))
* vcf-report tests ([#215](https://www.github.com/rust-bio/rust-bio-tools/issues/215)) ([b3e45b3](https://www.github.com/rust-bio/rust-bio-tools/commit/b3e45b3a290bfea7d1553ff34e186f1f870c6e02))


### Miscellaneous Chores

* release 0.33.2 ([10d2446](https://www.github.com/rust-bio/rust-bio-tools/commit/10d24469bfc5ea72ea385e82537e1fb52ca8f364))

## [0.34.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.33.1...v0.34.0) (2021-12-11)

### Features

* Add --annotation-field option to vcf-report ([#211](https://github.com/rust-bio/rust-bio-tools/pull/211)) ([29719d5](https://github.com/rust-bio/rust-bio-tools/commit/29719d51f0169de39570c6715af5b45336691ebb))


### Bug Fixes

* add missing js-file to vcf-report ([#212](https://github.com/rust-bio/rust-bio-tools/pull/212)) ([55b3166](https://github.com/rust-bio/rust-bio-tools/commit/55b3166516b6f369577df74f8cc020fbe9629721))

## [0.33.1](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.33.0...v0.33.1) (2021-12-07)


### Features

* Improve gene search for vcf-report ([#209](https://www.github.com/rust-bio/rust-bio-tools/issues/209)) ([cbeaff5](https://www.github.com/rust-bio/rust-bio-tools/commit/cbeaff5a79bf19fb745fa0e15505d5c41a9a6dd5))

## [0.33.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.32.0...v0.33.0) (2021-12-06)


### Features

* Updates and fixes for vcf-report ([#207](https://www.github.com/rust-bio/rust-bio-tools/issues/207)) ([e19b4fb](https://www.github.com/rust-bio/rust-bio-tools/commit/e19b4fb4e3d935bfa5f41796dc23669a903563dc))


### Bug Fixes

* Add error messages for missing annotations ([#206](https://www.github.com/rust-bio/rust-bio-tools/issues/206)) ([6210dd2](https://www.github.com/rust-bio/rust-bio-tools/commit/6210dd22749c4cff818724f1a24dc6d8f104c927))

## [0.32.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.31.0...v0.32.0) (2021-12-01)


### Features

* Update rust-bio ([#202](https://www.github.com/rust-bio/rust-bio-tools/issues/202)) ([707d0ca](https://www.github.com/rust-bio/rust-bio-tools/commit/707d0ca294c56f79c770a8d8f7ff530642868b39))


### Bug Fixes

* Fix read processing when building consensus reads ([#199](https://www.github.com/rust-bio/rust-bio-tools/issues/199)) ([3fdcb3f](https://www.github.com/rust-bio/rust-bio-tools/commit/3fdcb3ffd586db94e40c182e980110334280962f))

## [0.31.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.30.2...v0.31.0) (2021-11-16)


### Features

* rbt plot-bam: show position of reference base in tooltip ([#196](https://www.github.com/rust-bio/rust-bio-tools/issues/196)) ([3711c09](https://www.github.com/rust-bio/rust-bio-tools/commit/3711c0963b1fb95a341449aba8fe259b9faafcea))


### Bug Fixes

* rbt plot-bam: Disconnect neighboring insertion markers ([#197](https://www.github.com/rust-bio/rust-bio-tools/issues/197)) ([d322e6a](https://www.github.com/rust-bio/rust-bio-tools/commit/d322e6a9a64b870fe60b9ba21c4ec1b8e6c6c67a))

### [0.30.2](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.30.1...v0.30.2) (2021-11-09)


### Bug Fixes

* various fixes for consensus calling ([#191](https://www.github.com/rust-bio/rust-bio-tools/issues/191)) ([997d920](https://www.github.com/rust-bio/rust-bio-tools/commit/997d92087542f60de7818a4315266141d180b7ef))

### [0.30.1](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.30.0...v0.30.1) (2021-10-22)


### Bug Fixes

* Refine MAPQ color scheme and scale type ([#189](https://www.github.com/rust-bio/rust-bio-tools/issues/189)) ([6f47c68](https://www.github.com/rust-bio/rust-bio-tools/commit/6f47c685716e50501b13ae7d2ec9a23c8973f67b))

## [0.30.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.29.0...v0.30.0) (2021-10-22)


### Features

* Various alignment plot and CLI updates for vcf-report and plot-bam ([#187](https://www.github.com/rust-bio/rust-bio-tools/issues/187)) ([8739ef1](https://www.github.com/rust-bio/rust-bio-tools/commit/8739ef18acc8c357522c4619c294298c4f32091f))

## [0.29.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.28.0...v0.29.0) (2021-10-04)


### Features

* anonymize bam files ([#181](https://www.github.com/rust-bio/rust-bio-tools/issues/181)) ([22b884f](https://www.github.com/rust-bio/rust-bio-tools/commit/22b884fdf5cc8e617ad5c1657c85bff3429c128e))

## [0.28.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.27.0...v0.28.0) (2021-09-20)


### Features

* Revise third stage of vcf-report ([#182](https://www.github.com/rust-bio/rust-bio-tools/issues/182)) ([e34a466](https://www.github.com/rust-bio/rust-bio-tools/commit/e34a466f93058310b9ad6308af69e8369b980870))

## [0.27.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.26.1...v0.27.0) (2021-08-25)


### Features

* increase base quality scale  ([#179](https://www.github.com/rust-bio/rust-bio-tools/issues/179)) ([aba3a49](https://www.github.com/rust-bio/rust-bio-tools/commit/aba3a49c016f88e3688172e5e0bd387a5c30eeb9))


### Bug Fixes

* remove custom classification of duplicates and inversions from vcf-report

### [0.26.1](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.26.0...v0.26.1) (2021-08-19)


### Bug Fixes

* Fix handling of softclips in cigars
* Respect target CPU compat via rocksdb ^0.17 ([#175](https://www.github.com/rust-bio/rust-bio-tools/issues/175)) ([b2428f5](https://www.github.com/rust-bio/rust-bio-tools/commit/b2428f517a3ecd4a85596676f0adfc49ec24d268))

## [0.26.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.25.0...v0.26.0) (2021-08-13)


### Features

* Add tooltip for vcf-report bar plots ([#173](https://www.github.com/rust-bio/rust-bio-tools/issues/173)) ([fb91211](https://www.github.com/rust-bio/rust-bio-tools/commit/fb912112ba6ca1ccdfa5bc6be4fd741e46dd958b))

## [0.25.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.24.1...v0.25.0) (2021-08-09)


### Features

* add linkout for COSV entries to rbt vcf-report. ([224a076](https://www.github.com/rust-bio/rust-bio-tools/commit/224a0767835ead9c50ad3f0ab5c17f26aa8e9276))

### [0.24.1](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.24.0...v0.24.1) (2021-07-16)


### Bug Fixes

* Adjust FASTA read interval ([#167](https://www.github.com/rust-bio/rust-bio-tools/issues/167)) ([eccb3f6](https://www.github.com/rust-bio/rust-bio-tools/commit/eccb3f65491b1eef508359bfb08e4a6aea714451))

## [0.24.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.23.0...v0.24.0) (2021-07-16)


### Features

* Add ellipsis to csv report ([#166](https://www.github.com/rust-bio/rust-bio-tools/issues/166)) ([08ccad5](https://www.github.com/rust-bio/rust-bio-tools/commit/08ccad5491feca51f0b5c33da3380d92dca10683))
* Mark repeats in reference with lower opacity ([#163](https://www.github.com/rust-bio/rust-bio-tools/issues/163)) ([07235cf](https://www.github.com/rust-bio/rust-bio-tools/commit/07235cfec2bfe8553c735c2ab06d0813a8dbdfc2))


### Bug Fixes

* Missing lowercase reference bases in reports ([#160](https://www.github.com/rust-bio/rust-bio-tools/issues/160)) ([a661da1](https://www.github.com/rust-bio/rust-bio-tools/commit/a661da1448f1977a08e9b69f512a76bfb2337d7e))

## [0.23.0](https://www.github.com/rust-bio/rust-bio-tools/compare/v0.22.0...v0.23.0) (2021-07-07)


### Features

* New subcommand for plotting bam files ([#156](https://www.github.com/rust-bio/rust-bio-tools/issues/156)) ([707ef77](https://www.github.com/rust-bio/rust-bio-tools/commit/707ef772e08eb951ef631f54748bcf372994fb21))

## [0.22.0] - 2021-06-08
### Changed
- Revised `rbt collapse-reads-to-fragments bam`, CIGAR strings will now considered when merging reads  (@FelixMoelder).

## [0.21.1] - 2021-05-20
### Changed
- Further fixes for `rbt vcf-report` (@fxwiegand).

## [0.21.0] - 2021-04-30
### Changed
- Small fixes for `rbt vcf-report`, e.g. for handling unexpected multiple canonical transcripts and eliminating the potential risk of receiving plots without any reads overlaying the variant when using `--max-read-depth`.

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
