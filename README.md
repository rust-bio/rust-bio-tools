# Rust-Bio-Tools

[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/rust-bio/rust-bio-tools/CI)](https://github.com/rust-bio/rust-bio-tools/actions)

A set of ultra fast and robust command line utilities for bioinformatics tasks based on Rust-Bio.
Rust-Bio-Tools provides a command `rbt`, which currently supports the following operations:

* a linear time implementation for fuzzy matching of two vcf/bcf files (`rbt vcf-match`)
* a vcf/bcf to txt converter, that flexibly allows to select tags and properly handles multiallelic sites (`rbt vcf-to-txt`)
* a linear time round-robin FASTQ splitter that splits a given FASTQ files into a given number of chunks (`rbt fastq-split`)
* a linear time extraction of depth information from BAMs at given loci (`rbt bam-depth`)
* a utility to quickly filter records from a FASTQ file (`rbt fastq-filter`)
* a tool to merge BAM or FASTQ reads using marked duplicates respectively unique molecular identifiers (UMIs) (`rbt collapse-reads-to-fragments bam|fastq`)
* a tool to generate interactive HTML based reports that offer multiple plots visualizing the provided genomics data in VCF and BAM format (`rbt vcf-report`)
* a tool for splitting VCF/BCF files into N equal chunks, including BND support (`rbt vcf-split`)

Further functionality is added as it is needed by the authors. Check out the [Contributing](#Contributing) section if you want contribute anything yourself.
For a list of changes, take a look at the [CHANGELOG](CHANGELOG.md).


## Installation

### Bioconda

Rust-Bio-Tools is available via [Bioconda](https://bioconda.github.io).
With Bioconda set up, installation is as easy as

    conda install rust-bio-tools

### Cargo

If the [Rust](https://www.rust-lang.org/tools/install) compiler and associated [Cargo](https://github.com/rust-lang/cargo/) are installed, Rust-Bio-Tools may be installed via

    cargo install rust-bio-tools

### Source

Download the source code and within the root directory of source run

    cargo install

## Usage and Documentation

Rust-Bio-Tools installs a command line utility `rbt`. Issue

    rbt --help

for a summary of all options and tools.

## Contributing

Any contributions are highly welcome. If you plan to contribute we suggest installing pre-commit hooks. To do so:
1. Install `pre-commit` as explained [here](https://pre-commit.com/#installation)
2. Run `pre-commit install` in the rust-bio-tools base directory

This should format, check and lint your code when committing.

## Authors

* Johannes Köster (https://koesterlab.github.io)
* Felix Mölder
* Henning Timm
