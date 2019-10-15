# Rust-Bio-Tools

[![Travis](https://img.shields.io/travis/rust-bio/rust-bio-tools/master.svg?style=flat-square)](https://travis-ci.org/rust-bio/rust-bio-tools)

A set of ultra fast and robust command line utilities for bioinformatics tasks based on Rust-Bio.
Rust-Bio-Tools provides a command `rbt`, which currently supports the following operations:

* a linear time implementation for fuzzy matching of two vcf/bcf files (`rbt vcf-match`)
* a vcf/bcf to txt converter, that flexibly allows to select tags and properly handles multiallelic sites (`rbt vcf-to-txt`)
* a linear time round-robin FASTQ splitter that splits a given FASTQ files into a given number of chunks (`rbt fastq-split`)
* a linear time extraction of depth information from BAMs at given loci (`rbt bam-depth`)
* a utility to quickly filter records from a FASTQ file (`rbt fastq-filter`)
* a tool to merge BAM or FASTQ reads using marked duplicates respectively unique molecular identifiers (UMIs) (`rbt call-consensus-reads bam|fastq`)

Further functionality is added as it is needed by the authors. Any contributions are highly welcome.
For a list of changes, take a look at the [CHANGELOG](CHANGELOG.md).


## Installation

### Bioconda

Rust-Bio-Tools is available via [Bioconda](https://bioconda.github.io).
With Bioconda set up, installation is as easy as

    conda install rust-bio-tools

### Manual

Download the source code and issue

    cargo install

from the root directory of the source.

## Usage and Documentation

Rust-Bio-Tools installs a command line utility `rbt`. Issue

    rbt --help

for a summary of all options and tools.



## Authors

* Johannes Köster (https://koesterlab.github.io)
* Felix Mölder
* Henning Timm

