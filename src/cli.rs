use std::path::PathBuf;
use structopt::StructOpt;

#[derive(StructOpt)]
#[structopt(
    about = "A set of ultra-fast command line utilities for bioinformatics tasks based on Rust-Bio.",
    author = "Johannes Köster <johannes.koester@tu-dortmund.de>",
    name = "Rust-Bio-Tools"
)]
pub(crate) struct Rbt {
    #[structopt(long, short, help = "Verbose output.")]
    pub(crate) verbose: bool,

    #[structopt(subcommand)]
    pub(crate) cmd: Command,
}

#[derive(StructOpt)]
pub(crate) enum Command {
    /// Split FASTQ file from STDIN into N chunks.
    ///
    /// Example:
    /// rbt fastq-split A.fastq B.fastq < test.fastq
    FastqSplit {
        #[structopt(parse(from_os_str), help = "File name(s) for the chunks to create.")]
        chunks: Vec<PathBuf>,
    },
    /// Remove records from a FASTQ file (from STDIN), output to STDOUT.
    ///
    /// Example:
    /// rbt fastq-filter ids.txt < test.fastq > filtered.fastq
    FastqFilter {
        #[structopt(parse(from_os_str))]
        /// File with list of record IDs to remove, one per line.
        ids: PathBuf,
    },

    /// Print depth of BAM or CRAM file at given positions from STDIN (tab separated: chrom, pos).
    ///
    /// Usage:
    /// $ rbt bam-depth test.bam < pos.txt > depth.txt
    ///
    /// The positions file contains the name of one reference sequence and one position per line (tab separated).
    /// Example:
    ///
    /// 16	1
    /// 17	38
    /// 17	39
    ///
    /// Depths are written to stdout as tab-separated lines, similar to the positions input.
    /// Example:
    ///
    /// 16	1	0
    /// 17	38	14
    /// 17	39	13
    BamDepth {
        /// Path to indexed BAM file.
        #[structopt(parse(from_os_str))]
        bam_path: PathBuf,

        /// Maximum read length to consider. This affects the speed of the involved pileup.
        /// Reads longer than this length can be missed when calculating the depth.
        #[structopt(long, short, default_value = "1000")]
        max_read_length: u32,

        /// Skip reads with mask bits unset [].
        #[structopt(long = "incl-flags", short, default_value = "0")]
        include_flags: u16,

        /// Skip reads with mask bits set [UNMAP, SECONDARY, QCFAIL, DUP].
        #[structopt(long = "excl-flags", short, default_value = "1796")]
        exclude_flags: u16,

        /// Minimum mapping quality.
        #[structopt(long, short = "q", default_value = "0")]
        min_mapq: u8,
    },

    /// Convert any IUPAC codes in alleles into Ns (in order to comply with VCF 4 specs).
    /// Reads VCF/BCF from STDIN and writes BCF to STDOUT.
    ///
    /// Example:
    /// rbt vcf-fix-iupac-alleles < test.vcf > fixed.bcf
    VcfFixIupacAlleles {},

    /// Convert VCF/BCF file from STDIN to tab-separated TXT file at STDOUT.
    /// INFO and FORMAT tags have to be selected explicitly.
    ///
    /// Example:
    /// rbt vcf-to-txt --genotypes --fmt S --info T X SOMATIC < test.vcf > variant-table.txt
    ///
    /// The resulting table can be e.g. parsed with PANDAS in Python:
    ///
    /// pd.read_table("variants.txt", header=[0, 1])
    VcfToTxt {
        /// Select INFO tags
        #[structopt(long, short, value_name = "NAME")]
        info: Vec<String>,

        /// Select FORMAT tags.
        #[structopt(long = "fmt", short, value_name = "NAME")]
        format: Vec<String>,

        /// Display genotypes.
        #[structopt(long, short)]
        genotypes: bool,
    },

    /// Annotate for each variant in a VCF/BCF at STDIN whether it is contained in a
    /// given second VCF/BCF. The matching is fuzzy for indels and exact for SNVs.
    /// Results are printed as BCF to STDOUT, with an additional INFO tag MATCHING.
    /// The two vcfs do not have to be sorted.
    ///
    /// Example:
    /// rbt vcf-match dbsnp.vcf < calls.vcf | bcftools view
    VcfMatch {
        /// VCF/BCF file to match against.
        #[structopt(parse(from_os_str))]
        vcf: PathBuf,

        /// Maximum distance between centres of two indels considered to match.
        #[structopt(long, short = "d", value_name = "INT", default_value = "20")]
        max_dist: u32,

        /// Maximum difference between lengths of two indels.
        #[structopt(long, short = "l", value_name = "INT", default_value = "10")]
        max_len_diff: u32,
    },

    /// Annotate b-allele frequency for each single nucleotide variant and sample.
    ///
    /// Example:
    /// rbt vcf-baf < calls.bcf > annotated.bcf
    VcfBaf {},

    /// Looks for interacting drugs in DGIdb and annotates them for every gene in every record.
    ///
    /// Example:
    /// rbt vcf-annotate-dgidb input.vcf > output.vcf
    VcfAnnotateDgidb {
        /// VCF/BCF file to be extended by dgidb drug entries
        #[structopt()]
        vcf: String,

        /// Url prefix for requesting interaction drugs by gene names.
        #[structopt(
            long,
            short = "p",
            default_value = "http://dgidb.org/api/v2/interactions.json?genes="
        )]
        api_path: String,

        /// Info field name to be used for annotation.
        #[structopt(long, short = "f", default_value = "dgiDB_drugs")]
        field: String,

        /// A list of data sources included in query. If omitted all sources are considered.
        /// A list of all sources can be found at http://dgidb.org/api/v2/interaction_sources.json
        #[structopt(long, short = "s", value_name = "STR")]
        datasources: Option<Vec<String>>,

        /// Number of genes to submit per api request. A lower value increases the number of api requests in return.
        /// Too many requests could be rejected by the DGIdb server.
        #[structopt(long, short = "g", default_value = "500")]
        genes_per_request: usize,
    },
}
