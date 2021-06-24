use crate::common::Region;
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
    #[structopt(author = "Johannes Köster <johannes.koester@tu-dortmund.de>")]
    FastqSplit {
        #[structopt(parse(from_os_str), help = "File name(s) for the chunks to create.")]
        chunks: Vec<PathBuf>,
    },
    /// Remove records from a FASTQ file (from STDIN), output to STDOUT.
    ///
    /// Example:
    /// rbt fastq-filter ids.txt < test.fastq > filtered.fastq
    #[structopt(author = "Erik Clarke <ecl@pennmedicine.upenn.edu>")]
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
    /// 16    1
    /// 17    38
    /// 17    39
    ///
    /// Depths are written to stdout as tab-separated lines, similar to the positions input.
    /// Example:
    ///
    /// 16    1    0
    /// 17    38    14
    /// 17    39    13
    #[structopt(author = "Johannes Köster <johannes.koester@tu-dortmund.de>")]
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
    #[structopt(author = "Johannes Köster <johannes.koester@tu-dortmund.de>")]
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
    #[structopt(author = "Johannes Köster <johannes.koester@tu-dortmund.de>")]
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
    #[structopt(author = "Johannes Köster <johannes.koester@tu-dortmund.de>")]
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
    #[structopt(
        author = "Johannes Köster <johannes.koester@uni-due.de>, Jan Forster <j.forster@dkfz.de>"
    )]
    VcfBaf {},

    /// Looks for interacting drugs in DGIdb and annotates them for every gene in every record.
    ///
    /// Example:
    /// rbt vcf-annotate-dgidb input.vcf > output.vcf
    #[structopt(author = "Felix Mölder <felix.moelder@uni-due.de>")]
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

    /// Creates report from a given csv file containing a table with the given data
    /// Examples:
    /// With current directory as default ouput path:
    /// rbt csv-report path/to/table.csv --rows-per-page 100 --sort-column "p-value" --sort-order ascending
    #[structopt(author = "Felix Wiegand <felix.wiegand@tu-dortmund.de>")]
    CsvReport {
        /// CSV file including the data for the report.
        #[structopt()]
        csv_path: String,

        /// Sets the numbers of rows of each table per page. Default is 100.
        #[structopt(long, short = "r", default_value = "100")]
        rows_per_page: u32,

        /// Column that the data should be sorted by.
        #[structopt(long, short = "c")]
        sort_column: Option<String>,

        /// Order the data ascending or descending. Default is descending.
        #[structopt(long, short = "o", default_value = "descending", possible_values = &["ascending","descending"])]
        sort_order: String,

        /// Change the separator of the csv file to tab or anything else. Default is ",".
        #[structopt(long, short = "s", default_value = ",")]
        separator: char,

        /// Configure a custom formatter function for each column by providing a file containing a javascript object with csv column title as the key and a format function as the value.
        /// More information on the formatting functions and how to use them here: https://bootstrap-table.com/docs/api/column-options/#formatter.
        #[structopt(long, short = "f")]
        formatter: Option<String>,

        /// Pins the table until the given column such that scrolling to the right does not hide the given column and those before.
        #[structopt(long, short = "p")]
        pin_until: Option<String>,

        /// Relative output path for the report files. Default value is the current directory.
        #[structopt(default_value = ".")]
        output_path: String,
    },

    /// Creates a html file with a vega visualization of the given bam region
    /// Example:
    /// rbt plot-bam -b input.bam -g 2:132424-132924 -r input.fa > plot.html
    #[structopt(author = "Felix Wiegand <felix.wiegand@tu-dortmund.de>")]
    PlotBam {
        /// BAM file to be visualized.
        #[structopt(long, short = "b", parse(from_os_str))]
        bam_path: Vec<PathBuf>,

        /// Path to the reference fasta file.
        #[structopt(long, short = "r", parse(from_os_str))]
        reference: PathBuf,

        /// Chromosome and region for the visualization. Example: 2:132424-132924
        #[structopt(long, short = "g")]
        region: Region,

        /// Set the maximum rows that will be shown in the alignment plots.
        #[structopt(long, short = "d", default_value = "500")]
        max_read_depth: u32,
    },

    /// Creates report from a given VCF file including a visual plot
    /// for every variant with the given BAM and FASTA file.
    /// The VCF file has to be annotated with VEP, using the options --hgvs and --hgvsg.
    ///
    /// Examples:
    /// With current directory as default ouput path:
    /// rbt vcf-report fasta.fa --vcfs a=a.vcf b=b.vcf --bams a:sample1=a.bam b:sample1=b.bam
    /// With custom directory as default ouput path:
    /// rbt vcf-report fasta.fa --vcfs a=a.vcf b=b.vcf --bams a:sample1=a.bam b:sample1=b.bam -- my/output/path/
    /// With custom info tags in table report:
    /// rbt vcf-report fasta.fa --vcfs a=a.vcf b=b.vcf --bams a:sample1=a.bam b:sample1=b.bam --info PROB_SOMATIC PROB_GERMLINE
    #[structopt(
        author = "Johannes Köster <johannes.koester@uni-due.de>, Felix Wiegand <felix.wiegand@tu-dortmund.de>"
    )]
    VcfReport {
        /// FASTA file containing the reference genome for the visual plot
        #[structopt()]
        fasta: String,

        /// VCF files to include (multi-sample). Group is the name that will be used in the oncoprint. There needs to be one corresponding BAM file for each sample of a VCF/BCF file. Please only use VCF/BCF files annotated by VEP.
        #[structopt(long, short = "v", value_name = "GROUP=VCF_FILE")]
        vcfs: Vec<String>,

        /// VCF files to include (multi-sample). Group is the name that will be used in the oncoprint. There needs to be one corresponding BAM file for each sample of a VCF/BCF file. Please only use VCF/BCF files annotated by VEP.
        #[structopt(long, short = "b", value_name = "GROUP:SAMPLE=BAM_FILE")]
        bams: Vec<String>,

        /// Set the maximum number of cells in the oncoprint per page. Lowering max-cells should improve the performance of the plots in the browser. Default value is 1000.
        #[structopt(long, short = "c", default_value = "1000")]
        cells: u32,

        /// Set the maximum lines of reads that will be shown in the alignment plots. Default value is 500.
        #[structopt(long, short = "d", default_value = "500")]
        max_read_depth: u32,

        /// Add custom values from the info field to each variant as a data attribute to access them via the custom javascript. Multiple fields starting with the same prefix can be added by placing '*' at the end of a prefix.
        #[structopt(long, short = "i", value_name = "INFO_TAG")]
        infos: Option<Vec<String>>,

        /// Add custom values from the format field to each variant as a data attribute to access them via the custom javascript. All given format values will also be inserted into the main table.
        #[structopt(long, short = "f", value_name = "FORMAT_TAG")]
        formats: Option<Vec<String>>,

        /// Add multiple keys from the info field of your vcf to the plots of the first and second stage of the report.
        #[structopt(long, value_name = "PLOT_INFO")]
        plot_info: Option<Vec<String>>,

        /// Change the default javascript file for the table-report to a custom one to add own plots or tables to the sidebar by appending these to an empty div in the HTML template.
        #[structopt(long, short = "j", value_name = "JS_FILE_PATH")]
        custom_js_template: Option<String>,

        /// Add one or multiple js file (e.g. libraries) for usage in the custom-js-file. The ordering of the arguments will be the same as they will be imported.
        #[structopt(long, short = "l", value_name = "JS_FILE_PATH")]
        custom_js_files: Option<Vec<String>>,

        /// Add a TSV file that contains one or multiple custom values for each sample for the oncoprint. First column has to be the sample name, followed by one or more columns with custom values. Make sure you include one row for each given sample.
        #[structopt(long, short = "t", value_name = "TSV_FILE_PATH")]
        tsv: Option<String>,

        /// Sets the number of threads used to build the table reports.
        #[structopt(long, default_value = "0")]
        threads: usize,

        /// Relative output path for the report files. Default value is the current directory.
        #[structopt(default_value = ".")]
        output_path: String,
    },

    /// Split a given VCF/BCF file into N chunks of approximately the same size. Breakends are kept together.
    /// Output type is always BCF.
    ///
    /// Example:
    /// rbt vcf-split input.bcf output1.bcf output2.bcf output3.bcf ... outputN.bcf
    #[structopt(author = "Johannes Köster <johannes.koester@uni-due.de>")]
    VcfSplit {
        #[structopt(parse(from_os_str), help = "Input VCF/BCF that shall be splitted.")]
        input: PathBuf,

        #[structopt(
            parse(from_os_str),
            help = "BCF files to split into. Breakends are kept together. Each file will contain approximately the same number of records."
        )]
        output: Vec<PathBuf>,
    },

    /// Tool to predict maximum likelihood fragment sequence from FASTQ or BAM files.
    ///
    /// Requirements:
    /// - starcode
    #[structopt(
        author = "Johannes Köster <johannes.koester@uni-due.de>, Henning Timm <henning.timm@tu-dortmund.de>, Felix Mölder <felix.moelder@uni-due.de>"
    )]
    CollapseReadsToFragments {
        #[structopt(subcommand)]
        cmd: CollapseReadsToFragmentsSubcommand,
    },

    /// Tool to compute stats on sequence file (from STDIN), output is in YAML with fields:
    /// - min: length of shortest sequence
    /// - max: length of longest sequence
    /// - average: average length of sequence
    /// - median: median length of sequence
    /// - nb_reads: number of reads
    /// - nb_bases: number of bases
    /// - n50: N50 of sequences
    ///
    /// Example:
    /// rbt sequence-stats < test.fasta
    /// rbt sequence-stats -q < test.fastq
    #[structopt(author = "Pierre Marijon <pmarijon@mpi-inf.mpg.de>")]
    SequenceStats {
        #[structopt(
            long,
            short = "q",
            help = "Flag to indicate the sequence in stdin is in fastq format."
        )]
        fastq: bool,
    },
}

#[derive(StructOpt)]
pub enum CollapseReadsToFragmentsSubcommand {
    /// Tool to merge sets of reads from paired FASTQ files that share the UMI and have similar read sequence. The result is a maximum likelihood fragment sequence per set with base quality scores improved accordingly.
    ///
    /// Takes two FASTQ files (forward and reverse) and returns two FASTQ files in which all PCR duplicates have been merged into a consensus read.
    /// Duplicates are identified by a Unique Molecular Identifier (UMI).
    ///
    /// Assumptions:
    /// - Reads are of equal length
    /// - UMI is the prefix of the reads
    ///
    /// Example:
    /// rbt collapse-reads-to-fragments fastq \
    /// reads_1.fq reads_2.fq \    # input files
    /// merged_1.fq merged_2.fq \  # output files
    /// -l 13 \                    # length of UMI
    /// -d 1 \                     # max hamming distance of UMIs within a cluster
    /// -D 2 \                     # max hamming distance of sequences within a cluster
    /// --umi-on-reverse           # UMI is the prefix of the reverse read
    #[structopt(
        author = "Johannes Köster <johannes.koester@uni-due.de>, Henning Timm <henning.timm@tu-dortmund.de>, Felix Mölder <felix.moelder@uni-due.de>"
    )]
    Fastq {
        #[structopt(parse(from_os_str), help = "Input FASTQ file with forward reads.")]
        fq1: PathBuf,

        #[structopt(parse(from_os_str), help = "Input FASTQ file with reverse reads.")]
        fq2: PathBuf,

        #[structopt(parse(from_os_str), help = "Output FASTQ file with forward reads")]
        consensus_fq1: PathBuf,

        #[structopt(parse(from_os_str), help = "Output FASTQ file with reverse reads")]
        consensus_fq2: PathBuf,

        #[structopt(
            parse(from_os_str),
            requires_all(&["insert-size", "std-dev"]),
            help = "Output FASTQ file for overlapping consensus reads  (Required for calculating overlapping consensus only)"
        )]
        consensus_fq3: Option<PathBuf>,

        #[structopt(
            long,
            short = "d",
            default_value = "1",
            help = "Maximum hamming distance between the UMIs of any pair of reads in the same cluster."
        )]
        max_umi_dist: usize,

        #[structopt(
            long,
            short = "l",
            default_value = "8",
            help = "Length of UMI in read."
        )]
        umi_len: usize,

        #[structopt(long, short = "D", possible_values = &["1","2","3","4","5","6","7","8"], default_value = "2", help = "Maximum hamming distance between the sequences of any pair of reads in the same cluster.")]
        max_seq_dist: usize,

        #[structopt(long, short = "u", help = "Set if UMI is on reverse read")]
        umi_on_reverse: bool,

        #[structopt(
            long,
            help = "Add list of reads that were merged for each consensus read. Note that this can yield very long FASTQ name lines which cannot be handled by some tools."
        )]
        verbose_read_names: bool,

        #[structopt(
            long,
            short = "i",
            requires = "consensus-fq3",
            help = "Expected insert size of sequenced fragment (Required for calculating overlapping consensus only)"
        )]
        insert_size: Option<usize>,

        #[structopt(
            long,
            short = "s",
            requires = "consensus-fq3",
            help = "Standard deviation of expected insert size. Defines search space of the most likely overlap. (Required for calculating overlapping consensus only)"
        )]
        std_dev: Option<usize>,
    },

    /// Tool to merge sets of PCR duplicate reads from a BAM file into one maximum likelihood fragment sequence each with accordingly improved base quality scores.
    ///
    /// Takes a BAM file and returns a BAM file in which all PCR duplicates have been merged into a consensus read.
    /// Duplicates must be marked by Picard Tools using the TAG_DUPLICATE_SET_MEMBERS option.
    ///
    /// Assumptions:
    /// - Reads are of equal length
    /// - Reads are marked by Picard Tools
    #[structopt(author = "Felix Mölder <felix.moelder@uni-due.de>")]
    Bam {
        #[structopt(parse(from_os_str), help = "Input BAM file with marked duplicates")]
        bam: PathBuf,

        #[structopt(parse(from_os_str), help = "Output FASTQ file with forward reads")]
        consensus_fq1: PathBuf,

        #[structopt(parse(from_os_str), help = "Output FASTQ file with reverse reads")]
        consensus_fq2: PathBuf,

        #[structopt(
            parse(from_os_str),
            help = "Output FASTQ file for overlapping consensus reads."
        )]
        consensus_fq_se: PathBuf,

        #[structopt(
            parse(from_os_str),
            help = "Output FASTQ file for overlapping consensus reads."
        )]
        skipped_bam: PathBuf,

        #[structopt(
            long,
            help = "Add list of reads that were merged for each consensus read. Note that this can yield very long FASTQ name lines which cannot be handled by some tools."
        )]
        verbose_read_names: bool,
    },
}
