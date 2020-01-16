use crate::common::CalcConsensus;
use bio::io::fastq;
use bio::stats::probs::LogProb;
use derive_new::new;
use itertools::Itertools;
use regex;

const ALLELES: &'static [u8] = b"ACGT";

/// Compute a consensus sequence for a collection of FASTQ reads.
///
/// For each position, compute the likelihood of each allele and
/// choose the most likely one. Write the most likely allele i.e. base
/// as sequence into the consensus sequence. The quality value is the
/// likelihood for this allele, encoded in PHRED+33.
/// //TODO Generalize as this is identical to BAM except Offset and cigar/writing to record
#[derive(new)]
pub struct CalcNonOverlappingConsensus<'a> {
    recs: &'a [fastq::Record],
    seqids: &'a [usize],
    uuid: &'a str,
    verbose_read_names: bool,
}
impl<'a> CalcNonOverlappingConsensus<'a> {
    pub fn calc_consensus(&self) -> (fastq::Record, LogProb) {
        let seq_len = self.recs()[0].seq().len();
        let mut consensus_seq: Vec<u8> = Vec::with_capacity(seq_len);
        let mut consensus_qual: Vec<u8> = Vec::with_capacity(seq_len);

        // assert that all reads have the same length here
        assert_eq!(
            Self::validate_read_lengths(self.recs()),
            true,
            "Read length of FASTQ records {:?} differ. Cannot compute consensus sequence.",
            self.seqids()
        );

        // Potential workflow for different read lengths
        // compute consensus of all reads with max len
        // compute offset of all shorter reads
        // pad shorter reads
        // drop first consensus, compute consensus of full length reads and padded reads
        // ignore padded bases for consensus computation

        let mut consensus_lh = LogProb::ln_one();

        for i in 0..seq_len {
            // Maximum a-posteriori estimate for the consensus base.
            // Find the allele (theta \in ACGT) with the highest likelihood
            // given the bases at this position, weighted with their quality values
            let likelihoods = ALLELES
                .iter()
                .map(|a| Self::overall_allele_likelihood(&self, a, i))
                .collect_vec(); //Check this. See below
            Self::build_consensus_sequence(
                likelihoods,
                &mut consensus_lh,
                &mut consensus_seq,
                &mut consensus_qual,
                33.0,
            );
        }

        let name = match self.verbose_read_names {
            true => {

                // read_from:'Individual 73' at_locus:'0'
                let pattern = regex::Regex::new(r"at_locus:'(\d+)'").unwrap();
                let mut collected_stuff = Vec::new();
                for name in self.recs.iter().map(|x| x.desc().unwrap()) {
                    // dbg!(name);
                    // eprintln!("Found a read name: {}", name);
                    match pattern.captures(name) {
                        Some(m) => {
                            collected_stuff.push(
                                format!(
                                    "Locus_{}",
                                    m.get(1).map_or("", |m| m.as_str()),
                                    )
                            );
                        },
                        None => panic!("invalid read name"),
                    }
                }
                
                format!(
                    "{}_consensus-read-from:{}",
                    self.uuid(),
                    // self.seqids().iter().map(|i| format!("{}", i)).join(",")
                    collected_stuff.iter().map(|i| format!("{}", i)).join(",")
                )
            },
            false => format!(
                "{}_consensus-read-from:{}_reads",
                self.uuid(),
                self.seqids().len(),
            ),
        };

        (
            fastq::Record::with_attrs(&name, None, &consensus_seq, &consensus_qual),
            consensus_lh,
        )
    }

    pub fn recs(&self) -> &[fastq::Record] {
        self.recs
    }
}

//TODO Generalized as it is identical to BAM except Offset
impl<'a> CalcConsensus<'a, fastq::Record> for CalcNonOverlappingConsensus<'a> {
    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb {
        let mut lh = LogProb::ln_one(); // posterior: log(P(theta)) = 1
        for rec in self.recs() {
            lh += Self::allele_likelihood_in_rec(allele, rec.seq(), rec.qual(), i, 33);
        }
        lh
    }

    fn seqids(&self) -> &'a [usize] {
        self.seqids
    }

    fn uuid(&self) -> &'a str {
        self.uuid
    }
}

/// Compute a consensus sequence for a collection of paired-end FASTQ
/// reads taking overlap into account.
///
/// For each position, compute the likelihood of each allele and
/// choose the most likely one. Write the most likely allele i.e. base
/// as sequence into the consensus sequence. The quality value is the
/// likelihood for this allele, encoded in PHRED+33.
#[derive(new)]
pub struct CalcOverlappingConsensus<'a> {
    recs1: &'a [fastq::Record],
    recs2: &'a [fastq::Record],
    overlap: usize,
    seqids: &'a [usize],
    uuid: &'a str,
    verbose_read_names: bool,
}

//TODO Generalize as this is identical to BAM except Offset and cigar/writing to record
impl<'a> CalcOverlappingConsensus<'a> {
    pub fn calc_consensus(&self) -> (fastq::Record, LogProb) {
        let seq_len = self.recs1()[0].seq().len() + self.recs2()[0].seq().len() - self.overlap();
        let mut consensus_seq: Vec<u8> = Vec::with_capacity(seq_len);
        let mut consensus_qual: Vec<u8> = Vec::with_capacity(seq_len);

        // assert that all reads have the same length here
        assert_eq!(
            Self::validate_read_lengths(self.recs1()),
            true,
            "Read length of FASTQ forward records {:?} differ. Cannot compute consensus sequence.",
            self.seqids()
        );

        assert_eq!(
            Self::validate_read_lengths(self.recs2()),
            true,
            "Read length of FASTQ reverse records {:?} differ. Cannot compute consensus sequence.",
            self.seqids()
        );
        let mut consensus_lh = LogProb::ln_one();

        for i in 0..seq_len {
            let likelihoods = ALLELES
                .iter()
                .map(|a| Self::overall_allele_likelihood(&self, a, i))
                .collect_vec(); //This will be calculated every iteration
            Self::build_consensus_sequence(
                likelihoods,
                &mut consensus_lh,
                &mut consensus_seq,
                &mut consensus_qual,
                33.0,
            );
        }
        let name = match self.verbose_read_names {
            true => format!(
                "{}_consensus-read-from:{}",
                self.uuid(),
                self.seqids().iter().map(|i| format!("{}", i)).join(",")
            ),
            false => format!(
                "{}_consensus-read-from:{}_reads",
                self.uuid(),
                self.seqids().len(),
            ),
        };
        (
            fastq::Record::with_attrs(&name, None, &consensus_seq, &consensus_qual),
            consensus_lh,
        )
    }

    fn recs1(&self) -> &[fastq::Record] {
        self.recs1
    }

    fn recs2(&self) -> &[fastq::Record] {
        self.recs2
    }

    fn overlap(&self) -> usize {
        self.overlap
    }
}

impl<'a> CalcConsensus<'a, fastq::Record> for CalcOverlappingConsensus<'a> {
    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb {
        let mut lh = LogProb::ln_one();
        for (rec1, rec2) in self.recs1().into_iter().zip(self.recs2()) {
            if i < rec1.seq().len() {
                lh += Self::allele_likelihood_in_rec(allele, rec1.seq(), rec1.qual(), i, 33);
            };
            if i >= rec1.seq().len() - self.overlap() {
                let rec2_i = i - (rec1.seq().len() - self.overlap());
                let rec2_seq = bio::alphabets::dna::revcomp(rec2.seq());
                let rec2_qual: Vec<u8> = rec2.qual().iter().rev().cloned().collect();
                lh += Self::allele_likelihood_in_rec(allele, &rec2_seq, &rec2_qual, rec2_i, 33);
            };
        }
        lh
    }

    fn seqids(&self) -> &'a [usize] {
        self.seqids
    }

    fn uuid(&self) -> &'a str {
        self.uuid
    }
}
