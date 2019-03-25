use bio::io::fastq;
use bio::stats::probs::{LogProb, PHREDProb};
use itertools::Itertools;
use ordered_float::NotNaN;
use std::cmp;

const PROB_CONFUSION: LogProb = LogProb(-1.0986122886681098); // (1 / 3).ln()

const ALLELES: &'static [u8] = b"ACGT";

pub trait CalcConsensus<'a> {
    fn validate_read_lengths(recs: &[fastq::Record]) -> bool {
        let reference_length = recs[0].seq().len();
        recs.iter()
            .map(|rec| rec.seq().len())
            .all(|len| len == reference_length)
    }
    #[allow(unused_doc_comments)]
    /// Compute the likelihood for the given allele and read position.
    /// The allele (A, C, G, or T) is an explicit parameter,
    /// the position i is captured by the closure.
    ///
    /// Likelihoods are managed in log space.
    /// A matching base is scored with (1 - PHRED score), a mismatch
    /// with PHRED score + confusion constant.
    fn allele_likelihood_in_rec(allele: &u8, seq: &[u8], qual: &[u8], i: usize) -> LogProb {
        let q = LogProb::from(PHREDProb::from((qual[i] - 33) as f64));
        if *allele == seq[i].to_ascii_uppercase() {
            q.ln_one_minus_exp()
        } else {
            q + PROB_CONFUSION
        }
    }
    fn build_consensus_sequence(
        likelihoods: Vec<LogProb>,
        consensus_lh: &mut LogProb,
        consensus_seq: &mut Vec<u8>,
        consensus_qual: &mut Vec<u8>,
    ) {
        let (max_posterior, allele_lh) = likelihoods
            .iter()
            .enumerate()
            .max_by_key(|&(_, &lh)| NotNaN::new(*lh).unwrap())
            .unwrap();
        *consensus_lh += *allele_lh;
        let marginal = LogProb::ln_sum_exp(&likelihoods);
        // new base: MAP
        consensus_seq.push(ALLELES[max_posterior]);
        // new qual: (1 - MAP)
        let qual = (likelihoods[max_posterior] - marginal).ln_one_minus_exp();
        // Assume the maximal quality, if the likelihood is infinite
        let truncated_quality: f64;
        if (*PHREDProb::from(qual)).is_infinite() {
            truncated_quality = 41.0;
        } else {
            truncated_quality = *PHREDProb::from(qual);
        }
        // Truncate quality values to PHRED+33 range
        consensus_qual.push(cmp::min(74, (truncated_quality + 33.0) as u64) as u8);
    }

    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb;
    fn seqids(&self) -> &'a [usize];
    fn uuid(&self) -> &'a str;
}

/// Compute a consensus sequence for a collection of FASTQ reads.
///
/// For each position, compute the likelihood of each allele and
/// choose the most likely one. Write the most likely allele i.e. base
/// as sequence into the consensus sequence. The quality value is the
/// likelihood for this allele, encoded in PHRED+33.
pub struct CalcNonOverlappingConsensus<'a> {
    recs: &'a [fastq::Record],
    seqids: &'a [usize],
    uuid: &'a str,
}

impl<'a> CalcNonOverlappingConsensus<'a> {
    pub fn new(recs: &'a [fastq::Record], seqids: &'a [usize], uuid: &'a str) -> Self {
        CalcNonOverlappingConsensus { recs, seqids, uuid }
    }
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
            );
        }
        let name = format!(
            "{}_consensus-read-from:{}",
            self.uuid(),
            self.seqids().iter().map(|i| format!("{}", i)).join(",")
        );
        (
            fastq::Record::with_attrs(&name, None, &consensus_seq, &consensus_qual),
            consensus_lh,
        )
    }
    pub fn recs(&self) -> &[fastq::Record] {
        self.recs
    }
}

impl<'a> CalcConsensus<'a> for CalcNonOverlappingConsensus<'a> {
    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb {
        let mut lh = LogProb::ln_one(); // posterior: log(P(theta)) = 1
        for rec in self.recs() {
            lh += Self::allele_likelihood_in_rec(allele, rec.seq(), rec.qual(), i);
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
pub struct CalcOverlappingConsensus<'a> {
    recs1: &'a [fastq::Record],
    recs2: &'a [fastq::Record],
    overlap: usize,
    seqids: &'a [usize],
    uuid: &'a str,
}

impl<'a> CalcOverlappingConsensus<'a> {
    pub fn new(
        recs1: &'a [fastq::Record],
        recs2: &'a [fastq::Record],
        overlap: usize,
        seqids: &'a [usize],
        uuid: &'a str,
    ) -> Self {
        CalcOverlappingConsensus {
            recs1,
            recs2,
            overlap,
            seqids,
            uuid,
        }
    }
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
            );
        }
        let name = format!(
            "{}_consensus-read-from:{}",
            self.uuid(),
            self.seqids().iter().map(|i| format!("{}", i)).join(",")
        );
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

impl<'a> CalcConsensus<'a> for CalcOverlappingConsensus<'a> {
    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb {
        let mut lh = LogProb::ln_one();
        for (rec1, rec2) in self.recs1().into_iter().zip(self.recs2()) {
            if i < rec1.seq().len() {
                lh += Self::allele_likelihood_in_rec(allele, rec1.seq(), rec1.qual(), i);
            };
            if i >= rec1.seq().len() - self.overlap() {
                let rec2_i = i - (rec1.seq().len() - self.overlap());
                let rec2_seq = bio::alphabets::dna::revcomp(rec2.seq());
                let rec2_qual: Vec<u8> = rec2.qual().iter().rev().cloned().collect();
                lh += Self::allele_likelihood_in_rec(allele, &rec2_seq, &rec2_qual, rec2_i);
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
