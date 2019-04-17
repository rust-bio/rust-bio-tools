use bio::stats::probs::{LogProb, PHREDProb};
use itertools::Itertools;
use ordered_float::NotNaN;
use rust_htslib::bam;
use std::cmp;

const PROB_CONFUSION: LogProb = LogProb(-1.0986122886681098); // (1 / 3).ln()

const ALLELES: &'static [u8] = b"ACGT";

pub struct CalcOverlappingConsensus<'a> {
    recs1: &'a [bam::Record],
    recs2: &'a [bam::Record],
    overlap: usize,
    uuid: &'a str,
}

impl<'a> CalcOverlappingConsensus<'a> {
    pub fn new(
        recs1: &'a [bam::Record],
        recs2: &'a [bam::Record],
        overlap: usize,
        uuid: &'a str,
    ) -> Self {
        CalcOverlappingConsensus {
            recs1,
            recs2,
            overlap,
            uuid,
        }
    }
    pub fn calc_consensus(&self) -> (bam::Record, LogProb) {
        let seq_len = self.recs1()[0].insert_size() as usize;
        let mut consensus_seq: Vec<u8> = Vec::with_capacity(seq_len);
        let mut consensus_qual: Vec<u8> = Vec::with_capacity(seq_len);

        // assert that all reads have the same length here
        assert_eq!(
            Self::validate_read_lengths(self.recs1()),
            true,
            "Read length of BAM forward records differ. Cannot compute consensus sequence.",
        );

        assert_eq!(
            Self::validate_read_lengths(self.recs2()),
            true,
            "Read length of BAM reverse records differ. Cannot compute consensus sequence."
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
        //TODO Add seq ids
        let name = format!("{}_consensus-read", self.uuid());
        let mut cigar_string = seq_len.to_string();
        cigar_string.push('M');
        let cigar = bam::record::CigarString::from_str(cigar_string.as_str()).unwrap();
        let mut consensus_rec = bam::Record::new();
        consensus_rec.set(name.as_bytes(), &cigar, &consensus_seq, &consensus_qual);
        (consensus_rec, consensus_lh)
    }
    fn recs1(&self) -> &[bam::Record] {
        self.recs1
    }
    fn recs2(&self) -> &[bam::Record] {
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
                lh +=
                    Self::allele_likelihood_in_rec(allele, &rec1.seq().as_bytes(), rec1.qual(), i);
            };
            if i >= rec1.seq().len() - self.overlap() {
                let rec2_i = i - (rec1.seq().len() - self.overlap());
                lh += Self::allele_likelihood_in_rec(
                    allele,
                    &rec2.seq().as_bytes(),
                    &rec2.qual(),
                    rec2_i,
                );
            };
        }
        lh
    }
    fn uuid(&self) -> &'a str {
        self.uuid
    }
}

pub struct CalcNonOverlappingConsensus<'a> {
    recs: &'a [bam::Record],
    uuid: &'a str,
}

impl<'a> CalcNonOverlappingConsensus<'a> {
    pub fn new(recs: &'a [bam::Record], uuid: &'a str) -> Self {
        CalcNonOverlappingConsensus { recs, uuid }
    }
    pub fn calc_consensus(&self) -> (bam::Record, LogProb) {
        let seq_len = self.recs()[0].seq().len();
        let mut consensus_seq: Vec<u8> = Vec::with_capacity(seq_len);
        let mut consensus_qual: Vec<u8> = Vec::with_capacity(seq_len);

        // assert that all reads have the same length here
        assert_eq!(
            Self::validate_read_lengths(self.recs()),
            true,
            "Read length of BAM records differ. Cannot compute consensus sequence.",
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
        //TODO Add seq ids
        let name = format!("{}_consensus-read", self.uuid());
        let mut cigar_string = seq_len.to_string();
        cigar_string.push('M');
        let cigar = bam::record::CigarString::from_str(cigar_string.as_str()).unwrap();
        let mut consensus_rec = bam::Record::new();
        consensus_rec.set(name.as_bytes(), &cigar, &consensus_seq, &consensus_qual);
        (consensus_rec, consensus_lh)
    }
    pub fn recs(&self) -> &[bam::Record] {
        self.recs
    }
}

impl<'a> CalcConsensus<'a> for CalcNonOverlappingConsensus<'a> {
    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb {
        let mut lh = LogProb::ln_one(); // posterior: log(P(theta)) = 1
        for rec in self.recs() {
            lh += Self::allele_likelihood_in_rec(allele, &rec.seq().as_bytes(), rec.qual(), i);
        }
        lh
    }
    fn uuid(&self) -> &'a str {
        self.uuid
    }
}

pub trait CalcConsensus<'a> {
    fn validate_read_lengths(recs: &[bam::Record]) -> bool {
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
        let q = LogProb::from(PHREDProb::from((qual[i]) as f64));
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
        consensus_qual.push(cmp::min(41, (truncated_quality) as u64) as u8);
    }
    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb;
    fn uuid(&self) -> &'a str;
}
