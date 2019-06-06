use bio::stats::probs::{LogProb, PHREDProb};
use bio_types::sequence::SequenceRead;
use ordered_float::NotNaN;
use std::cmp;

const PROB_CONFUSION: LogProb = LogProb(-1.0986122886681098); // (1 / 3).ln()
const ALLELES: &'static [u8] = b"ACGT";

pub trait CalcConsensus<'a, R: SequenceRead> {
    fn validate_read_lengths(recs: &[R]) -> bool {
        let reference_length = recs[0].len();
        recs.iter()
            .map(|rec| rec.len())
            .all(|len| len == reference_length)
    }
    /// Compute the likelihood for the given allele and read position.
    /// The allele (A, C, G, or T) is an explicit parameter,
    /// the position i is captured by the closure.
    ///
    /// Likelihoods are managed in log space.
    /// A matching base is scored with (1 - PHRED score), a mismatch
    /// with PHRED score + confusion constant.
    fn allele_likelihood_in_rec(
        allele: &u8,
        seq: &[u8],
        qual: &[u8],
        i: usize,
        offset: u8,
    ) -> LogProb {
        let q = LogProb::from(PHREDProb::from((qual[i] - offset) as f64));
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
        offset: f64,
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
        consensus_qual
            .push(cmp::min(41 + offset as u64, (truncated_quality + offset) as u64) as u8);
    }

    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb;
    fn seqids(&self) -> &'a [usize];
    fn uuid(&self) -> &'a str;
}
