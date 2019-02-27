use bio::io::fastq;
use bio::stats::probs::{LogProb, PHREDProb};
use itertools::Itertools;
use ordered_float::NotNaN;
use std::cmp;

const PROB_CONFUSION: LogProb = LogProb(-1.0986122886681098); // (1 / 3).ln()

const ALLELES: &'static [u8] = b"ACGT";

/// Compute a consensus sequence for a collection of FASTQ reads.
///
/// For each position, compute the likelihood of each allele and
/// choose the most likely one. Write the most likely allele i.e. base
/// as sequence into the consensus sequence. The quality value is the
/// likelihood for this allele, encoded in PHRED+33.
pub fn calc_consensus(recs: &[fastq::Record], seqids: &[usize], uuid: &str) -> fastq::Record {
    let seq_len = recs[0].seq().len();
    let mut consensus_seq = Vec::with_capacity(seq_len);
    let mut consensus_qual = Vec::with_capacity(seq_len);

    // assert that all reads have the same length here
    let identical_lengths = || {
        let reference_length = recs[0].seq().len();
        recs.iter()
            .map(|rec| rec.seq().len())
            .all(|len| len == reference_length)
    };
    assert_eq!(
        identical_lengths(),
        true,
        "Read length of FASTQ records {:?} differ. Cannot compute consensus sequence.",
        seqids
    );
    // Potential workflow for different read lengths
    // compute consensus of all reads with max len
    // compute offset of all shorter reads
    // pad shorter reads
    // drop first consensus, compute consensus of full length reads and padded reads
    // ignore padded bases for consensus computation

    for i in 0..seq_len {
        #[allow(unused_doc_comments)]
        /// Compute the likelihood for the given allele and read position.
        /// The allele (A, C, G, or T) is an explicit parameter,
        /// the position i is captured by the closure.
        ///
        /// Likelihoods are managed in log space.
        /// A matching base is scored with (1 - PHRED score), a mismatch
        /// with PHRED score + confusion constant.
        let likelihood = |allele: &u8| {
            let mut lh = LogProb::ln_one(); // posterior: log(P(theta)) = 1
            for rec in recs {
                let q = LogProb::from(PHREDProb::from((rec.qual()[i] - 33) as f64));
                lh += if *allele == rec.seq()[i].to_ascii_uppercase() {
                    q.ln_one_minus_exp()
                } else {
                    q + PROB_CONFUSION
                };
            }
            lh
        };

        // Maximum a-posteriori estimate for the consensus base.
        // Find the allele (theta \in ACGT) with the highest likelihood
        // given the bases at this position, weighted with their quality values
        let likelihoods = ALLELES.iter().map(&likelihood).collect_vec();
        let max_posterior = likelihoods
            .iter()
            .enumerate()
            .max_by_key(|&(_, &lh)| NotNaN::new(*lh).unwrap()) // argmax of MAP
            .unwrap()
            .0;

        //
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
    let name = format!(
        "{}_consensus-read-from:{}",
        uuid,
        seqids.iter().map(|i| format!("{}", i)).join(",")
    );
    fastq::Record::with_attrs(&name, None, &consensus_seq, &consensus_qual)
}
