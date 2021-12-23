use anyhow::Context;
use approx::relative_eq;
use bio::stats::probs::{LogProb, PHREDProb};
use bio_types::sequence::SequenceRead;
use itertools::Itertools;
use ordered_float::NotNaN;
use std::cmp;
use std::collections::HashMap;
use std::str::FromStr;

const PROB_CONFUSION: LogProb = LogProb(-1.0986122886681098); // (1 / 3).ln()
const ALLELES: &[u8] = b"ACGT";

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
        if relative_eq!(*likelihoods[0], *likelihoods[1])
            && relative_eq!(*likelihoods[1], *likelihoods[2])
            && relative_eq!(*likelihoods[2], *likelihoods[3])
        {
            consensus_seq.push(b'N');
            consensus_qual.push(offset as u8);
        } else {
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
            let truncated_quality: f64 = if (*PHREDProb::from(qual)).is_infinite() {
                93.0
            } else {
                *PHREDProb::from(qual)
            };
            // Truncate quality values to PHRED+33 range
            consensus_qual
                .push(cmp::min(93 + offset as u64, (truncated_quality + offset) as u64) as u8);
        }
    }
    fn build_verbose_read_name(
        uuid: &str,
        seq_ids: &[usize],
        read_ids: &Option<HashMap<usize, Vec<u8>>>,
    ) -> String {
        format!(
            "{}_consensus-read-from:{}",
            uuid,
            seq_ids
                .iter()
                .map(|i| String::from_utf8(
                    read_ids
                        .as_ref()
                        .map(|x| x.get(i).unwrap())
                        .unwrap()
                        .to_vec()
                )
                .unwrap())
                .join(",")
        )
    }

    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb;
    fn seqids(&self) -> &'a [usize];
    fn uuid(&self) -> &'a str;
}

#[derive(Debug, Clone)]
pub struct Region {
    pub(crate) target: String,
    pub(crate) start: u64,
    pub(crate) end: u64,
}

impl FromStr for Region {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (target, range) = s.split_once(':').context("No ':' in region string")?;
        let (start, end) = range.split_once('-').context("No '-' in region string")?;
        let start = start.parse::<u64>()?;
        let end = end.parse::<u64>()?;
        Ok(Region {
            target: target.into(),
            start,
            end,
        })
    }
}
