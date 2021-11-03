use crate::common::CalcConsensus;
use bio::io::fastq;
use bio::stats::probs::LogProb;
use bio_types::sequence::SequenceRead;
use bio_types::sequence::SequenceReadPairOrientation;
use derive_new::new;
use itertools::Itertools;
use rust_htslib::bam;
use std::collections::{HashMap, HashSet};
use std::ops::BitOrAssign;

const ALLELES: &[u8] = b"ACGT";

#[derive(Eq, PartialEq)]
enum StrandObservation {
    None,
    Forward,
    Reverse,
    Both,
}

impl BitOrAssign for StrandObservation {
    fn bitor_assign(&mut self, rhs: Self) {
        if let StrandObservation::None = self {
            *self = rhs;
        } else if *self != rhs {
            *self = StrandObservation::Both;
        }
    }
}

#[derive(new)]
pub struct CalcOverlappingConsensus<'a> {
    recs1: &'a [bam::Record],
    recs2: &'a [bam::Record],
    r1_vec: &'a [bool],
    r2_vec: &'a [bool],
    seqids: &'a [usize],
    uuid: &'a str,
    read_ids: &'a mut Option<HashMap<usize, Vec<u8>>>,
}

impl<'a> CalcOverlappingConsensus<'a> {
    pub fn calc_consensus(&self) -> (fastq::Record, LogProb) {
        let seq_len = self.r1_vec().len();
        let mut consensus_seq: Vec<u8> = Vec::with_capacity(seq_len);
        let mut consensus_qual: Vec<u8> = Vec::with_capacity(seq_len);
        let mut consensus_strand = b"SI:Z:".to_vec();
        let read_orientations_opt = self.build_read_orientation_string();
        let mut consensus_lh = LogProb::ln_one();
        for i in 0..seq_len {
            match (
                self.recs1().len() == 1,
                self.map_read_pos(i, self.r1_vec()),
                self.map_read_pos(i, self.r2_vec()),
            ) {
                (true, Some(base_pos), None) => {
                    let base = self.recs1()[0].seq().as_bytes()[base_pos];
                    consensus_seq.push(base);
                    consensus_qual.push(self.recs1()[0].qual()[base_pos] + 33);
                    consensus_lh += Self::overall_allele_likelihood(self, &base, i);
                }
                (true, None, Some(base_pos)) => {
                    let base = self.recs2()[0].seq().as_bytes()[base_pos];
                    consensus_seq.push(base);
                    consensus_qual.push(self.recs2()[0].qual()[base_pos] + 33);
                    consensus_lh += Self::overall_allele_likelihood(self, &base, i);
                }
                _ => {
                    let likelihoods = ALLELES
                        .iter()
                        .map(|a| Self::overall_allele_likelihood(self, a, i))
                        .collect_vec();
                    Self::build_consensus_sequence(
                        likelihoods,
                        &mut consensus_lh,
                        &mut consensus_seq,
                        &mut consensus_qual,
                        33.0,
                    );
                }
            };
            self.build_consensus_strand(&mut consensus_strand, consensus_seq[i], i);
        }
        let name = if self.read_ids.is_some() {
            format!(
                "{}_consensus-read-from:{}",
                self.uuid(),
                self.seqids()
                    .iter()
                    .map(|i| String::from_utf8(
                        self.read_ids
                            .as_ref()
                            .map(|x| x.get(i).unwrap())
                            .unwrap()
                            .to_vec()
                    )
                    .unwrap())
                    .join(",")
            )
        } else {
            format!(
                "{}_consensus-read-from:{}_reads",
                self.uuid(),
                self.seqids().len(),
            )
        };
        if let Some(mut read_orientations) = read_orientations_opt {
            consensus_strand.append(&mut read_orientations)
        }
        let consensus_rec = fastq::Record::with_attrs(
            &name,
            Some(&String::from_utf8(consensus_strand).unwrap()),
            &consensus_seq,
            &consensus_qual,
        );
        (consensus_rec, consensus_lh)
    }

    fn recs1(&self) -> &[bam::Record] {
        self.recs1
    }

    fn recs2(&self) -> &[bam::Record] {
        self.recs2
    }

    fn r1_vec(&self) -> &[bool] {
        self.r1_vec
    }

    fn r2_vec(&self) -> &[bool] {
        self.r2_vec
    }

    fn build_consensus_strand(&self, consensus_strand: &mut Vec<u8>, ref_base: u8, pos: usize) {
        let mut strand = StrandObservation::None;
        let rec1_pos = self.map_read_pos(pos, self.r1_vec());
        let rec2_pos = self.map_read_pos(pos, self.r2_vec());
        let mut strand_observation = |recs: &[bam::Record], rec_pos: Option<usize>| {
            if let Some(pos) = rec_pos {
                recs.iter().for_each(|rec| {
                    if rec.base(pos) == ref_base {
                        match rec.is_reverse() {
                            true => strand |= StrandObservation::Reverse,
                            false => strand |= StrandObservation::Forward,
                        };
                    }
                });
            }
        };
        strand_observation(self.recs1(), rec1_pos);
        strand_observation(self.recs2(), rec2_pos);
        match strand {
            StrandObservation::Forward => consensus_strand.push(b'+'),
            StrandObservation::Reverse => consensus_strand.push(b'-'),
            StrandObservation::Both => consensus_strand.push(b'*'),
            StrandObservation::None => {
                unreachable!()
            }
        }
    }
    fn build_read_orientation_string(&self) -> Option<Vec<u8>> {
        let mut read_orientations_set: HashSet<_> = self
            .recs1()
            .iter()
            .filter_map(|rec| match rec.read_pair_orientation() {
                SequenceReadPairOrientation::F2F1 => Some(b"F2F1,"),
                SequenceReadPairOrientation::F2R1 => Some(b"F2R1,"),
                SequenceReadPairOrientation::F1F2 => Some(b"F1F2,"),
                SequenceReadPairOrientation::R2F1 => Some(b"R2F1,"),
                SequenceReadPairOrientation::F1R2 => Some(b"F1R2,"),
                SequenceReadPairOrientation::R2R1 => Some(b"R2R1,"),
                SequenceReadPairOrientation::R1F2 => Some(b"R1F2,"),
                SequenceReadPairOrientation::R1R2 => Some(b"R1R2,"),
                SequenceReadPairOrientation::None => None,
            })
            .collect();
        let mut read_orientations_string = b" RO:Z:".to_vec();
        read_orientations_set
            .drain()
            .for_each(|entry| read_orientations_string.extend_from_slice(entry));
        match read_orientations_string.pop() {
            Some(b',') => Some(read_orientations_string),
            Some(b':') => None,
            Some(_) => unreachable!(),
            None => unreachable!(),
        }
    }
    fn map_read_pos(&self, consensus_pos: usize, alignment_vec: &[bool]) -> Option<usize> {
        match alignment_vec[consensus_pos] {
            true => Some(
                alignment_vec[0..(consensus_pos + 1)]
                    .iter()
                    .filter(|&v| *v)
                    .count()
                    - 1,
            ),
            false => None,
        }
    }
}

impl<'a> CalcConsensus<'a, bam::Record> for CalcOverlappingConsensus<'a> {
    fn overall_allele_likelihood(&self, allele: &u8, pos: usize) -> LogProb {
        let mut lh = LogProb::ln_one();
        let rec1_pos = self.map_read_pos(pos, self.r1_vec());
        let rec2_pos = self.map_read_pos(pos, self.r2_vec());
        for (rec1, rec2) in self.recs1().iter().zip(self.recs2()) {
            if let Some(pos) = rec1_pos {
                lh += Self::allele_likelihood_in_rec(
                    allele,
                    &rec1.seq().as_bytes(),
                    rec1.qual(),
                    pos,
                    0,
                );
            };
            if let Some(pos) = rec2_pos {
                lh += Self::allele_likelihood_in_rec(
                    allele,
                    &rec2.seq().as_bytes(),
                    rec2.qual(),
                    pos,
                    0,
                );
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

#[derive(new)]
pub struct CalcNonOverlappingConsensus<'a> {
    recs: &'a [bam::Record],
    seqids: &'a [usize],
    uuid: &'a str,
    read_ids: &'a mut Option<HashMap<usize, Vec<u8>>>,
}

impl<'a> CalcNonOverlappingConsensus<'a> {
    pub fn calc_consensus(&self) -> (fastq::Record, LogProb) {
        let seq_len = self.recs()[0].seq().len();
        let mut consensus_seq: Vec<u8> = Vec::with_capacity(seq_len);
        let mut consensus_qual: Vec<u8> = Vec::with_capacity(seq_len);
        let mut consensus_strand = b"SI:Z:".to_vec();
        let mut cigar_map = HashMap::new();
        for record in self.recs() {
            let cached_cigar = record.raw_cigar();
            if !cigar_map.contains_key(cached_cigar) {
                cigar_map.insert(cached_cigar, Vec::new());
            }
            cigar_map.get_mut(cached_cigar).unwrap().push(record);
        }

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
                .map(|a| Self::overall_allele_likelihood(self, a, i))
                .collect_vec(); //Check this. See below
            Self::build_consensus_sequence(
                likelihoods,
                &mut consensus_lh,
                &mut consensus_seq,
                &mut consensus_qual,
                33.0,
            );
            self.build_consensus_strand(&mut consensus_strand, consensus_seq[i], i);
        }
        let name = if self.read_ids.is_some() {
            format!(
                "{}_consensus-read-from:{}",
                self.uuid(),
                self.seqids()
                    .iter()
                    .map(|i| String::from_utf8(
                        self.read_ids
                            .as_ref()
                            .map(|x| x.get(i).unwrap())
                            .unwrap()
                            .to_vec()
                    )
                    .unwrap())
                    .join(",")
            )
        } else {
            format!(
                "{}_consensus-read-from:{}_reads",
                self.uuid(),
                self.seqids().len(),
            )
        };
        let consensus_rec = fastq::Record::with_attrs(
            &name,
            Some(&String::from_utf8(consensus_strand).unwrap()),
            &consensus_seq,
            &consensus_qual,
        );
        (consensus_rec, consensus_lh)
    }
    pub fn recs(&self) -> &[bam::Record] {
        self.recs
    }
    fn build_consensus_strand(
        &self,
        consensus_strand: &mut Vec<u8>,
        ref_base: u8,
        current_pos: usize,
    ) {
        let mut strand = StrandObservation::None;
        self.recs().iter().for_each(|rec| {
            if rec.base(current_pos) == ref_base {
                match rec.is_reverse() {
                    true => strand |= StrandObservation::Reverse,
                    false => strand |= StrandObservation::Forward,
                };
            }
        });
        match strand {
            StrandObservation::Forward => consensus_strand.push(b'+'),
            StrandObservation::Reverse => consensus_strand.push(b'-'),
            StrandObservation::Both => consensus_strand.push(b'*'),
            StrandObservation::None => consensus_strand.push(b'.'),
        }
    }
}

impl<'a> CalcConsensus<'a, bam::Record> for CalcNonOverlappingConsensus<'a> {
    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb {
        let mut lh = LogProb::ln_one(); // posterior: log(P(theta)) = 1
        for rec in self.recs() {
            lh += Self::allele_likelihood_in_rec(allele, &rec.seq().as_bytes(), rec.qual(), i, 0);
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
