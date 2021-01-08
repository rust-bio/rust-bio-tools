use crate::common::CalcConsensus;
use bio::io::fastq;
use bio::stats::probs::LogProb;
use bio_types::sequence::SequenceRead;
use bio_types::sequence::SequenceReadPairOrientation;
use derive_new::new;
use itertools::Itertools;
use rust_htslib::bam;
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
    overlap: usize,
    seqids: &'a [usize],
    uuid: &'a str,
    verbose_read_names: bool,
}

impl<'a> CalcOverlappingConsensus<'a> {
    pub fn calc_consensus(&self) -> (fastq::Record, LogProb) {
        let seq_len = self.recs1()[0].seq().len() + self.recs2()[0].seq().len() - self.overlap();
        let mut consensus_seq: Vec<u8> = Vec::with_capacity(seq_len);
        let mut consensus_qual: Vec<u8> = Vec::with_capacity(seq_len);
        let mut consensus_strand = b"SI:Z:".to_vec();
        let read_orientations_opt = self.build_read_orientation_string();
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
            self.build_consensus_strand(&mut consensus_strand, consensus_seq[i], i);
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

    fn overlap(&self) -> usize {
        self.overlap
    }
    fn build_consensus_strand(
        &self,
        consensus_strand: &mut Vec<u8>,
        ref_base: u8,
        base_pos: usize,
    ) {
        let mut strand = StrandObservation::None;
        let first_end_pos = self.recs1()[0].len();
        let second_start_pos = first_end_pos - self.overlap();
        if base_pos < first_end_pos {
            self.recs1().iter().for_each(|rec| {
                if rec.base(base_pos) == ref_base {
                    match rec.is_reverse() {
                        true => strand |= StrandObservation::Reverse,
                        false => strand |= StrandObservation::Forward,
                    };
                }
            });
        }
        if base_pos >= second_start_pos {
            let second_base_pos = base_pos - second_start_pos;
            self.recs2().iter().for_each(|rec| {
                if rec.base(second_base_pos) == ref_base {
                    match rec.is_reverse() {
                        true => strand |= StrandObservation::Reverse,
                        false => strand |= StrandObservation::Forward,
                    };
                }
            });
        }
        match strand {
            StrandObservation::Forward => consensus_strand.push(b'+'),
            StrandObservation::Reverse => consensus_strand.push(b'-'),
            StrandObservation::Both => consensus_strand.push(b'*'),
            StrandObservation::None => {
                dbg!(&base_pos);
                dbg!(String::from_utf8([ref_base].to_vec()).unwrap());
                dbg!(&self.overlap());
                dbg!(&self.recs1()[0].pos());
                dbg!(&self.recs1()[0].len());
                dbg!(&self.recs1()[0].cigar_cached().unwrap().end_pos());
                dbg!(&self.recs1()[0].cigar_cached().unwrap().into_iter());
                dbg!(String::from_utf8([self.recs1()[0].base(base_pos)].to_vec()).unwrap());
                dbg!(&self.recs2()[0].pos());
                dbg!(&self.recs2()[0].len());
                dbg!(&self.recs2()[0].cigar_cached().unwrap().end_pos());
                dbg!(&self.recs2()[0].cigar_cached().unwrap().into_iter());
                dbg!(&self.recs2()[0].seq());
                unreachable!()
            }
        }
    }
    fn build_read_orientation_string(&self) -> Option<Vec<u8>> {
        let mut read_orientations = b" RO:Z:".to_vec();
        self.recs1()
            .iter()
            .for_each(|rec| match rec.read_pair_orientation() {
                SequenceReadPairOrientation::F2F1 => read_orientations.extend_from_slice(b"F2F1,"),
                SequenceReadPairOrientation::F2R1 => read_orientations.extend_from_slice(b"F2R1,"),
                SequenceReadPairOrientation::F1F2 => read_orientations.extend_from_slice(b"F1F2,"),
                SequenceReadPairOrientation::R2F1 => read_orientations.extend_from_slice(b"R2F1,"),
                SequenceReadPairOrientation::F1R2 => read_orientations.extend_from_slice(b"F1R2,"),
                SequenceReadPairOrientation::R2R1 => read_orientations.extend_from_slice(b"R2R1,"),
                SequenceReadPairOrientation::R1F2 => read_orientations.extend_from_slice(b"R1F2,"),
                SequenceReadPairOrientation::R1R2 => read_orientations.extend_from_slice(b"R1R2,"),
                SequenceReadPairOrientation::None => unreachable!(),
            });
        match read_orientations.pop() {
            Some(b',') => Some(read_orientations),
            Some(b':') => None,
            Some(_) => unreachable!(),
            None => unreachable!(),
        }
    }
}

impl<'a> CalcConsensus<'a, bam::Record> for CalcOverlappingConsensus<'a> {
    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb {
        let mut lh = LogProb::ln_one();
        for (rec1, rec2) in self.recs1().iter().zip(self.recs2()) {
            if i < rec1.seq().len() {
                lh += Self::allele_likelihood_in_rec(
                    allele,
                    &rec1.seq().as_bytes(),
                    rec1.qual(),
                    i,
                    0,
                );
            };
            if i >= rec1.seq().len() - self.overlap() {
                let rec2_i = i - (rec1.seq().len() - self.overlap());
                lh += Self::allele_likelihood_in_rec(
                    allele,
                    &rec2.seq().as_bytes(),
                    &rec2.qual(),
                    rec2_i,
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
}

impl<'a> CalcNonOverlappingConsensus<'a> {
    pub fn calc_consensus(&self) -> (fastq::Record, LogProb) {
        let seq_len = self.recs()[0].seq().len();
        let mut consensus_seq: Vec<u8> = Vec::with_capacity(seq_len);
        let mut consensus_qual: Vec<u8> = Vec::with_capacity(seq_len);
        let mut consensus_strand = b"SI:Z:".to_vec();
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
            self.build_consensus_strand(&mut consensus_strand, consensus_seq[i], i);
        }
        let consensus_rec = fastq::Record::with_attrs(
            &self.uuid(),
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
            StrandObservation::None => unreachable!(),
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
