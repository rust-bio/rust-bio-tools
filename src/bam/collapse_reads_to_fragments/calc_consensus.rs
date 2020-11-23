use crate::common::CalcConsensus;
use bio::io::fastq;
use bio::stats::probs::LogProb;
use bio_types::sequence::SequenceRead;
use derive_new::new;
use itertools::Itertools;
use rust_htslib::bam;
use std::collections::HashSet;

const ALLELES: &[u8] = b"ACGT";

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
        let mut consensus_strand = b"S:Z:".to_vec();
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
        let mut strands = HashSet::new();
        let fwd_end_pos = self.recs1()[0].len();
        let rev_start_pos = fwd_end_pos - self.overlap();
        if base_pos < fwd_end_pos {
            self.recs1().iter().for_each(|rec| {
                if rec.base(base_pos) == ref_base {
                    match rec.is_reverse() {
                        true => strands.insert(b'-'),
                        false => strands.insert(b'+'),
                    };
                }
            });
        }
        if base_pos >= rev_start_pos {
            let rev_base_pos = base_pos - rev_start_pos;
            self.recs2().iter().for_each(|rec| {
                if rec.base(rev_base_pos) == ref_base {
                    match rec.is_reverse() {
                        true => strands.insert(b'-'),
                        false => strands.insert(b'+'),
                    };
                }
            });
        }
        match strands.len() == 1 {
            true => consensus_strand.push(strands.take(&(b'-')).unwrap_or(b'+')),
            false => consensus_strand.push(b'*'),
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
        let mut consensus_strand = b"S:Z:".to_vec();
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
        let mut strands = HashSet::new();
        self.recs().iter().for_each(|rec| {
            if rec.base(current_pos) == ref_base {
                match rec.is_reverse() {
                    true => strands.insert(b'-'),
                    false => strands.insert(b'+'),
                };
            }
        });
        match strands.len() == 1 {
            true => consensus_strand.push(strands.take(&(b'-')).unwrap_or(b'+')),
            false => consensus_strand.push(b'*'),
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
