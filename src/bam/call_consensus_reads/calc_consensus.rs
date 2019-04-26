use crate::common_functions::CalcConsensus;
use bio::stats::probs::LogProb;
use itertools::Itertools;
use rust_htslib::bam;


const ALLELES: &'static [u8] = b"ACGT";

pub struct CalcOverlappingConsensus<'a> {
    recs1: &'a [bam::Record],
    recs2: &'a [bam::Record],
    overlap: usize,
    seqids: &'a [usize],
    uuid: &'a str,
}

impl<'a> CalcOverlappingConsensus<'a> {
    pub fn new(
        recs1: &'a [bam::Record],
        recs2: &'a [bam::Record],
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
                0.0,
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

impl<'a> CalcConsensus<'a, bam::Record> for CalcOverlappingConsensus<'a> {
    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb {
        let mut lh = LogProb::ln_one();
        for (rec1, rec2) in self.recs1().into_iter().zip(self.recs2()) {
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
    fn uuid(&self) -> &'a str {
        self.uuid
    }
    fn seqids(&self) -> &'a [usize] {
        self.seqids
    }
}

pub struct CalcNonOverlappingConsensus<'a> {
    recs: &'a [bam::Record],
    seqids: &'a [usize],
    uuid: &'a str,
}

impl<'a> CalcNonOverlappingConsensus<'a> {
    pub fn new(recs: &'a [bam::Record], seqids: &'a [usize], uuid: &'a str) -> Self {
        CalcNonOverlappingConsensus { recs, seqids, uuid }
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
                0.0,
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

impl<'a> CalcConsensus<'a, bam::Record> for CalcNonOverlappingConsensus<'a> {
    fn overall_allele_likelihood(&self, allele: &u8, i: usize) -> LogProb {
        let mut lh = LogProb::ln_one(); // posterior: log(P(theta)) = 1
        for rec in self.recs() {
            lh += Self::allele_likelihood_in_rec(allele, &rec.seq().as_bytes(), rec.qual(), i, 0);
        }
        lh
    }
    fn uuid(&self) -> &'a str {
        self.uuid
    }
    fn seqids(&self) -> &'a [usize] {
        self.seqids
    }
}
