use anyhow::Result;
use bio::alphabets::dna::n_alphabet;
use itertools::Itertools;
use rust_htslib::bcf::{self, Format, Read};

pub fn fix_iupac_alleles() -> Result<()> {
    let mut inbcf = bcf::Reader::from_stdin()?;
    let mut outbcf = bcf::Writer::from_stdout(
        &bcf::Header::from_template(inbcf.header()),
        false,
        Format::Bcf,
    )?;
    let valid_alphabet = n_alphabet();

    for res in inbcf.records() {
        let mut rec = res?;

        let alleles = rec.alleles();
        if !alleles.iter().all(|allele| valid_alphabet.is_word(*allele)) {
            let fixed = alleles
                .into_iter()
                .map(|allele| {
                    allele
                        .iter()
                        .map(|base| {
                            if valid_alphabet.is_word(&[*base]) {
                                *base
                            } else {
                                b'N'
                            }
                        })
                        .collect_vec()
                })
                .collect_vec();

            rec.set_alleles(&fixed.iter().map(|allele| allele.as_slice()).collect_vec())?;
        }

        outbcf.write(&rec)?;
    }

    Ok(())
}
