//! Compute the B-allele frequencies for a given VCF file.
//!
//! ## Usage:
//! ```bash
//! $ rbt vcf-baf < tests/test-freebayes.vcf > tests/baf.
//! ```
//!
use itertools::repeat_n;
use itertools::Itertools;
use rust_htslib::bcf;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::Read;
use std::f32;

use crate::errors;
use snafu::ResultExt;

pub fn calculate_baf() -> errors::Result<()> {
    let mut reader = bcf::Reader::from_stdin().context(errors::BCFReaderStdinError {})?;

    let mut header = bcf::Header::from_template(reader.header());
    header.push_record(b"##FORMAT=<ID=BAF,Number=A,Type=Float,Description=\"b-allele frequency\">");

    let mut writer = bcf::Writer::from_stdout(&header, false, bcf::Format::BCF).context(
        errors::BCFWriterStdoutError {
            header: format!("{:?}", header),
        },
    )?;

    for record in reader.records() {
        let mut record = record.context(errors::BCFReadError {
            header: format!("{:?}", header),
        })?;

        let allele_lens = record.alleles().iter().map(|a| a.len()).collect_vec();
        let mut bafs = Vec::new();
        {
            let ref_depths = record
                .format(b"RO")
                .integer()
                .context(errors::BCFFormatReadError)?
                .into_iter()
                .map(|d| d.to_owned())
                .collect_vec();
            let alt_depths = record
                .format(b"AO")
                .integer()
                .context(errors::BCFFormatReadError)?;

            for (sample_ref_depth, sample_alt_depth) in ref_depths.iter().zip(alt_depths.iter()) {
                if allele_lens[0] != 1 || sample_ref_depth[0].is_missing() {
                    bafs.extend(repeat_n(f32::missing(), allele_lens.len() - 1));
                } else {
                    let total_depth =
                        sample_ref_depth[0] + sample_alt_depth.into_iter().sum::<i32>();
                    bafs.extend(allele_lens[1..].iter().zip(sample_alt_depth.iter()).map(
                        |(alen, d)| {
                            if *alen == 1 {
                                *d as f32 / total_depth as f32
                            } else {
                                f32::missing()
                            }
                        },
                    ));
                };
            }
        }

        writer.translate(&mut record);
        record
            .push_format_float(b"BAF", &bafs)
            .context(errors::BCFTagWriteError {
                data: format!("{:?}", bafs),
                fd: String::from("BAF"),
            })?;
        writer.write(&record).context(errors::BCFWriteError)?;
    }

    Ok(())
}
