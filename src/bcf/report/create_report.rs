use crate::bcf::report::fasta_reader::{get_fasta_length, read_fasta};
use crate::bcf::report::static_reader::{get_static_reads, Variant};
use rust_htslib::bcf::Read;
use rustc_serialize::json::Json;
use serde::Serialize;
use serde_json::{json, Value};
use std::error::Error;
use std::path::Path;

#[derive(Serialize, Clone, Debug, PartialEq)]
pub enum VariantType {
    Deletion,
    Insertion,
    Duplicate,
    Inversion,
    Variant,
}

#[derive(Serialize, Clone, Debug, PartialEq)]
pub struct Report {
    id: String,
    name: String,
    position: i64,
    reference: String,
    var_type: VariantType,
    alternatives: Option<String>,
    ann: Option<Vec<Vec<String>>>,
    vis: String,
}

pub(crate) fn make_report(
    vcf_path: &Path,
    fasta_path: &Path,
    bam_path: &Path,
) -> Result<Vec<Report>, Box<dyn Error>> {
    let mut vcf = rust_htslib::bcf::Reader::from_path(&vcf_path).unwrap();
    let header = vcf.header().clone();

    let mut reports = Vec::new();

    for v in vcf.records() {
        let mut variant = v.unwrap();

        let n = header.rid2name(variant.rid().unwrap()).unwrap().to_owned();
        let i = variant.id();

        let chrom = String::from_utf8(n).unwrap();

        let id = String::from_utf8(i).unwrap();

        let pos = variant.pos();
        let end_pos = match variant.info(b"END").integer() {
            Ok(Some(end_pos)) => {
                // Subtraction of 0.5 because of the 0-based positioning in the whole plot
                let end_pos = end_pos[0] as f64 - 0.5; // -1 due to 0-basing, + 0.5 du to end pos
                Some(end_pos)
            }
            _ => None,
        };

        let mut ann_strings = Vec::new();
        let ann = variant.info(b"ANN").string()?;

        if let Some(ann) = ann {
            for entry in ann {
                let fields: Vec<_> = entry.split(|c| *c == b'|').collect();
                let mut s = Vec::new();
                for f in fields {
                    let attr = String::from_utf8(f.to_owned()).unwrap();
                    s.push(attr);
                }
                ann_strings.push(s);
            }
        }

        let alleles = variant.alleles();

        if alleles.len() > 0 {
            let ref_vec = alleles[0].to_owned();
            let ref_allele = String::from_utf8(ref_vec).unwrap();

            let len: u8 = ref_allele.len() as u8;

            for i in 1..alleles.len() {
                let alt = alleles[i];
                let var_string = String::from("Variant");
                let var_type: VariantType;
                let alternatives: Option<String>;
                let end_position: f64;
                let plot_start_position;

                match alt {
                    b"<DEL>" => {
                        var_type = VariantType::Deletion;
                        alternatives = None;
                        end_position = end_pos.unwrap();
                        plot_start_position = pos as f64 - 0.5;
                    }
                    b"<INV>" => {
                        var_type = VariantType::Inversion;
                        let rev: String = ref_allele.chars().rev().collect();
                        alternatives = Some(rev.clone());
                        end_position = end_pos.unwrap();
                        plot_start_position = pos as f64 - 0.5;
                    }
                    b"<DUP>" => {
                        var_type = VariantType::Duplicate;
                        let dup: String = [ref_allele.clone(), ref_allele.clone()].concat();
                        alternatives = Some(dup.clone());
                        end_position = end_pos.unwrap();
                        plot_start_position = pos as f64 - 0.5;
                    }
                    _ => {
                        let mut alt_allele = String::from("");

                        for c in alt {
                            if *c as char != '<' && *c as char != '>' {
                                alt_allele.push(*c as char);
                            }
                        }

                        match alt_allele.len() {
                            a if a < ref_allele.len() => {
                                plot_start_position = pos as f64 + 0.5; // start position + 1 due to alignment with deletions from bam (example: ref: ACTT alt: A  -> deletion is just CTT)
                                end_position = pos as f64 - 0.5 + len as f64;
                                var_type = VariantType::Deletion;
                                alternatives = Some(alt_allele.clone());
                            }
                            a if a > ref_allele.len() => {
                                plot_start_position = pos as f64;
                                end_position = pos as f64 + len as f64;
                                var_type = VariantType::Insertion;
                                alternatives = Some(alt_allele.clone());
                            }
                            _ => {
                                plot_start_position = pos as f64 - 0.5;
                                end_position = pos as f64 - 0.5 + len as f64;
                                var_type = VariantType::Variant;
                                alternatives = Some(alt_allele.clone());
                            }
                        }
                    }
                }

                let var = Variant {
                    marker_type: var_string,
                    reference: ref_allele.clone(),
                    alternatives: alternatives,
                    start_position: plot_start_position,
                    end_position: end_position,
                    row: -1,
                    var_type: var_type,
                };

                let visualization: Value;
                let fasta_length = get_fasta_length(fasta_path);

                if variant.pos() < 75 {
                    let content = create_report_data(
                        fasta_path,
                        var.clone(),
                        bam_path,
                        chrom.clone(),
                        0,
                        end_position as u64 + 75,
                    );
                    visualization = manipulate_json(content, 0, end_position as u64 + 75);
                } else if variant.pos() + 75 >= fasta_length as i64 {
                    let content = create_report_data(
                        fasta_path,
                        var.clone(),
                        bam_path,
                        chrom.clone(),
                        variant.pos() as u64 - 75,
                        fasta_length - 1,
                    );
                    visualization =
                        manipulate_json(content, variant.pos() as u64 - 75, fasta_length - 1);
                } else {
                    let content = create_report_data(
                        fasta_path,
                        var.clone(),
                        bam_path,
                        chrom.clone(),
                        variant.pos() as u64 - 75,
                        end_position as u64 + 75,
                    );
                    visualization = manipulate_json(
                        content,
                        variant.pos() as u64 - 75,
                        end_position as u64 + 75,
                    );
                }

                let r = Report {
                    id: id.clone(),
                    name: chrom.clone(),
                    position: variant.pos(),
                    reference: ref_allele.clone(),
                    var_type: var.var_type,
                    alternatives: var.alternatives,
                    ann: Some(ann_strings.clone()),
                    vis: visualization.to_string(),
                };

                reports.push(r);
            }
        }
    }

    Ok(reports.clone())
}

fn create_report_data(
    fasta_path: &Path,
    variant: Variant,
    bam_path: &Path,
    chrom: String,
    from: u64,
    to: u64,
) -> Json {
    let mut data = Vec::new();

    for f in read_fasta(fasta_path.clone(), chrom.clone(), from, to) {
        let nucleobase = json!(f);
        data.push(nucleobase);
    }

    let (bases, matches) = get_static_reads(bam_path, fasta_path, chrom.clone(), from, to);

    for b in bases {
        let base = json!(b);
        data.push(base);
    }

    for m in matches {
        let mat = json!(m);
        data.push(mat);
    }

    data.push(json!(variant));

    let values = Json::from_str(&json!(data).to_string()).unwrap();

    values
}

/// Inserts the json containing the genome data into the vega specs.
/// It also changes keys and values of the json data for the vega plot to look better.
fn manipulate_json(data: Json, from: u64, to: u64) -> Value {
    let json_string = include_str!("vegaSpecs.json");

    let mut vega_specs: Value = serde_json::from_str(&json_string).unwrap();
    let values: Value = serde_json::from_str(&data.to_string()).unwrap();
    let mut values = json!({"values": values, "name": "fasta"});

    let v = values["values"].as_array().unwrap().clone();

    for i in 0..v.len() {
        let k = v[i]["marker_type"].clone().as_str().unwrap().to_owned();

        if k == "A" || k == "T" || k == "G" || k == "C" || k == "U" {
            values["values"][i]["base"] = values["values"][i]["marker_type"].clone();
        } else if k == "Deletion"
            || k == "Match"
            || k == "Pairing"
            || k == "Duplicate"
            || k == "Inversion"
        {
            values["values"][i]["typ"] = values["values"][i]["marker_type"].clone();
        } else if k == "Insertion" {
            values["values"][i]["typ"] = values["values"][i]["marker_type"].clone();
            values["values"][i]["inserts"] = values["values"][i]["bases"].clone();
        }
    }

    vega_specs["width"] = json!(700);
    let domain = json!([from, to]);

    vega_specs["scales"][0]["domain"] = domain;
    vega_specs["data"][1] = values;

    vega_specs
}
