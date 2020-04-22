use std::collections::HashMap;
use std::error::Error;
use std::io::stdout;
use std::io::Write;
use std::str;

use derive_new::new;
use itertools::Itertools;
use serde_derive::Serialize;
use serde_json;
use tera::{self, Context, Tera};

use rust_htslib::bcf::{self, Read};

pub fn oncoprint(sample_calls: &HashMap<String, String>) -> Result<(), Box<dyn Error>> {
    let mut data = Vec::new();
    for (sample, path) in sample_calls.iter().sorted() {
        let mut genes = HashMap::new();
        let mut bcf_reader = bcf::Reader::from_path(path)?;

        for res in bcf_reader.records() {
            let mut record = res?;
            let alleles = record
                .alleles()
                .into_iter()
                .map(|allele| allele.to_owned())
                .collect_vec();
            let alt_alleles = &alleles[1..];
            let ref_allele = alleles[0].to_owned();
            let ann = record.info(b"ANN").string()?;

            if let Some(ann) = ann {
                for alt_allele in alt_alleles {
                    let variant = if alt_allele == b"<DEL>" {
                        "DEL"
                    } else if alt_allele == b"<BND>" {
                        "BND"
                    } else if alt_allele == b"<INV>" {
                        "INV"
                    } else if alt_allele == b"<DUP>" {
                        "DUP"
                    } else if alt_allele.len() == 1 && ref_allele.len() == 1 {
                        "SNV"
                    } else if alt_allele.len() == ref_allele.len() {
                        "MNV"
                    } else {
                        "Complex"
                    };

                    for entry in &ann {
                        let fields: Vec<_> = entry.split(|c| *c == b'|').collect();
                        let alt = &fields[0];
                        if alt != &alt_allele.as_slice() {
                            continue;
                        }

                        let gene = str::from_utf8(fields[3])?;
                        let dna_alteration = str::from_utf8(fields[9])?;
                        let protein_alteration = str::from_utf8(fields[10])?;

                        let rec = genes
                            .entry(gene.to_owned())
                            .or_insert_with(|| Record::new(sample.to_owned(), gene.to_owned()));
                        rec.dna_alterations.push(dna_alteration.to_owned());
                        if protein_alteration.len() > 0 {
                            rec.protein_alterations.push(protein_alteration.to_owned());
                        }
                        rec.variants.push(variant.to_owned());
                    }
                }
            }
        }

        for gene in genes.keys().sorted() {
            let record = genes.get(gene).unwrap();
            data.push(FinalRecord::from(record));
        }
    }

    let mut templates = Tera::new("templates/*.html.tera")?;
    templates.register_filter("embed_source", embed_source);
    let mut context = Context::new();
    let data = serde_json::to_string(&data)?;
    context.insert("data", &data);

    let html = templates.render("oncoprint.html.tera", &context)?;

    stdout().write(html.as_bytes())?;

    Ok(())
}

fn embed_source(
    value: &tera::Value,
    _: &HashMap<String, tera::Value>,
) -> tera::Result<tera::Value> {
    let url = tera::try_get_value!("upper", "value", String, value);
    let source = reqwest::get(&url).unwrap().text().unwrap();
    Ok(tera::to_value(source).unwrap())
}

#[derive(new)]
struct Record {
    sample: String,
    gene: String,
    #[new(default)]
    dna_alterations: Vec<String>,
    #[new(default)]
    protein_alterations: Vec<String>,
    #[new(default)]
    variants: Vec<String>,
}

#[derive(Serialize)]
struct FinalRecord {
    sample: String,
    gene: String,
    dna_alterations: String,
    protein_alterations: String,
    variants: String,
}

impl From<&Record> for FinalRecord {
    fn from(record: &Record) -> Self {
        FinalRecord {
            sample: record.sample.to_owned(),
            gene: record.gene.to_owned(),
            dna_alterations: record
                .dna_alterations
                .iter()
                .sorted()
                .iter()
                .unique()
                .join(", "),
            protein_alterations: record
                .protein_alterations
                .iter()
                .sorted()
                .iter()
                .unique()
                .join(", "),
            variants: record.variants.iter().sorted().iter().unique().join("/"),
        }
    }
}
