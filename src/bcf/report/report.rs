use std::collections::HashMap;
use std::error::Error;
use std::io::Write;
use std::{fs, str};

use derive_new::new;
use itertools::Itertools;
use serde_derive::Serialize;
use serde_json;
use tera::{self, Context, Tera};

use chrono::{DateTime, Local};
use jsonm::packer::{PackOptions, Packer};
use rust_htslib::bcf::{self, Read};
use serde_json::{json, Value};
use std::fs::File;
use std::iter::FromIterator;
use std::path::Path;

pub fn oncoprint(
    sample_calls: &HashMap<String, String>,
    output_path: &str,
    vep_annotation: bool,
) -> Result<(), Box<dyn Error>> {
    let mut data = HashMap::new();
    let mut gene_data = HashMap::new();
    let mut unique_genes = HashMap::new();
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

                        let impact = str::from_utf8(fields[2])?;
                        let gene = str::from_utf8(fields[3])?;
                        let dna_alteration = if vep_annotation {
                            str::from_utf8(fields[10])?
                        } else {
                            str::from_utf8(fields[9])?
                        };
                        let protein_alteration = if vep_annotation {
                            str::from_utf8(fields[11])?
                        } else {
                            str::from_utf8(fields[10])?
                        };

                        let gene_rec = unique_genes.entry(gene.to_owned()).or_insert_with(|| 0);
                        *gene_rec += 1;

                        let rec = genes
                            .entry(gene.to_owned())
                            .or_insert_with(|| Record::new(sample.to_owned(), gene.to_owned()));
                        rec.dna_alterations.push(dna_alteration.to_owned());
                        if protein_alteration.len() > 0 {
                            rec.protein_alterations.push(protein_alteration.to_owned());
                        }
                        rec.variants.push(variant.to_owned());
                        rec.impact.push((variant.to_owned(), impact.to_owned()))
                    }
                }
            }
        }

        for gene in genes.keys().sorted() {
            let record = genes.get(gene).unwrap();
            // data for first stage
            let entry = data.entry(gene.to_owned()).or_insert_with(&Vec::new);
            entry.push(FinalRecord::from(record));

            // data for second stage
            let gene_entry = gene_data.entry(gene.to_owned()).or_insert_with(&Vec::new);
            let gene_records: Vec<GeneRecord> = Vec::<GeneRecord>::from(record);
            for rec in gene_records {
                gene_entry.push(rec);
            }
        }
    }

    let gene_path = output_path.to_owned() + "/genes/";
    fs::create_dir(Path::new(&gene_path))?;
    let mut gene_templates = Tera::default();
    gene_templates.register_filter("embed_source", embed_source);
    gene_templates.add_raw_template("genes.html.tera", include_str!("genes.html.tera"))?;

    let gene_specs: Value = serde_json::from_str(include_str!("gene_specs.json")).unwrap();

    for (gene, gene_data) in gene_data {
        let mut specs = gene_specs.clone();
        let values = json!({ "values": gene_data });
        specs["data"] = values;
        let mut packer = Packer::new();
        let options = PackOptions::new();
        let packed_gene_specs = packer.pack(&specs, &options).unwrap();
        let mut context = Context::new();
        context.insert("genespecs", &serde_json::to_string(&packed_gene_specs)?);
        context.insert("gene", &gene);
        let local: DateTime<Local> = Local::now();
        context.insert("time", &local.format("%a %b %e %T %Y").to_string());
        context.insert("version", &env!("CARGO_PKG_VERSION"));
        let html = gene_templates.render("genes.html.tera", &context)?;
        let filepath = String::from(gene_path.clone()) + &gene + ".html";
        let mut file = File::create(filepath)?;
        file.write_all(html.as_bytes())?;
    }

    // only keep recurrent entries
    let data: Vec<_> = data
        .values()
        .filter(|entry| entry.len() >= 1)
        .flatten()
        .sorted()
        .collect();

    let page_size = 100;

    let mut v = Vec::from_iter(unique_genes);
    v.sort_by(|&(_, a), &(_, b)| b.cmp(&a));

    let pages = v.len() / page_size;

    for i in 0..pages + 1 {
        let current_genes = if i != pages {
            &v[(i * page_size)..((i + 1) * page_size)] // get genes for current page
        } else {
            &v[(i * page_size)..] // get genes for last page
        };

        let mut sorted_genes = Vec::new();
        for (g, _) in current_genes {
            sorted_genes.push(g);
        }

        let page = i + 1;

        let page_data: Vec<_> = data
            .iter()
            .filter(|entry| sorted_genes.contains(&&entry.gene))
            .collect();

        let mut vl_specs: Value = serde_json::from_str(include_str!("report_specs.json")).unwrap();
        let values = json!({ "values": page_data });
        vl_specs["data"] = values;

        let mut packer = Packer::new();
        let options = PackOptions::new();
        let packed_specs = packer.pack(&vl_specs, &options).unwrap();
        let mut templates = Tera::default();
        templates.register_filter("embed_source", embed_source);
        templates.add_raw_template("report.html.tera", include_str!("report.html.tera"))?;
        let mut context = Context::new();
        let data = serde_json::to_string(&packed_specs)?;
        context.insert("oncoprint", &data);
        context.insert("pages", &(pages + 1));
        let local: DateTime<Local> = Local::now();
        context.insert("time", &local.format("%a %b %e %T %Y").to_string());
        context.insert("version", &env!("CARGO_PKG_VERSION"));

        let html = templates.render("report.html.tera", &context)?;

        let index = output_path.to_owned() + "/index" + &page.to_string() + ".html";
        let mut file = File::create(index)?;
        file.write_all(html.as_bytes())?;
    }

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

#[derive(new, Debug)]
struct Record {
    sample: String,
    gene: String,
    #[new(default)]
    dna_alterations: Vec<String>,
    #[new(default)]
    protein_alterations: Vec<String>,
    #[new(default)]
    variants: Vec<String>,
    #[new(default)]
    impact: Vec<(String, String)>,
}

// Record for the first stage of the report
#[derive(Serialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct FinalRecord {
    sample: String,
    gene: String,
    dna_alterations: String,
    protein_alterations: String,
    variants: String,
    impact: String,
}

// Record for the second stage of the report
#[derive(Serialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct GeneRecord {
    sample: String,
    gene: String,
    alteration: String,
    variants: String,
}

impl From<&Record> for Vec<GeneRecord> {
    fn from(record: &Record) -> Self {
        let mut gene_vec = Vec::new();
        for dna_alt in record.dna_alterations.iter().sorted().unique() {
            let alt = GeneRecord {
                sample: record.sample.to_owned(),
                gene: record.gene.to_owned(),
                alteration: dna_alt.to_owned(),
                variants: record.variants.iter().sorted().unique().join("/"),
            };
            gene_vec.push(alt);
        }

        for protein_alt in record.protein_alterations.iter().sorted().unique() {
            let alt = GeneRecord {
                sample: record.sample.to_owned(),
                gene: record.gene.to_owned(),
                alteration: protein_alt.to_owned(),
                variants: record.variants.iter().sorted().unique().join("/"),
            };
            gene_vec.push(alt);
        }

        gene_vec
    }
}

impl From<&Record> for FinalRecord {
    fn from(record: &Record) -> Self {
        let mut imp_map = HashMap::new();
        let mut impact_vec = Vec::new();
        for (variant, imp) in &record.impact {
            let rec = imp_map.entry(variant.to_owned()).or_insert_with(|| 0);
            if imp.to_uppercase() == "HIGH" {
                *rec += 1;
            }
        }

        for (variant, impact) in imp_map {
            let i = variant + "=" + &impact.to_string();
            impact_vec.push(i);
        }

        FinalRecord {
            sample: record.sample.to_owned(),
            gene: record.gene.to_owned(),
            dna_alterations: record.dna_alterations.iter().sorted().unique().join(", "),
            protein_alterations: record
                .protein_alterations
                .iter()
                .sorted()
                .unique()
                .join(", "),
            variants: record.variants.iter().sorted().unique().join("/"),
            impact: impact_vec.iter().join("/"),
        }
    }
}
