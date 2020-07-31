use std::collections::HashMap;
use std::error::Error;
use std::io::Write;
use std::{fs, str};

use derive_new::new;
use itertools::Itertools;
use serde_derive::Serialize;
use serde_json;
use tera::{self, Context, Tera};

use crate::bcf::report::table_report::create_report_table::get_ann_description;
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
) -> Result<(), Box<dyn Error>> {
    let mut data = HashMap::new();
    let mut gene_data = HashMap::new();
    let mut impact_data = Vec::new();
    let mut clin_sig_data = Vec::new();
    let mut unique_genes = HashMap::new();
    for (sample, path) in sample_calls.iter().sorted() {
        let mut genes = HashMap::new();
        let mut impacts = HashMap::new();
        let mut clin_sigs = HashMap::new();
        let mut ann_indices = HashMap::new();
        let mut bcf_reader = bcf::Reader::from_path(path)?;
        let header = bcf_reader.header().clone();
        let header_records = header.header_records();
        let ann_fields: Vec<_> = get_ann_description(header_records).unwrap();

        for (i, field) in ann_fields.iter().enumerate() {
            ann_indices.insert(field, i);
        }

        for res in bcf_reader.records() {
            let mut record = res?;
            let alleles = record
                .alleles()
                .into_iter()
                .map(|allele| allele.to_owned())
                .collect_vec();
            let alt_alleles = &alleles[1..];
            let ref_allele = alleles[0].to_owned();

            let _af = record.format(b"AF").float()?[0]
                .into_iter()
                .fold(f32::INFINITY, |a, &b| a.min(b));

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

                        let impact = str::from_utf8(
                            fields[*ann_indices.get(&String::from("IMPACT")).expect("No field named IMPACT found. Please only use VEP-annotated VCF-files.")],
                        )?;
                        let clin_sig = str::from_utf8(
                            fields[*ann_indices.get(&String::from("CLIN_SIG")).expect("No field named CLIN_SIG found. Please only use VEP-annotated VCF-files.")],
                        )?;
                        let gene = str::from_utf8(
                            fields[*ann_indices.get(&String::from("Gene")).expect("No field named Gene found. Please only use VEP-annotated VCF-files.")],
                        )?;
                        let dna_alteration = str::from_utf8(
                            fields[*ann_indices.get(&String::from("HGVSc")).expect("No field named HGVSc found. Please only use VEP-annotated VCF-files.")],
                        )?;
                        let protein_alteration = str::from_utf8(
                            fields[*ann_indices.get(&String::from("HGVSp")).expect("No field named HGVSp found. Please only use VEP-annotated VCF-files.")],
                        )?;

                        let gene_rec = unique_genes
                            .entry(gene.to_owned())
                            .or_insert_with(|| Vec::new());
                        gene_rec.push(sample.to_owned());

                        let rec = genes
                            .entry(gene.to_owned())
                            .or_insert_with(|| Record::new(sample.to_owned(), gene.to_owned()));
                        rec.dna_alterations.push(dna_alteration.to_owned());
                        if protein_alteration.len() > 0 {
                            rec.protein_alterations.push(protein_alteration.to_owned());
                        }
                        rec.variants.push(variant.to_owned());

                        let imp_rec = impacts
                            .entry(gene.to_owned())
                            .or_insert_with(|| Impact::new(sample.to_owned(), gene.to_owned()));
                        imp_rec.impact.push(impact.to_owned());
                        let clin_rec = clin_sigs
                            .entry(gene.to_owned())
                            .or_insert_with(|| ClinSig::new(sample.to_owned(), gene.to_owned()));
                        let sigs: Vec<_> = clin_sig.split('&').collect();
                        for s in sigs {
                            clin_rec.clin_sig.push(s.to_owned());
                        }
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

            // data for impact
            let impact = impacts.get(gene).unwrap();
            let final_impacts = Vec::<FinalImpact>::from(impact);
            impact_data.push(final_impacts);

            // data for clin_sig
            let cls = clin_sigs.get(gene).unwrap();
            let final_clin_sig = Vec::<FinalClinSig>::from(cls);
            clin_sig_data.push(final_clin_sig);
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
    let mut sort_genes = HashMap::new();

    // remove duplicate samples and calculate order for oncoprint
    for (gene, mut samples) in unique_genes {
        samples.sort();
        samples.dedup();
        sort_genes.insert(gene.to_owned(), samples.len());
    }

    let impact_data: Vec<_> = impact_data.iter().flatten().collect();
    let clin_sig_data: Vec<_> = clin_sig_data.iter().flatten().collect();

    let page_size = 100;

    let order = json!(sort_genes);

    let mut v = Vec::from_iter(sort_genes);
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

        let impact_page_data: Vec<_> = impact_data
            .iter()
            .filter(|entry| sorted_genes.contains(&&entry.gene))
            .collect();

        let clin_sig_page_data: Vec<_> = clin_sig_data
            .iter()
            .filter(|entry| sorted_genes.contains(&&entry.gene))
            .collect();

        let mut vl_specs: Value = serde_json::from_str(include_str!("report_specs.json")).unwrap();
        let values = json!({"main": page_data , "impact": impact_page_data , "clin_sig": clin_sig_page_data});

        vl_specs["datasets"] = values;

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
        context.insert("order", &serde_json::to_string(&order)?);
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
}

#[derive(new, Debug)]
struct Impact {
    sample: String,
    gene: String,
    #[new(default)]
    impact: Vec<String>,
}

#[derive(Serialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct FinalImpact {
    sample: String,
    gene: String,
    count: u32,
    impact: String,
    sort: usize,
}

#[derive(new, Debug)]
struct ClinSig {
    sample: String,
    gene: String,
    #[new(default)]
    clin_sig: Vec<String>,
}

#[derive(Serialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct FinalClinSig {
    sample: String,
    gene: String,
    count: u32,
    clin_sig: String,
    sort: usize,
}

// Record for the first stage of the report
#[derive(Serialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct FinalRecord {
    sample: String,
    gene: String,
    dna_alterations: String,
    protein_alterations: String,
    variants: String,
    sort: usize,
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
        gene_vec
    }
}

impl From<&Record> for FinalRecord {
    fn from(record: &Record) -> Self {
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
            sort: 0,
        }
    }
}

impl From<&Impact> for Vec<FinalImpact> {
    fn from(impact: &Impact) -> Self {
        let mut imp_map = HashMap::new();
        for i in &impact.impact {
            let rec = imp_map.entry(i.to_owned()).or_insert_with(|| 0);
            *rec += 1;
        }

        let mut res = Vec::new();

        for (imp, count) in imp_map {
            let record = FinalImpact {
                sample: impact.sample.clone(),
                gene: impact.gene.clone(),
                count: count,
                impact: imp,
                sort: 0,
            };

            res.push(record);
        }

        res
    }
}

impl From<&ClinSig> for Vec<FinalClinSig> {
    fn from(clin: &ClinSig) -> Self {
        let mut clin_map = HashMap::new();
        for i in &clin.clin_sig {
            let rec = clin_map.entry(i.to_owned()).or_insert_with(|| 0);
            *rec += 1;
        }

        let mut res = Vec::new();

        for (c, count) in clin_map {
            let record = FinalClinSig {
                sample: clin.sample.clone(),
                gene: clin.gene.clone(),
                count: count,
                clin_sig: c,
                sort: 0,
            };

            res.push(record);
        }

        res
    }
}
