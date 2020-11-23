use itertools::Itertools;
use regex::Regex;
use rust_htslib::bcf;
use rust_htslib::bcf::{Format, Read};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::str;

#[derive(Serialize, Deserialize, Debug)]
struct Dgidb {
    #[serde(rename = "matchedTerms")]
    matched_terms: Vec<MatchedTerm>,
}

#[derive(Serialize, Deserialize, Debug)]
struct MatchedTerm {
    #[serde(rename = "geneName")]
    gene_name: String,
    interactions: Vec<Interaction>,
}

#[derive(Serialize, Deserialize, Debug)]
struct Interaction {
    #[serde(rename = "drugName")]
    drug_name: String,
    #[serde(rename = "interactionTypes")]
    interaction_types: Vec<String>,
}

pub fn annotate_dgidb(
    vcf_path: &str,
    api_path: String,
    field_name: &str,
    datasources: Option<Vec<&str>>,
    genes_per_request: usize,
) -> Result<(), Box<dyn Error>> {
    let gene_drug_interactions =
        request_interaction_drugs(vcf_path, api_path, datasources, genes_per_request)?;
    modify_vcf_entries(vcf_path, gene_drug_interactions, field_name)
}

type Interactions = HashMap<String, Vec<(String, Vec<String>)>>;

fn request_interaction_drugs(
    vcf_path: &str,
    api_path: String,
    datasources_opt: Option<Vec<&str>>,
    genes_per_request: usize,
) -> Result<Option<Interactions>, Box<dyn Error>> {
    let mut genes = collect_genes(vcf_path)?;
    let datasources = if let Some(entries) = datasources_opt {
        let mut b = String::from("&interaction_sources=");
        b.push_str(entries.join(",").as_str());
        b
    } else {
        String::new()
    };
    if genes.is_empty() {
        return Ok(None);
    }
    let mut gene_drug_interactions: HashMap<String, Vec<(String, Vec<String>)>> = HashMap::new();
    for gene_slice in genes.drain().collect_vec().chunks(genes_per_request) {
        let mut slice_api_path = api_path.clone();
        slice_api_path.push_str(gene_slice.join(",").as_str());
        slice_api_path.push_str(datasources.as_str());
        let res: Dgidb = reqwest::get(&slice_api_path)?.json()?;

        for term in res.matched_terms {
            if !term.interactions.is_empty() {
                gene_drug_interactions.insert(
                    term.gene_name,
                    term.interactions
                        .iter()
                        .map(|interaction| {
                            (
                                interaction.drug_name.clone(),
                                interaction.interaction_types.clone(),
                            )
                        })
                        .collect(),
                );
            }
        }
    }
    Ok(Some(gene_drug_interactions))
}

fn collect_genes(vcf_path: &str) -> Result<HashSet<String>, Box<dyn Error>> {
    let mut total_genes = HashSet::new();
    let mut reader = bcf::Reader::from_path(vcf_path)?;
    for result in reader.records() {
        let mut rec = result?;
        let genes_opt = extract_genes(&mut rec)?;
        if let Some(genes) = genes_opt {
            for gene in genes {
                total_genes.insert(gene);
            }
        }
    }
    Ok(total_genes)
}

fn extract_genes(
    rec: &mut bcf::Record,
) -> Result<Option<impl Iterator<Item = String> + '_>, Box<dyn Error>> {
    let annotation = rec.info(b"ANN").string()?;
    match annotation {
        Some(transcripts) => Ok(Some(transcripts.clone().into_iter().map(|transcript| {
            str::from_utf8(transcript.split(|c| *c == b'|').nth(3).unwrap())
                .unwrap()
                .to_owned()
        }))),
        None => Ok(None),
    }
}

fn modify_vcf_entries(
    vcf_path: &str,
    gene_drug_interactions_opt: Option<Interactions>,
    field_name: &str,
) -> Result<(), Box<dyn Error>> {
    let mut reader = bcf::Reader::from_path(vcf_path)?;
    let mut header = bcf::header::Header::from_template(reader.header());
    header.push_record(format!("##INFO=<ID={},Number=.,Type=String,Description=\"Combination of gene, drug, interaction types extracted from dgiDB. Each combination is pipe-seperated annotated as GENE|DRUG|TYPE\">", field_name).as_bytes());
    let mut writer = bcf::Writer::from_stdout(&header, true, Format::BCF)?;
    match gene_drug_interactions_opt {
        None => {
            for result in reader.records() {
                let mut rec = result?;
                writer.translate(&mut rec);
                writer.write(&rec)?;
            }
        }
        Some(gene_drug_interactions) => {
            for result in reader.records() {
                let mut rec = result?;
                writer.translate(&mut rec);
                let genes = extract_genes(&mut rec)?.map(|genes| genes.collect_vec());
                if let Some(mut genes) = genes {
                    genes.sort();
                    genes.dedup();
                    let field_entries = build_dgidb_field(&gene_drug_interactions, genes)?;
                    let field_entries: Vec<&[u8]> =
                        field_entries.iter().map(|v| v.as_slice()).collect();
                    rec.push_info_string(field_name.as_bytes(), &field_entries[..])?;
                }
                writer.write(&rec)?;
            }
        }
    }
    Ok(())
}

fn build_dgidb_field(
    gene_drug_interactions: &HashMap<String, Vec<(String, Vec<String>)>>,
    genes: Vec<String>,
) -> Result<Vec<Vec<u8>>, Box<dyn Error>> {
    let mut field_entries: Vec<Vec<u8>> = Vec::new();
    let re = Regex::new(r"\s\(\w+\)").unwrap();
    for gene in genes.iter() {
        match gene_drug_interactions.get(gene) {
            Some(drug_interactions) => {
                drug_interactions
                    .iter()
                    .for_each(|(drug, interaction_types)| {
                        if !interaction_types.is_empty() {
                            interaction_types.iter().for_each(|interaction_type| {
                                field_entries.push(
                                    format!(
                                        "{g}|{d}|{t}",
                                        g = gene,
                                        d = re.replace(drug, ""),
                                        t = interaction_type
                                    )
                                    .as_bytes()
                                    .to_vec(),
                                )
                            })
                        } else {
                            field_entries.push(
                                format!("{g}|{d}|.", g = gene, d = re.replace(drug, ""))
                                    .as_bytes()
                                    .to_vec(),
                            )
                        }
                    });
            }
            None => field_entries.push(format!("{g}|.|.", g = gene).as_bytes().to_vec()),
        }
    }
    Ok(field_entries)
}
