use itertools::Itertools;
use reqwest;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
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
}

pub fn annotate_dgidb(
    vcf_path: &str,
    api_path: String,
    field_name: &str,
) -> Result<(), Box<dyn Error>> {
    let gene_drug_interactions = request_interaction_drugs(vcf_path, api_path)?;
    modify_vcf_entries(vcf_path, gene_drug_interactions, field_name)
}

fn request_interaction_drugs(
    vcf_path: &str,
    mut api_path: String,
) -> Result<HashMap<String, Vec<String>>, Box<dyn Error>> {
    let mut genes = collect_genes(vcf_path)?;
    api_path.push_str(genes.drain().join(",").as_str());
    let res: Dgidb = reqwest::get(&api_path)?.json()?;

    let mut gene_drug_interactions: HashMap<String, Vec<String>> = HashMap::new();
    for term in res.matched_terms {
        if !term.interactions.is_empty() {
            gene_drug_interactions.insert(
                term.gene_name,
                term.interactions
                    .iter()
                    .map(|interaction| interaction.drug_name.clone())
                    .collect(),
            );
        }
    }
    Ok(gene_drug_interactions)
}
//TODO Remove split by ',' when updating to latest htslib release
fn collect_genes(vcf_path: &str) -> Result<HashSet<String>, Box<dyn Error>> {
    let mut total_genes = HashSet::new();
    let mut reader = bcf::Reader::from_path(vcf_path)?;
    //let all_genes = HashSet::new();
    for result in reader.records() {
        let mut rec = result?;
        let annotation = rec.info("ANN".as_bytes()).string()?;
        match annotation {
            Some(transcripts) => {
                str::from_utf8(transcripts[0])
                    .unwrap()
                    .split(",")
                    .for_each(|transcript| {
                        let gene = transcript.split("|").nth(3).unwrap().to_string();
                        //TODO This should be the code for the latest htslib implementation (Needs to be tested)
                        /*
                        transcripts.iter().for_each(|transcript| {
                        let gene = str::from_utf8(transcript)
                            .unwrap()
                            .split("|")
                            .nth(3)
                            .unwrap()
                            .to_string();
                         */
                        total_genes.insert(gene);
                    })
            }
            None => {}
        };
    }
    Ok(total_genes)
}

//TODO Remove split by ',' when updating to latest htslib release
fn modify_vcf_entries(
    vcf_path: &str,
    gene_drug_interactions: HashMap<String, Vec<String>>,
    field_name: &str,
) -> Result<(), Box<dyn Error>> {
    let mut reader = bcf::Reader::from_path(vcf_path)?;
    let mut header = bcf::header::Header::from_template(reader.header());
    header.push_record(format!("##INFO=<ID={},Number=.,Type=String,Description=\"Interacting drugs for each gene extracted from dgiDB. Multiple drugs for one gene are pipe-seperated.\">", field_name).as_bytes());
    let mut writer = bcf::Writer::from_stdout(&header, true, true)?;

    for result in reader.records() {
        let mut rec = result?;
        let annotation = rec.info("ANN".as_bytes()).string()?;
        match annotation {
            Some(transcripts) => {
                let genes = str::from_utf8(transcripts[0])
                    .unwrap()
                    .split(',')
                    .map(|transcript| transcript.split("|").nth(3).unwrap().to_string())
                    .collect();
                //TODO This should be the code for the latest htslib implementation (Needs to be tested)
                /*let genes = transcripts.iter().map(|transcript| {
                    str::from_utf8(transcript)
                        .unwrap()
                        .split("|")
                        .nth(3)
                        .unwrap()
                        .to_string()
                }).collect();
                */
                let field_entries = build_dgidb_field(&gene_drug_interactions, genes)?;
                let field_entries: Vec<&[u8]> =
                    field_entries.iter().map(|v| v.as_slice()).collect();
                writer.translate(&mut rec); //That's bullshit
                rec.push_info_string(field_name.as_bytes(), &field_entries[..])?;
                writer.write(&rec)?;
            }
            None => {}
        };
    }
    Ok(())
}

fn build_dgidb_field(
    gene_drug_interactions: &HashMap<String, Vec<String>>,
    genes: Vec<String>,
) -> Result<Vec<Vec<u8>>, Box<dyn Error>> {
    let mut field_entries: Vec<Vec<u8>> = Vec::new();
    for gene in genes.iter() {
        let drugs_opt = gene_drug_interactions.get(gene);
        match drugs_opt {
            Some(drugs) => field_entries.push(drugs.join("|").as_bytes().to_vec()),
            None => field_entries.push([b'.'].to_vec()),
        }
    }
    Ok(field_entries)
}
