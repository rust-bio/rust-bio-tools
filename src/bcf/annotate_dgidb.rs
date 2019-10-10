use itertools::Itertools;
use reqwest;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::str;

#[derive(Serialize, Deserialize, Debug)]
struct Dgidb {
    matchedTerms: Vec<MatchedTerm>,
}

#[derive(Serialize, Deserialize, Debug)]
struct MatchedTerm {
    geneName: String,
    interactions: Vec<Interaction>,
}

#[derive(Serialize, Deserialize, Debug)]
struct Interaction {
    drugName: String,
}

pub fn annotate_dgidb(vcf_path: &str) -> Result<(), Box<dyn Error>> {
    let gene_drug_interactions = request_interaction_drugs(vcf_path)?;
    modify_vcf_entries(vcf_path, gene_drug_interactions)
}

fn request_interaction_drugs(
    vcf_path: &str,
) -> Result<HashMap<String, Vec<String>>, Box<dyn Error>> {
    let mut genes = collect_genes(vcf_path)?;
    let mut url = "http://dgidb.org/api/v2/interactions.json?genes=".to_string();
    url.push_str(genes.drain().join(",").as_str());
    let res: Dgidb = reqwest::get(&url)?.json()?;

    let mut gene_drug_interactions: HashMap<String, Vec<String>> = HashMap::new();
    for term in res.matchedTerms {
        if !term.interactions.is_empty() {
            gene_drug_interactions.insert(
                term.geneName,
                term.interactions
                    .iter()
                    .map(|interaction| interaction.drugName.clone())
                    .collect(),
            );
        }
    }
    Ok(gene_drug_interactions)
}
//TODO Remove split by ',' when updating to latest htslib release
//TODO Only use first transcripts entry and split transcripts by ','
fn collect_genes(vcf_path: &str) -> Result<HashSet<String>, Box<dyn Error>> {
    let mut total_genes = HashSet::new();
    let mut reader = bcf::Reader::from_path(vcf_path)?;
    //let all_genes = HashSet::new();
    for result in reader.records() {
        let mut rec = result?;
        let annotation = rec.info("ANN".as_bytes()).string()?;
        match annotation {
            Some(transcripts) => transcripts.iter().for_each(|transcript| {
                let gene = str::from_utf8(transcript)
                    .unwrap()
                    .split("|")
                    .nth(3)
                    .unwrap()
                    .to_string();
                total_genes.insert(gene);
            }),
            None => {}
        };
    }
    Ok(total_genes)
}

//TODO Only use first transcripts entry and split transcripts by ','
fn modify_vcf_entries(
    vcf_path: &str,
    gene_drug_interactions: HashMap<String, Vec<String>>,
) -> Result<(), Box<dyn Error>> {
    let mut reader = bcf::Reader::from_path(vcf_path)?;
    let mut writer =
        bcf::Writer::from_stdout(&bcf::Header::from_template(reader.header()), true, true)?;

    for result in reader.records() {
        let mut rec = result?;
        let annotation = rec.info("ANN".as_bytes()).string()?;
        match annotation {
            Some(transcripts) => {
                dbg!(transcripts.len());
                let genes: Vec<String> = transcripts.iter().map(|transcript| {
                    str::from_utf8(transcript)
                        .unwrap()
                        .split("|")
                        .nth(3)
                        .unwrap()
                        .to_string()
                }).collect();
                dbg!(&genes);
            }
            None => {}
        }
        //TODO Reenable after implementation is completed
        //writer.write(&rec)?;
    }
    Ok(())
}
