use rust_htslib::bcf;
use std::error::Error;
use std::collections::{HashSet, HashMap};
use rust_htslib::bcf::Read;
use std::str;
use itertools::Itertools;
use reqwest;
use rustc_serialize::json::ToJson;

pub fn annotate_dgidb(vcf_path: &str) -> Result<(), Box<dyn Error>> {
    request_interaction_drugs(vcf_path)?;
    modify_vcf_entries(vcf_path)
}

fn request_interaction_drugs(vcf_path: &str) -> Result<(), Box<dyn Error>> {
    let mut genes = collect_genes(vcf_path)?;
    let mut url = "http://dgidb.org/api/v2/interactions.json?genes=".to_string();
    url.push_str( genes.drain().join(",").as_str());
    let response = reqwest::get(url.as_str())?.text()?.to_json();
    dbg!(&response);
    Ok(())

}

fn collect_genes(vcf_path: &str) -> Result<HashSet<String>, Box<dyn Error>> {
    let mut total_genes = HashSet::new();
    let mut reader = bcf::Reader::from_path(vcf_path)?;
    //let all_genes = HashSet::new();
    for result in reader.records() {
        let mut rec = result?;
        let annotation = rec.info("ANN".as_bytes()).string()?;
        match annotation {
            Some(transcripts) => {
                transcripts.iter().for_each(|transcript| {
                    let gene = str::from_utf8(transcript).unwrap().split("|").nth(3).unwrap().to_string();
                    total_genes.insert(gene);
                })
            },
            None => {}
        };
    }
    Ok(total_genes)
}

fn modify_vcf_entries(vcf_path: &str) -> Result<(), Box<dyn Error>> {
    let mut reader = bcf::Reader::from_path(vcf_path)?;
    let mut writer = bcf::Writer::from_stdout(&bcf::Header::from_template(reader.header()), true, true)?;


    for result in reader.records() {
        let rec = result?;
        writer.write(&rec)?;
    }
    Ok(())
}

