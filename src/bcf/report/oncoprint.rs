use std::collections::HashMap;
use std::error::Error;
use std::io::Write;
use std::{fs, str};

use derive_new::new;
use itertools::Itertools;
use serde_derive::Serialize;
use tera::{self, Context, Tera};

use crate::bcf::report::table_report::create_report_table::get_ann_description;
use chrono::{DateTime, Local};
use jsonm::packer::{PackOptions, Packer};
use rust_htslib::bcf::{self, Read};
use serde_json::{json, Value};
use std::fs::File;
use std::path::Path;
use std::str::FromStr;

pub fn oncoprint(
    sample_calls: &HashMap<String, String>,
    output_path: &str,
    max_cells: u32,
) -> Result<(), Box<dyn Error>> {
    let mut data = HashMap::new();
    let mut gene_data = HashMap::new();
    let mut impact_data = Vec::new();
    let mut gene_impact_data = HashMap::new();
    let mut consequence_data = Vec::new();
    let mut gene_consequence_data = HashMap::new();
    let mut clin_sig_data = Vec::new();
    let mut gene_clin_sig_data = HashMap::new();
    let mut af_data = Vec::new();
    let mut gene_af_data = HashMap::new();
    let mut unique_genes = HashMap::new();
    for (sample, path) in sample_calls.iter().sorted() {
        let mut genes = HashMap::new();
        let mut impacts = HashMap::new();
        let mut gene_impacts = HashMap::new();
        let mut consequences = HashMap::new();
        let mut gene_consequences = HashMap::new();
        let mut clin_sigs = HashMap::new();
        let mut gene_clin_sigs = HashMap::new();
        let mut ann_indices = HashMap::new();
        let mut bcf_reader = bcf::Reader::from_path(path)?;
        let header = bcf_reader.header().clone();
        let mut sample_names = Vec::new();
        for s in header.samples() {
            sample_names.push(String::from_utf8(s.to_owned()).unwrap());
        }
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

            let allel_frequencies = record.format(b"AF").float()?[0].to_vec();

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

                        let get_field = |field: &str| {
                            str::from_utf8(
                                fields[*ann_indices.get(&field.to_owned()).unwrap_or_else(|| panic!("No field named {} found. Please only use VEP-annotated VCF-files.", field))],
                            )
                        };

                        let mut impact = get_field("IMPACT")?;
                        let clin_sig = get_field("CLIN_SIG")?;
                        let gene = get_field("SYMBOL")?;
                        let dna_alteration = get_field("HGVSc")?;
                        let protein_alteration = get_field("HGVSp")?;
                        let consequence = get_field("Consequence")?;

                        let gene_rec = unique_genes.entry(gene.to_owned()).or_insert_with(Vec::new);
                        gene_rec.push(sample.to_owned());

                        let rec = genes
                            .entry(gene.to_owned())
                            .or_insert_with(|| Record::new(sample.to_owned(), gene.to_owned()));
                        rec.dna_alterations.push(dna_alteration.to_owned());
                        if !protein_alteration.is_empty() {
                            rec.protein_alterations.push(protein_alteration.to_owned());
                        }
                        rec.variants.push(variant.to_owned());

                        let imp_rec = impacts.entry(gene.to_owned()).or_insert_with(Vec::new);
                        if impact == "" {
                            impact = "unknown";
                        }
                        imp_rec.push(BarPlotRecord::new(gene.to_owned(), impact.to_owned()));

                        let alt = if protein_alteration.is_empty() {
                            dna_alteration
                        } else {
                            protein_alteration
                        };
                        let gene_imp_rec =
                            gene_impacts.entry(gene.to_owned()).or_insert_with(Vec::new);
                        gene_imp_rec.push(BarPlotRecord::new(alt.to_owned(), impact.to_owned()));

                        let split_consequences: Vec<_> = consequence.split('&').collect();

                        let cons_rec = consequences.entry(gene.to_owned()).or_insert_with(Vec::new);

                        let gene_cons_rec = gene_consequences
                            .entry(gene.to_owned())
                            .or_insert_with(Vec::new);

                        for mut c in split_consequences {
                            if c == "" {
                                c = "unknown";
                            }
                            cons_rec.push(BarPlotRecord::new(gene.to_owned(), c.to_owned()));
                            gene_cons_rec.push(BarPlotRecord::new(alt.to_owned(), c.to_owned()));
                        }

                        let gene_clin_sig_rec = gene_clin_sigs
                            .entry(gene.to_owned())
                            .or_insert_with(Vec::new);

                        let clin_rec = clin_sigs.entry(gene.to_owned()).or_insert_with(Vec::new);
                        let sigs: Vec<_> = clin_sig.split('&').collect();
                        for mut s in sigs {
                            if s == "" {
                                s = "unknown";
                            }
                            s = s.trim_matches('_');
                            clin_rec.push(BarPlotRecord::new(gene.to_owned(), s.to_owned()));
                            gene_clin_sig_rec
                                .push(BarPlotRecord::new(alt.to_owned(), s.to_owned()));
                        }

                        for (i, name) in sample_names.iter().enumerate() {
                            let af = AlleleFrequency {
                                sample: sample.to_owned() + ":" + name,
                                key: gene.to_owned(),
                                allel_frequency: allel_frequencies[i],
                            };

                            af_data.push(af);

                            let gene_af = AlleleFrequency {
                                sample: sample.to_owned() + ":" + name,
                                key: alt.to_owned(),
                                allel_frequency: allel_frequencies[i],
                            };

                            let f = gene_af_data.entry(gene.to_owned()).or_insert_with(Vec::new);
                            f.push(gene_af);
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

            // data for first stage impact
            let impact = impacts.get(gene).unwrap();
            let final_impacts = make_final_bar_plot_records(impact);
            impact_data.push(final_impacts);

            let consequence = consequences.get(gene).unwrap();
            let final_consequences = make_final_bar_plot_records(consequence);
            consequence_data.push(final_consequences);

            // data for first stage clin_sig
            let cls = clin_sigs.get(gene).unwrap();
            let final_clin_sig = make_final_bar_plot_records(cls);
            clin_sig_data.push(final_clin_sig);

            // data for second stage
            let gene_entry = gene_data.entry(gene.to_owned()).or_insert_with(&Vec::new);
            let gene_records: Vec<SecondStageRecord> = Vec::<SecondStageRecord>::from(record);
            for rec in gene_records {
                gene_entry.push(rec);
            }

            // data for second stage impact
            let gene_impact = gene_impacts.get(gene).unwrap();
            let final_gene_impacts = make_final_bar_plot_records(gene_impact);
            let e = gene_impact_data
                .entry(gene.to_owned())
                .or_insert_with(Vec::new);
            e.push(final_gene_impacts);

            // data for second stage consequence
            let gene_consequence = gene_consequences.get(gene).unwrap();
            let final_gene_consequences = make_final_bar_plot_records(gene_consequence);
            let g = gene_consequence_data
                .entry(gene.to_owned())
                .or_insert_with(Vec::new);
            g.push(final_gene_consequences);

            // data for second stage clin_sig
            let gene_clin_sig = gene_clin_sigs.get(gene).unwrap();
            let final_gene_clin_sigs = make_final_bar_plot_records(gene_clin_sig);
            let f = gene_clin_sig_data
                .entry(gene.to_owned())
                .or_insert_with(Vec::new);
            f.push(final_gene_clin_sigs);
        }
    }

    let gene_path = output_path.to_owned() + "/genes/";
    fs::create_dir(Path::new(&gene_path))?;
    let mut gene_templates = Tera::default();
    gene_templates.add_raw_template("genes.html.tera", include_str!("genes.html.tera"))?;

    let gene_specs: Value = serde_json::from_str(include_str!("gene_specs.json")).unwrap();

    let page_size = max_cells as usize / sample_calls.len();

    // create html for second stage
    for (gene, data) in gene_data {
        let gene_data: Vec<_> = data.iter().sorted().collect();
        let mut alterations = HashMap::new();
        for rec in &gene_data {
            let entry = alterations.entry(&rec.alteration).or_insert_with(Vec::new);
            entry.push(rec.sample.to_owned());
        }
        let mut sort_alterations = HashMap::new();
        for (alteration, mut samples) in alterations {
            samples.sort();
            samples.dedup();
            sort_alterations.insert(alteration, samples.len());
        }

        let mut specs = gene_specs.clone();
        let impact_data = gene_impact_data.get(&gene).unwrap();
        let final_impact: Vec<_> = impact_data.iter().flatten().sorted().collect();
        let consequence_data = gene_consequence_data.get(&gene).unwrap();
        let final_consequence: Vec<_> = consequence_data.iter().flatten().sorted().collect();
        let clin_sig_data = gene_clin_sig_data.get(&gene).unwrap();
        let final_clin_sig: Vec<_> = clin_sig_data.iter().flatten().sorted().collect();
        let allel_frequency_data = gene_af_data.get(&gene).unwrap();

        let sorted_impacts = order_by_impact(final_impact.clone());
        let sorted_clin_sigs = order_by_clin_sig(final_clin_sig.clone());
        let mut order = Vec::new();
        for (alt, sample_count) in sort_alterations {
            let impact_order = sorted_impacts.get(alt).unwrap();
            let clin_sig_order = sorted_clin_sigs.get(alt).unwrap();
            order.push((alt.to_owned(), sample_count, impact_order, clin_sig_order))
        }

        order.sort_by(|(g1, c1, i1, cs1), (g2, c2, i2, cs2)| {
            c2.cmp(c1) // first order by different sample occurrences
                .then(i2.cmp(i1)) // then by impact
                .then(cs2.cmp(cs1)) // then by clin_sig
                .then(g2.cmp(g1)) // lastly by alteration name for consistency
        });

        let ordered_alts: Vec<_> = order.iter().map(|(x, _, _, _)| x).collect();

        let pages = if order.len() % page_size == 0 {
            (order.len() / page_size) - 1
        } else {
            order.len() / page_size
        };

        for i in 0..pages + 1 {
            let current_alterations = if i != pages {
                &order[(i * page_size)..((i + 1) * page_size)] // get genes for current page
            } else {
                &order[(i * page_size)..] // get genes for last page
            };

            if !current_alterations.is_empty() {
                let mut sorted_alterations = Vec::new();
                for (g, _, _, _) in current_alterations {
                    sorted_alterations.push(g);
                }

                let page = i + 1;

                let page_data: Vec<_> = gene_data
                    .iter()
                    .filter(|entry| sorted_alterations.contains(&&&entry.alteration))
                    .sorted()
                    .collect();

                let impact_page_data: Vec<_> = final_impact
                    .iter()
                    .filter(|entry| sorted_alterations.contains(&&&entry.record.key))
                    .sorted()
                    .collect();

                let consequence_page_data: Vec<_> = final_consequence
                    .iter()
                    .filter(|entry| sorted_alterations.contains(&&&entry.record.key))
                    .sorted()
                    .collect();

                let clin_sig_page_data: Vec<_> = final_clin_sig
                    .iter()
                    .filter(|entry| sorted_alterations.contains(&&&entry.record.key))
                    .sorted()
                    .collect();

                let af_page_data: Vec<_> = allel_frequency_data
                    .iter()
                    .filter(|entry| sorted_alterations.contains(&&&entry.key))
                    .collect();

                let order: Vec<_> = ordered_alts
                    .iter()
                    .filter(|gene| sorted_alterations.contains(gene))
                    .collect();

                let values = json!({ "main": page_data, "impact": impact_page_data, "consequence": consequence_page_data, "clin_sig": clin_sig_page_data, "allel_frequency": af_page_data});
                specs["datasets"] = values;
                let mut packer = Packer::new();
                let options = PackOptions::new();
                let packed_gene_specs = packer.pack(&specs, &options).unwrap();
                let mut context = Context::new();
                context.insert("genespecs", &serde_json::to_string(&packed_gene_specs)?);
                context.insert("gene", &gene);
                context.insert("samples", &sample_calls.len());
                context.insert("current_page", &page);
                context.insert("pages", &(pages + 1));
                context.insert("order", &serde_json::to_string(&json!(order))?);
                let local: DateTime<Local> = Local::now();
                context.insert("time", &local.format("%a %b %e %T %Y").to_string());
                context.insert("version", &env!("CARGO_PKG_VERSION"));
                let html = gene_templates.render("genes.html.tera", &context)?;
                let filepath = gene_path.clone() + &gene + &page.to_string() + ".html";
                let mut file = File::create(filepath)?;
                file.write_all(html.as_bytes())?;
            }
        }
    }

    // only keep recurrent entries
    let data: Vec<_> = data
        .values()
        .filter(|entry| !entry.is_empty())
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
    let consequence_data: Vec<_> = consequence_data.iter().flatten().collect();
    let clin_sig_data: Vec<_> = clin_sig_data.iter().flatten().collect();

    let sorted_impacts = order_by_impact(impact_data.clone());
    let sorted_clin_sigs = order_by_clin_sig(clin_sig_data.clone());
    let mut v = Vec::new();
    for (gene, sample_count) in sort_genes {
        let impact_order = sorted_impacts.get(&gene).unwrap();
        let clin_sig_order = sorted_clin_sigs.get(&gene).unwrap();
        v.push((gene, sample_count, impact_order, clin_sig_order))
    }

    v.sort_by(|(g1, c1, i1, cs1), (g2, c2, i2, cs2)| {
        c2.cmp(c1) // first order by different sample occurrences
            .then(i2.cmp(i1)) // then by impact
            .then(cs2.cmp(cs1)) // then by clin_sig
            .then(g2.cmp(g1)) // lastly by gene name for consistency
    });
    let ordered_genes: Vec<_> = v.iter().map(|(x, _, _, _)| x).collect();

    let pages = if v.len() % page_size == 0 {
        (v.len() / page_size) - 1
    } else {
        v.len() / page_size
    };

    let index_path = output_path.to_owned() + "/indexes/";
    fs::create_dir(Path::new(&index_path))?;

    for i in 0..pages + 1 {
        let current_genes = if i != pages {
            &v[(i * page_size)..((i + 1) * page_size)] // get genes for current page
        } else {
            &v[(i * page_size)..] // get genes for last page
        };

        if !current_genes.is_empty() {
            let mut sorted_genes = Vec::new();
            for (g, _, _, _) in current_genes {
                sorted_genes.push(g);
            }

            let page = i + 1;

            let page_data: Vec<_> = data
                .iter()
                .filter(|entry| sorted_genes.contains(&&entry.gene))
                .sorted()
                .collect();

            let impact_page_data: Vec<_> = impact_data
                .iter()
                .filter(|entry| sorted_genes.contains(&&entry.record.key))
                .sorted()
                .collect();

            let consequence_page_data: Vec<_> = consequence_data
                .iter()
                .filter(|entry| sorted_genes.contains(&&entry.record.key))
                .sorted()
                .collect();

            let clin_sig_page_data: Vec<_> = clin_sig_data
                .iter()
                .filter(|entry| sorted_genes.contains(&&entry.record.key))
                .sorted()
                .collect();

            let af_page_data: Vec<_> = af_data
                .iter()
                .filter(|entry| sorted_genes.contains(&&entry.key))
                .collect();

            let order: Vec<_> = ordered_genes
                .iter()
                .filter(|gene| sorted_genes.contains(gene))
                .collect();

            let mut vl_specs: Value =
                serde_json::from_str(include_str!("report_specs.json")).unwrap();
            let values = json!({"main": page_data , "impact": impact_page_data, "consequence": consequence_page_data , "clin_sig": clin_sig_page_data, "allel_frequency": af_page_data});

            vl_specs["datasets"] = values;

            let mut packer = Packer::new();
            let options = PackOptions::new();
            let packed_specs = packer.pack(&vl_specs, &options).unwrap();
            let mut templates = Tera::default();
            templates.add_raw_template("report.html.tera", include_str!("report.html.tera"))?;
            let mut context = Context::new();
            let data = serde_json::to_string(&packed_specs)?;
            context.insert("oncoprint", &data);
            context.insert("current_page", &page);
            context.insert("pages", &(pages + 1));
            context.insert("order", &serde_json::to_string(&json!(order))?);
            context.insert("samples", &sample_calls.len());
            let local: DateTime<Local> = Local::now();
            context.insert("time", &local.format("%a %b %e %T %Y").to_string());
            context.insert("version", &env!("CARGO_PKG_VERSION"));

            let html = templates.render("report.html.tera", &context)?;

            let index = index_path.to_owned() + "/index" + &page.to_string() + ".html";
            let mut file = File::create(index)?;
            file.write_all(html.as_bytes())?;
        }
    }

    Ok(())
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

#[derive(Serialize, Debug, PartialEq, PartialOrd, Clone)]
struct AlleleFrequency {
    sample: String,
    key: String,
    allel_frequency: f32,
}

#[derive(new, Serialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct BarPlotRecord {
    key: String,
    value: String,
}

#[derive(Serialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct Counter {
    count: u32,
    #[serde(flatten)]
    record: BarPlotRecord,
}

// Record for the first stage of the report
#[derive(Serialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct FinalRecord {
    sample: String,
    gene: String,
    variants: String,
}

// Record for the second stage of the report
#[derive(Serialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct SecondStageRecord {
    sample: String,
    alteration: String,
    variants: String,
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug)]
enum Impact {
    Unknown,
    Low,
    Modifier,
    Moderate,
    High,
}

impl FromStr for Impact {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let impact;
        match s {
            "HIGH" => impact = Impact::High,
            "MODERATE" => impact = Impact::Moderate,
            "MODIFIER" => impact = Impact::Modifier,
            "LOW" => impact = Impact::Low,
            _ => impact = Impact::Unknown,
        }
        Ok(impact)
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug)]
enum ClinSig {
    Unknown,
    NotProvided,
    Other,
    Benign,
    BenignLikelyBenign,
    LikelyBenign,
    Protective,
    UncertainSignificance,
    Association,
    Affects,
    DrugResponse,
    RiskFactor,
    LikelyPathogenic,
    LikelyPathogenicPathogenic,
    Pathogenic,
}

impl FromStr for ClinSig {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let clin_sig;
        match s {
            "pathogenic" => clin_sig = ClinSig::Pathogenic,
            "likely_pathogenic/pathogenic"=> clin_sig = ClinSig::LikelyPathogenicPathogenic,
            "likely_pathogenic" => clin_sig = ClinSig::LikelyPathogenic,
            "risk_factor" => clin_sig = ClinSig::RiskFactor,
            "drug_response" => clin_sig = ClinSig::DrugResponse,
            "affects" => clin_sig = ClinSig::Affects,
            "association" => clin_sig = ClinSig::Association,
            "uncertain_significance" => clin_sig = ClinSig::UncertainSignificance,
            "protective" => clin_sig = ClinSig::Protective,
            "likely_benign" => clin_sig = ClinSig::LikelyBenign,
            "benign/likely_benign" => clin_sig = ClinSig::BenignLikelyBenign,
            "benign" => clin_sig = ClinSig::Benign,
            "other" => clin_sig = ClinSig::Other,
            "not_provided" => clin_sig = ClinSig::NotProvided,
            _ => clin_sig = ClinSig::Unknown,
        }
        Ok(clin_sig)
    }
}

impl From<&Record> for Vec<SecondStageRecord> {
    fn from(record: &Record) -> Self {
        let mut gene_vec = Vec::new();
        let dna_alterations: Vec<_> = record.dna_alterations.iter().unique().collect();
        let protein_alterations: Vec<_> = record.protein_alterations.iter().unique().collect();
        for (i, protein_alt) in protein_alterations.into_iter().enumerate() {
            let alt = if protein_alt.is_empty() {
                SecondStageRecord {
                    sample: record.sample.to_owned(),
                    alteration: dna_alterations[i].to_owned(),
                    variants: record.variants.iter().sorted().unique().join("/"),
                }
            } else {
                SecondStageRecord {
                    sample: record.sample.to_owned(),
                    alteration: protein_alt.to_owned(),
                    variants: record.variants.iter().sorted().unique().join("/"),
                }
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
            variants: record.variants.iter().sorted().unique().join("/"),
        }
    }
}

fn make_final_bar_plot_records(records: &[BarPlotRecord]) -> Vec<Counter> {
    let mut count_map = HashMap::new();
    for i in records {
        let r = count_map.entry((&i.key, &i.value)).or_insert_with(|| 0);
        *r += 1;
    }

    let mut res = Vec::new();

    for ((alt, imp), count) in count_map {
        let plot_rec = BarPlotRecord {
            key: alt.to_owned(),
            value: imp.to_owned(),
        };
        let record = Counter {
            count,
            record: plot_rec,
        };

        res.push(record);
    }

    res
}

fn order_by_impact(impacts: Vec<&Counter>) -> HashMap<String, Vec<Impact>> {
    let mut order = HashMap::new();
    let mut order_tuples = HashMap::new();
    for c in impacts {
        let impact = Impact::from_str(&c.record.value).unwrap();
        let rec = order_tuples
            .entry(c.record.key.to_owned())
            .or_insert_with(Vec::new);
        rec.push((impact, c.count))
    }

    for v in order_tuples.values_mut() {
        v.sort_by(|(_, a), (_, b)| a.cmp(&b))
    }

    for (k, v) in order_tuples {
        let removed_count = v.into_iter().map(|(x, _)| x).collect();
        order.insert(k, removed_count);
    }
    order
}

fn order_by_clin_sig(clin_sigs: Vec<&Counter>) -> HashMap<String, Vec<ClinSig>> {
    let mut order = HashMap::new();
    let mut order_tuples = HashMap::new();
    for c in clin_sigs {
        let impact = ClinSig::from_str(&c.record.value).unwrap();
        let rec = order_tuples
            .entry(c.record.key.to_owned())
            .or_insert_with(Vec::new);
        rec.push((impact, c.count))
    }

    for v in order_tuples.values_mut() {
        v.sort_by(|(_, a), (_, b)| a.cmp(&b))
    }

    for (k, v) in order_tuples {
        let removed_count = v.into_iter().map(|(x, _)| x).collect();
        order.insert(k, removed_count);
    }
    order
}
