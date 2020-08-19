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

                        let get_field = |str| {
                            str::from_utf8(
                                fields[*ann_indices.get(&String::from(str)).expect(
                                    &("No field named ".to_owned()
                                        + str
                                        + " found. Please only use VEP-annotated VCF-files."),
                                )],
                            )
                        };

                        let impact = get_field("IMPACT")?;
                        let clin_sig = get_field("CLIN_SIG")?;
                        let gene = get_field("Gene")?;
                        let dna_alteration = get_field("HGVSc")?;
                        let protein_alteration = get_field("HGVSp")?;
                        let consequence = get_field("Consequence")?;

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

                        let imp_rec = impacts.entry(gene.to_owned()).or_insert_with(|| Vec::new());
                        imp_rec.push(BarPlotRecord::new(gene.to_owned(), impact.to_owned()));

                        let cons_rec = consequences
                            .entry(gene.to_owned())
                            .or_insert_with(|| Vec::new());
                        cons_rec.push(BarPlotRecord::new(gene.to_owned(), consequence.to_owned()));

                        let alt = if protein_alteration.is_empty() {
                            dna_alteration
                        } else {
                            protein_alteration
                        };
                        let gene_imp_rec = gene_impacts
                            .entry(gene.to_owned())
                            .or_insert_with(|| Vec::new());
                        gene_imp_rec.push(BarPlotRecord::new(alt.to_owned(), impact.to_owned()));

                        let gene_cons_rec = gene_consequences
                            .entry(gene.to_owned())
                            .or_insert_with(|| Vec::new());
                        gene_cons_rec
                            .push(BarPlotRecord::new(alt.to_owned(), consequence.to_owned()));

                        let gene_clin_sig_rec = gene_clin_sigs
                            .entry(gene.to_owned())
                            .or_insert_with(|| Vec::new());

                        let clin_rec = clin_sigs
                            .entry(gene.to_owned())
                            .or_insert_with(|| Vec::new());
                        let sigs: Vec<_> = clin_sig.split('&').collect();
                        for s in sigs {
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

                            let f = gene_af_data
                                .entry(gene.to_owned())
                                .or_insert_with(|| Vec::new());
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
                .or_insert_with(|| Vec::new());
            e.push(final_gene_impacts);

            // data for second stage consequence
            let gene_consequence = gene_consequences.get(gene).unwrap();
            let final_gene_consequences = make_final_bar_plot_records(gene_consequence);
            let g = gene_consequence_data
                .entry(gene.to_owned())
                .or_insert_with(|| Vec::new());
            g.push(final_gene_consequences);

            // data for second stage clin_sig
            let gene_clin_sig = gene_clin_sigs.get(gene).unwrap();
            let final_gene_clin_sigs = make_final_bar_plot_records(gene_clin_sig);
            let f = gene_clin_sig_data
                .entry(gene.to_owned())
                .or_insert_with(|| Vec::new());
            f.push(final_gene_clin_sigs);
        }
    }

    let gene_path = output_path.to_owned() + "/genes/";
    fs::create_dir(Path::new(&gene_path))?;
    let mut gene_templates = Tera::default();
    gene_templates.register_filter("embed_source", embed_source);
    gene_templates.add_raw_template("genes.html.tera", include_str!("genes.html.tera"))?;

    let gene_specs: Value = serde_json::from_str(include_str!("gene_specs.json")).unwrap();

    // create html for second stage
    for (gene, data) in gene_data {
        create_second_stage(
            gene,
            data,
            &gene_impact_data,
            &gene_consequence_data,
            &gene_clin_sig_data,
            &gene_af_data,
            &gene_templates,
            &gene_specs,
            &gene_path,
        )?;
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
    let consequence_data: Vec<_> = consequence_data.iter().flatten().collect();
    let clin_sig_data: Vec<_> = clin_sig_data.iter().flatten().collect();

    let page_size = 100;

    let mut v = Vec::from_iter(sort_genes);
    v.sort_by(|(a, b), (c, d)| if b == d { a.cmp(&c) } else { d.cmp(&b) });
    let ordered_genes: Vec<_> = v.iter().map(|(x, _)| x).collect();
    let order = json!(ordered_genes);

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

        let mut vl_specs: Value = serde_json::from_str(include_str!("report_specs.json")).unwrap();
        let values = json!({"main": page_data , "impact": impact_page_data, "consequence": consequence_page_data , "clin_sig": clin_sig_page_data, "allel_frequency": af_page_data});

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

fn create_second_stage(
    gene: String,
    data: Vec<SecondStageRecord>,
    gene_impact_data: &HashMap<String, Vec<Vec<Counter>>>,
    gene_consequence_data: &HashMap<String, Vec<Vec<Counter>>>,
    gene_clin_sig_data: &HashMap<String, Vec<Vec<Counter>>>,
    gene_af_data: &HashMap<String, Vec<AlleleFrequency>>,
    gene_templates: &Tera,
    gene_specs: &Value,
    gene_path: &String,
) -> Result<(), Box<dyn Error>> {
    let gene_data: Vec<_> = data.iter().sorted().collect();
    let mut alterations = HashMap::new();
    for rec in &gene_data {
        let entry = alterations
            .entry(&rec.alteration)
            .or_insert_with(|| Vec::new());
        entry.push(rec.sample.to_owned());
    }
    let mut sort_alterations = HashMap::new();
    for (alteration, mut samples) in alterations {
        samples.sort();
        samples.dedup();
        sort_alterations.insert(alteration, samples.len());
    }
    let mut order = Vec::from_iter(sort_alterations);
    order.sort_by(|(a, b), (c, d)| if b == d { a.cmp(&c) } else { d.cmp(&b) });
    let order: Vec<_> = order.iter().map(|(x, _)| x).collect();
    let mut specs = gene_specs.clone();
    let impact_data = gene_impact_data.get(&gene).unwrap();
    let final_impact: Vec<_> = impact_data.iter().flatten().sorted().collect();
    let consequence_data = gene_consequence_data.get(&gene).unwrap();
    let final_consequence: Vec<_> = consequence_data.iter().flatten().sorted().collect();
    let clin_sig_data = gene_clin_sig_data.get(&gene).unwrap();
    let final_clin_sig: Vec<_> = clin_sig_data.iter().flatten().sorted().collect();
    let allel_frequency_data = gene_af_data.get(&gene).unwrap();
    let values = json!({ "main": gene_data, "impact": final_impact, "consequence": final_consequence, "clin_sig": final_clin_sig, "allel_frequency": allel_frequency_data});
    specs["datasets"] = values;
    let mut packer = Packer::new();
    let options = PackOptions::new();
    let packed_gene_specs = packer.pack(&specs, &options).unwrap();
    let mut context = Context::new();
    context.insert("genespecs", &serde_json::to_string(&packed_gene_specs)?);
    context.insert("gene", &gene);
    context.insert("order", &serde_json::to_string(&json!(order))?);
    let local: DateTime<Local> = Local::now();
    context.insert("time", &local.format("%a %b %e %T %Y").to_string());
    context.insert("version", &env!("CARGO_PKG_VERSION"));
    let html = gene_templates.render("genes.html.tera", &context)?;
    let filepath = String::from(gene_path.clone()) + &gene + ".html";
    let mut file = File::create(filepath)?;
    file.write_all(html.as_bytes())?;
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
    dna_alterations: String,
    protein_alterations: String,
    variants: String,
}

// Record for the second stage of the report
#[derive(Serialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct SecondStageRecord {
    sample: String,
    gene: String,
    alteration: String,
    variants: String,
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
                    gene: record.gene.to_owned(),
                    alteration: dna_alterations[i].to_owned(),
                    variants: record.variants.iter().sorted().unique().join("/"),
                }
            } else {
                SecondStageRecord {
                    sample: record.sample.to_owned(),
                    gene: record.gene.to_owned(),
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
            dna_alterations: record.dna_alterations.iter().sorted().unique().join(", "),
            protein_alterations: record
                .protein_alterations
                .iter()
                .sorted()
                .unique()
                .join(", "),
            variants: record.variants.iter().sorted().unique().join("/"),
        }
    }
}

fn make_final_bar_plot_records(records: &Vec<BarPlotRecord>) -> Vec<Counter> {
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
            count: count,
            record: plot_rec,
        };

        res.push(record);
    }

    res
}
