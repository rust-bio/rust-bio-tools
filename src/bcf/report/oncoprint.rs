use std::collections::HashMap;
use std::io::Write;
use std::{fs, str};

use derive_new::new;
use itertools::Itertools;
use lazy_static::lazy_static;
use regex::Regex;
use serde_derive::Serialize;
use tera::{self, Context, Tera};

use crate::bcf::report::table_report::create_report_table::get_ann_description;
use crate::bcf::report::table_report::create_report_table::read_tag_entries;
use anyhow::Context as AnyhowContext;
use anyhow::Result;
use chrono::{DateTime, Local};
use jsonm::packer::{PackOptions, Packer};
use lz_str::compress_to_utf16;
use rust_htslib::bcf::{self, Read};
use serde_json::{json, Value};
use std::fs::File;
use std::path::Path;
use std::str::FromStr;
use thiserror::Error;

lazy_static! {
    static ref HGVSP_PROTEIN_RE: Regex = Regex::new(r"ENSP[0-9]+(\.[0-9]+)?:").unwrap();
}

pub fn oncoprint(
    sample_calls: &HashMap<String, String>,
    output_path: &str,
    max_cells: u32,
    tsv_data_path: Option<&str>,
    plot_info: Option<Vec<String>>,
) -> Result<()> {
    let mut data = HashMap::new();
    let mut gene_data = HashMap::new();
    let mut impact_data = Vec::new();
    let mut gene_impact_data = HashMap::new();
    let mut existing_var_data = Vec::new();
    let mut gene_existing_var_data = HashMap::new();
    let mut consequence_data = Vec::new();
    let mut gene_consequence_data = HashMap::new();
    let mut clin_sig_data = Vec::new();
    let mut gene_clin_sig_data = HashMap::new();
    let mut af_data = Vec::new();
    let mut gene_af_data = HashMap::new();
    let mut unique_genes = HashMap::new();
    let mut plot_info_data = HashMap::new();
    let mut gene_plot_info_data = HashMap::new();
    let mut remove_existing_variation = true;

    let tsv_data = if let Some(tsv) = tsv_data_path {
        Some(make_tsv_records(tsv.to_owned())?)
    } else {
        None
    };

    // Check every VCF file for presence of CLIN_SIG
    let mut clin_sig_present = HashMap::new();
    for (sample, path) in sample_calls.iter().sorted() {
        let bcf_reader = bcf::Reader::from_path(path)?;
        let header_records = bcf_reader.header().header_records();
        let ann_fields: Vec<_> = get_ann_description(header_records).unwrap_or_else(|| {
            panic!("No ANN field found. Please only use VEP-annotated VCF-files.")
        });
        clin_sig_present.insert(
            sample.to_owned(),
            ann_fields.contains(&"CLIN_SIG".to_owned()),
        );
    }

    // Check wether any of the VCF files contain CLIN_SIG at all
    let cs_present_folded = clin_sig_present.iter().fold(false, |b, (_, c)| b || *c);

    for (sample, path) in sample_calls.iter().sorted() {
        let mut genes = HashMap::new();
        let mut impacts = HashMap::new();
        let mut gene_impacts = HashMap::new();
        let mut existing_variations = HashMap::new();
        let mut gene_existing_variations = HashMap::new();
        let mut consequences = HashMap::new();
        let mut gene_consequences = HashMap::new();
        let mut clin_sigs = HashMap::new();
        let mut gene_clin_sigs = HashMap::new();
        let mut ann_indices = HashMap::new();
        let mut pi_data = HashMap::new();
        let mut gene_pi_data = HashMap::new();

        let mut bcf_reader = bcf::Reader::from_path(path)?;
        let header = bcf_reader.header().clone();
        let mut sample_names = Vec::new();
        for s in header.samples() {
            sample_names.push(String::from_utf8(s.to_owned())?);
        }
        let header_records = header.header_records();
        let ann_fields: Vec<_> = get_ann_description(header_records).unwrap_or_else(|| {
            panic!("No ANN field found. Please only use VEP-annotated VCF-files.")
        });

        for (i, field) in ann_fields.iter().enumerate() {
            ann_indices.insert(field, i);
        }

        let clin_sig_pres = clin_sig_present.get(sample).unwrap();

        for res in bcf_reader.records() {
            let mut gene_data_per_record = HashMap::new();
            let mut record = res?;
            let alleles = record
                .alleles()
                .into_iter()
                .map(|allele| allele.to_owned())
                .collect_vec();
            let alt_alleles = &alleles[1..];
            let ref_allele = alleles[0].to_owned();
            let dup = [ref_allele.clone(), ref_allele.clone()].concat();
            let rev = ref_allele.clone().into_iter().rev().collect_vec();

            let mut info_map = HashMap::new();
            if plot_info.is_some() {
                for tag in &plot_info.clone().unwrap() {
                    read_tag_entries(&mut info_map, &mut record, &header, tag)?;
                }
            }

            let allel_frequencies = record
                .format(b"AF")
                .float()?
                .iter()
                .map(|s| s.to_vec())
                .collect_vec();

            let ann = record.info(b"ANN").string()?;
            if let Some(ann) = ann {
                for alt_allele in alt_alleles {
                    let variant = if alt_allele == b"<BND>" {
                        "BND"
                    } else if alt_allele == b"<DUP>" || alt_allele == &dup {
                        "DUP"
                    } else if alt_allele == b"<DEL>"
                        || alt_allele.len() == 1 && ref_allele.len() > 1
                    {
                        "DEL"
                    } else if alt_allele == b"<INS>"
                        || alt_allele.len() > 1 && ref_allele.len() == 1
                    {
                        "INS"
                    } else if alt_allele.len() == 1 && ref_allele.len() == 1 {
                        "SNV"
                    } else if alt_allele == b"<INV>" || alt_allele == &rev {
                        "INV"
                    } else if alt_allele.len() == ref_allele.len() {
                        "MNV"
                    } else {
                        "Replacement"
                    };

                    for entry in ann.iter() {
                        let fields: Vec<_> = entry.split(|c| *c == b'|').collect();

                        let get_field = |field: &str| {
                            str::from_utf8(
                                fields[*ann_indices.get(&field.to_owned()).unwrap_or_else(|| panic!("No field named {} found. Please only use VEP-annotated VCF-files.", field))],
                            )
                        };

                        let mut impact = get_field("IMPACT")?;
                        let clin_sig = if *clin_sig_pres {
                            get_field("CLIN_SIG")?
                        } else {
                            ""
                        };
                        let gene = if !get_field("SYMBOL")?.is_empty() {
                            get_field("SYMBOL")?
                        } else if !get_field("Gene")?.is_empty() {
                            get_field("Gene")?
                        } else if !get_field("HGVSg")?.is_empty() {
                            get_field("HGVSg")?
                        } else {
                            continue;
                        };
                        let dna_alteration = get_field("HGVSg")?;
                        let canonical =
                            if let Some(index) = ann_indices.get(&"CANONICAL".to_owned()) {
                                str::from_utf8(fields[*index])? == "YES"
                            } else {
                                false
                            };
                        let protein_alteration = get_field("HGVSp")?;
                        let consequence = get_field("Consequence")?;
                        let existing_var = get_field("Existing_variation")?;

                        let gene_rec = unique_genes.entry(gene.to_owned()).or_insert_with(Vec::new);
                        gene_rec.push(sample.to_owned());

                        let rec = genes
                            .entry(gene.to_owned())
                            .or_insert_with(|| Record::new(sample.to_owned(), gene.to_owned()));

                        // Data for second stage including whether the record is marked as canonical or not.
                        let gene_entry = gene_data_per_record
                            .entry(gene.to_owned())
                            .or_insert_with(&Vec::new);
                        gene_entry.push((
                            SecondStageRecord {
                                sample: rec.sample.clone(),
                                alteration: if protein_alteration.is_empty() {
                                    dna_alteration.to_owned()
                                } else {
                                    protein_alteration.to_owned()
                                },
                                variant: variant.to_owned(),
                            },
                            canonical,
                        ));

                        rec.variants.push(variant.to_owned());

                        let imp_rec = impacts.entry(gene.to_owned()).or_insert_with(Vec::new);
                        if impact.is_empty() {
                            impact = "unknown";
                        }
                        imp_rec.push(BarPlotRecord::new(gene.to_owned(), impact.to_owned()));

                        let ev_rec = existing_variations
                            .entry(gene.to_owned())
                            .or_insert_with(Vec::new);
                        let gene_ev_rec = gene_existing_variations
                            .entry(gene.to_owned())
                            .or_insert_with(Vec::new);

                        let alt = if protein_alteration.is_empty() {
                            dna_alteration
                        } else {
                            protein_alteration
                        };

                        if plot_info.is_some() {
                            for key in plot_info.clone().unwrap() {
                                let info = info_map.get(&key.clone());
                                if info.is_none() {
                                    let e = pi_data.entry(key.clone()).or_insert_with(HashMap::new);
                                    let rec = e.entry(gene.to_owned()).or_insert_with(Vec::new);
                                    rec.push(BarPlotRecord::new(
                                        gene.to_owned(),
                                        "unknown".to_string(),
                                    ));
                                    let e2 = gene_pi_data
                                        .entry(key.clone())
                                        .or_insert_with(HashMap::new);
                                    let rec2 = e2.entry(gene.to_owned()).or_insert_with(Vec::new);
                                    rec2.push(BarPlotRecord::new(
                                        alt.to_owned(),
                                        "unknown".to_string(),
                                    ));
                                } else {
                                    for val in info.unwrap() {
                                        let value = if val == &json!("") || val == &json!(".") {
                                            "unknown".to_string()
                                        } else {
                                            val.to_string().trim_matches('\"').to_owned()
                                        };
                                        let e =
                                            pi_data.entry(key.clone()).or_insert_with(HashMap::new);
                                        let rec = e.entry(gene.to_owned()).or_insert_with(Vec::new);
                                        rec.push(BarPlotRecord::new(
                                            gene.to_owned(),
                                            value.clone(),
                                        ));
                                        let e2 = gene_pi_data
                                            .entry(key.clone())
                                            .or_insert_with(HashMap::new);
                                        let rec2 =
                                            e2.entry(gene.to_owned()).or_insert_with(Vec::new);
                                        rec2.push(BarPlotRecord::new(
                                            alt.to_owned(),
                                            value.clone(),
                                        ));
                                    }
                                }
                            }
                        }

                        let split_ev = existing_var.split('&').collect_vec();
                        for ex_var in split_ev {
                            let mut ev: String =
                                ex_var.chars().filter(|c| !c.is_digit(10)).collect();
                            if ev.is_empty() {
                                ev = String::from("unknown");
                            } else {
                                remove_existing_variation = false;
                            }
                            ev_rec.push(BarPlotRecord::new(gene.to_owned(), ev.clone()));
                            gene_ev_rec.push(BarPlotRecord::new(alt.to_owned(), ev));
                        }

                        let gene_imp_rec =
                            gene_impacts.entry(gene.to_owned()).or_insert_with(Vec::new);
                        gene_imp_rec.push(BarPlotRecord::new(alt.to_owned(), impact.to_owned()));

                        let split_consequences: Vec<_> = consequence.split('&').collect();

                        let cons_rec = consequences.entry(gene.to_owned()).or_insert_with(Vec::new);

                        let gene_cons_rec = gene_consequences
                            .entry(gene.to_owned())
                            .or_insert_with(Vec::new);

                        for mut c in split_consequences {
                            if c.is_empty() {
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
                            if s.is_empty() {
                                s = "unknown";
                            }
                            s = s.trim_matches('_');
                            clin_rec.push(BarPlotRecord::new(gene.to_owned(), s.to_owned()));
                            gene_clin_sig_rec
                                .push(BarPlotRecord::new(alt.to_owned(), s.to_owned()));
                        }

                        for (i, name) in sample_names.iter().enumerate() {
                            for frequency in &allel_frequencies[i] {
                                let af = AlleleFrequency {
                                    sample: sample.to_owned() + ":" + name,
                                    key: gene.to_owned(),
                                    allel_frequency: *frequency,
                                };

                                af_data.push(af);

                                let gene_af = AlleleFrequency {
                                    sample: sample.to_owned() + ":" + name,
                                    key: alt.to_owned(),
                                    allel_frequency: *frequency,
                                };

                                let f =
                                    gene_af_data.entry(gene.to_owned()).or_insert_with(Vec::new);
                                f.push(gene_af);
                            }
                        }
                    }
                }
            }
            // Filter records marked with canonical. Keep all if no one is marked.
            for (k, record_tuples) in &gene_data_per_record {
                let rec = gene_data.entry(k.to_owned()).or_insert_with(&Vec::new);
                let filter_canonical = record_tuples
                    .iter()
                    .filter(|(_, canonical)| *canonical)
                    .collect_vec();
                match filter_canonical.len() {
                    0 => {
                        rec.extend(record_tuples.iter().map(|(r, _)| r.clone()));
                    }
                    1 => rec.extend(filter_canonical.iter().map(|(r, _)| r.clone())),
                    _ => panic!("Found more than one variant annotated as canonical!"),
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

            let ex_var = existing_variations.get(gene).unwrap();
            let final_evs = make_final_bar_plot_records(ex_var);
            existing_var_data.push(final_evs);

            for (tag, map) in pi_data.clone() {
                if let Some(t_data) = map.get(gene) {
                    let final_data = make_final_bar_plot_records(t_data);
                    let e = plot_info_data.entry(tag).or_insert_with(Vec::new);
                    e.push(final_data);
                }
            }

            let consequence = consequences.get(gene).unwrap();
            let final_consequences = make_final_bar_plot_records(consequence);
            consequence_data.push(final_consequences);

            // data for first stage clin_sig
            let cls = clin_sigs.get(gene).unwrap();
            let final_clin_sig = make_final_bar_plot_records(cls);
            clin_sig_data.push(final_clin_sig);

            // data for second stage impact
            let gene_impact = gene_impacts.get(gene).unwrap();
            let final_gene_impacts = make_final_bar_plot_records(gene_impact);
            let e = gene_impact_data
                .entry(gene.to_owned())
                .or_insert_with(Vec::new);
            e.push(final_gene_impacts);

            for (tag, map) in gene_pi_data.clone() {
                if let Some(t_data) = map.get(gene) {
                    let final_data = make_final_bar_plot_records(t_data);
                    let m = gene_plot_info_data.entry(tag).or_insert_with(HashMap::new);
                    let e = m.entry(gene.to_owned()).or_insert_with(Vec::new);
                    e.push(final_data);
                }
            }

            let gene_evs = gene_existing_variations.get(gene).unwrap();
            let final_gene_evs = make_final_bar_plot_records(gene_evs);
            let v = gene_existing_var_data
                .entry(gene.to_owned())
                .or_insert_with(Vec::new);
            v.push(final_gene_evs);

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
    fs::create_dir(Path::new(&gene_path)).context(WriteErr::CantCreateDir {
        dir_path: gene_path.to_owned(),
    })?;

    let gene_plots_path = output_path.to_owned() + "/genes/plots/";
    fs::create_dir(Path::new(&gene_plots_path)).context(WriteErr::CantCreateDir {
        dir_path: gene_plots_path.to_owned(),
    })?;

    let mut gene_templates = Tera::default();
    gene_templates.add_raw_template("genes.html.tera", include_str!("genes.html.tera"))?;
    gene_templates.add_raw_template("plots.js.tera", include_str!("plots.js.tera"))?;

    let gene_specs: Value = serde_json::from_str(include_str!("gene_specs.json"))?;

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

        let impact_data = gene_impact_data.get(&gene).unwrap();
        let final_impact: Vec<_> = impact_data.iter().flatten().sorted().collect();
        let existing_var_data = gene_existing_var_data.get(&gene).unwrap();
        let final_ev: Vec<_> = existing_var_data.iter().flatten().sorted().collect();
        let mut inf_plot_data = HashMap::new();
        for (tag, data) in &gene_plot_info_data {
            if let Some(d) = data.get(&gene) {
                let final_data: Vec<_> = d.iter().flatten().sorted().collect();
                inf_plot_data.insert(tag.to_owned(), final_data.clone());
            }
        }
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

                let ev_page_data: Vec<_> = final_ev
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

                let mut info_page_data = HashMap::new();
                if plot_info.is_some() {
                    for (tag, data) in inf_plot_data.clone() {
                        info_page_data.insert(
                            tag,
                            data.into_iter()
                                .filter(|entry| sorted_alterations.contains(&&&entry.record.key))
                                .collect_vec(),
                        );
                    }
                }

                let order: Vec<_> = ordered_alts
                    .iter()
                    .filter(|gene| sorted_alterations.contains(gene))
                    .collect();

                let samples: Vec<_> = page_data.iter().map(|r| r.sample.clone()).collect();
                let unique_samples: Vec<_> = samples.iter().unique().collect();

                let mut specs = gene_specs.clone();

                let mut values = if cs_present_folded {
                    json!({ "main": page_data, "impact": impact_page_data, "ev": ev_page_data, "consequence": consequence_page_data, "clin_sig": clin_sig_page_data, "allel_frequency": af_page_data})
                } else {
                    json!({ "main": page_data, "impact": impact_page_data, "ev": ev_page_data, "consequence": consequence_page_data, "allel_frequency": af_page_data})
                };

                if plot_info.is_some() {
                    for (tag, data) in info_page_data {
                        values[tag] = json!(data);
                    }
                }

                if let Some(ref tsv) = tsv_data {
                    for (title, data) in tsv {
                        values[title] = json!(data);
                    }
                }

                specs["datasets"] = values;
                if !cs_present_folded || remove_existing_variation {
                    let hconcat = specs["vconcat"][1]["hconcat"].as_array_mut().unwrap();
                    match (!cs_present_folded, remove_existing_variation) {
                        (true, true) => {
                            hconcat.remove(6);
                            hconcat.remove(4)
                        }
                        (true, false) => hconcat.remove(4),
                        (false, true) => hconcat.remove(6),
                        (_, _) => unreachable!(),
                    };
                    specs["vconcat"][1]["hconcat"] = json!(hconcat);
                }

                if plot_info.is_some() {
                    let info_specs: Value = serde_json::from_str(include_str!("info_specs.json"))?;
                    let highlight_specs: Value =
                        serde_json::from_str(include_str!("highlight_specs.json"))?;
                    let hconcat = specs["vconcat"][1]["hconcat"].as_array_mut().unwrap();
                    for tag in plot_info_data.keys() {
                        let mut tag_specs = info_specs.clone();
                        tag_specs["data"] = json!({ "name": tag });
                        let highlight_name = "highlight_".to_string() + tag;
                        tag_specs["selection"] = json!({ &highlight_name: highlight_specs });
                        tag_specs["encoding"]["color"]["title"] = json!(tag);
                        tag_specs["encoding"]["x"]["title"] = json!(tag);
                        tag_specs["encoding"]["fillOpacity"]["condition"]["selection"] =
                            json!(highlight_name);
                        hconcat.push(tag_specs);
                    }

                    specs["vconcat"][1]["hconcat"] = json!(hconcat);
                }

                if let Some(ref tsv) = tsv_data {
                    let tsv_specs: Value = serde_json::from_str(include_str!("tsv_specs.json"))?;

                    for title in tsv.keys() {
                        let mut tsv_plot = tsv_specs.clone();
                        tsv_plot["data"] = json!({ "name": title });
                        tsv_plot["encoding"]["color"]["title"] = json!(title);
                        let vconcat = specs["vconcat"].as_array_mut().unwrap();
                        vconcat.insert(1, tsv_plot);
                        specs["vconcat"] = json!(vconcat);
                    }
                }

                let mut packer = Packer::new();
                let options = PackOptions::new();
                let packed_gene_specs = packer.pack(&specs, &options)?;
                let mut context = Context::new();
                let oncoprint = json!(compress_to_utf16(&serde_json::to_string(
                    &packed_gene_specs
                )?))
                .to_string();
                context.insert("oncoprint", &oncoprint);
                context.insert("gene", &gene);
                context.insert("samples", &unique_samples.len());
                context.insert("current_page", &page);
                context.insert("pages", &(pages + 1));
                context.insert("order", &serde_json::to_string(&json!(order))?);
                let local: DateTime<Local> = Local::now();
                context.insert("time", &local.format("%a %b %e %T %Y").to_string());
                context.insert("version", &env!("CARGO_PKG_VERSION"));
                let html = gene_templates.render("genes.html.tera", &context)?;
                let js = gene_templates.render("plots.js.tera", &context)?;
                let filepath = gene_path.clone() + &gene + &page.to_string() + ".html";
                let js_filepath = gene_plots_path.clone() + &gene + &page.to_string() + ".js";
                let mut file = File::create(filepath)?;
                let mut js_file = File::create(js_filepath)?;
                file.write_all(html.as_bytes())?;
                js_file.write_all(js.as_bytes())?;
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
    let ev_data: Vec<_> = existing_var_data.iter().flatten().collect();
    let consequence_data: Vec<_> = consequence_data.iter().flatten().collect();
    let clin_sig_data: Vec<_> = clin_sig_data.iter().flatten().collect();
    let mut i_plot_data = HashMap::new();
    if plot_info.is_some() {
        for (tag, data) in plot_info_data.clone() {
            i_plot_data.insert(tag, data.into_iter().flatten().collect_vec());
        }
    }

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

    let index_path = output_path.to_owned() + "/indexes";
    fs::create_dir(Path::new(&index_path)).context(WriteErr::CantCreateDir {
        dir_path: index_path.to_owned(),
    })?;

    let prefixes = make_prefixes(ordered_genes.clone(), page_size);
    let prefix_path = output_path.to_owned() + "/prefixes/";
    fs::create_dir(Path::new(&prefix_path)).context(WriteErr::CantCreateDir {
        dir_path: prefix_path.to_owned(),
    })?;

    let mut templates = Tera::default();
    templates.add_raw_template(
        "prefix_table.html.tera",
        include_str!("prefix_table.html.tera"),
    )?;
    let mut context = Context::new();
    context.insert("table", &prefixes);
    let html = templates.render("prefix_table.html.tera", &context)?;

    let file_path = output_path.to_owned() + "/prefixes/prefixes.html";
    let mut file = File::create(file_path)?;
    file.write_all(html.as_bytes())?;

    let gene_path = prefix_path + "/genes/";
    fs::create_dir(Path::new(&gene_path)).context(WriteErr::CantCreateDir {
        dir_path: gene_path.to_owned(),
    })?;

    for (prefix, values) in prefixes {
        let mut templates = Tera::default();
        templates.add_raw_template(
            "lookup_table.html.tera",
            include_str!("lookup_table.html.tera"),
        )?;
        let mut context = Context::new();
        context.insert("values", &values);
        let html = templates.render("lookup_table.html.tera", &context)?;

        let file_path = gene_path.to_owned() + &prefix + ".html";
        let mut file = File::create(file_path)?;
        file.write_all(html.as_bytes())?;
    }

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

            let ev_page_data: Vec<_> = ev_data
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

            let mut info_page_data = HashMap::new();
            if plot_info.is_some() {
                for (tag, data) in i_plot_data.clone() {
                    info_page_data.insert(
                        tag,
                        data.into_iter()
                            .filter(|entry| sorted_genes.contains(&&entry.record.key))
                            .collect_vec(),
                    );
                }
            }

            let order: Vec<_> = ordered_genes
                .iter()
                .filter(|gene| sorted_genes.contains(gene))
                .collect();

            let samples: Vec<_> = page_data.iter().map(|r| r.sample.clone()).collect();
            let unique_samples: Vec<_> = samples.iter().unique().collect();

            let mut vl_specs: Value = serde_json::from_str(include_str!("report_specs.json"))?;

            let mut values = if cs_present_folded {
                json!({ "main": page_data, "impact": impact_page_data, "ev": ev_page_data, "consequence": consequence_page_data, "clin_sig": clin_sig_page_data, "allel_frequency": af_page_data})
            } else {
                json!({ "main": page_data, "impact": impact_page_data, "ev": ev_page_data, "consequence": consequence_page_data, "allel_frequency": af_page_data})
            };

            if plot_info.is_some() {
                for (tag, data) in info_page_data {
                    values[tag] = json!(data);
                }
            }

            if let Some(ref tsv) = tsv_data {
                for (title, data) in tsv {
                    values[title] = json!(data);
                }
            }

            vl_specs["datasets"] = values;
            if !cs_present_folded || remove_existing_variation {
                let hconcat = vl_specs["vconcat"][1]["hconcat"].as_array_mut().unwrap();
                match (!cs_present_folded, remove_existing_variation) {
                    (true, true) => {
                        hconcat.remove(6);
                        hconcat.remove(4)
                    }
                    (true, false) => hconcat.remove(4),
                    (false, true) => hconcat.remove(6),
                    (_, _) => unreachable!(),
                };
                vl_specs["vconcat"][1]["hconcat"] = json!(hconcat);
            }

            if plot_info.is_some() {
                let info_specs: Value = serde_json::from_str(include_str!("info_specs.json"))?;
                let highlight_specs: Value =
                    serde_json::from_str(include_str!("highlight_specs.json"))?;
                let hconcat = vl_specs["vconcat"][1]["hconcat"].as_array_mut().unwrap();
                for tag in plot_info_data.keys() {
                    let mut tag_specs = info_specs.clone();
                    tag_specs["data"] = json!({ "name": tag });
                    let highlight_name = "highlight_".to_string() + tag;
                    tag_specs["selection"] = json!({ &highlight_name: highlight_specs });
                    tag_specs["encoding"]["color"]["title"] = json!(tag);
                    tag_specs["encoding"]["x"]["title"] = json!(tag);
                    tag_specs["encoding"]["fillOpacity"]["condition"]["selection"] =
                        json!(highlight_name);
                    hconcat.push(tag_specs);
                }

                vl_specs["vconcat"][1]["hconcat"] = json!(hconcat);
            }

            if let Some(ref tsv) = tsv_data {
                let tsv_specs: Value = serde_json::from_str(include_str!("tsv_specs.json"))?;

                for title in tsv.keys() {
                    let mut tsv_plot = tsv_specs.clone();
                    tsv_plot["data"] = json!({ "name": title });
                    tsv_plot["encoding"]["color"]["title"] = json!(title);
                    let vconcat = vl_specs["vconcat"].as_array_mut().unwrap();
                    vconcat.insert(1, tsv_plot);
                    vl_specs["vconcat"] = json!(vconcat);
                }
            }

            let mut packer = Packer::new();
            let options = PackOptions::new();
            let packed_specs = packer.pack(&vl_specs, &options)?;
            let mut templates = Tera::default();
            templates.add_raw_template("report.html.tera", include_str!("report.html.tera"))?;
            templates.add_raw_template("plots.js.tera", include_str!("plots.js.tera"))?;
            let mut context = Context::new();
            let data = json!(compress_to_utf16(&serde_json::to_string(&packed_specs)?)).to_string();
            context.insert("oncoprint", &data);
            context.insert("current_page", &page);
            context.insert("pages", &(pages + 1));
            context.insert("order", &serde_json::to_string(&json!(order))?);
            context.insert("samples", &unique_samples.len());
            let local: DateTime<Local> = Local::now();
            context.insert("time", &local.format("%a %b %e %T %Y").to_string());
            context.insert("version", &env!("CARGO_PKG_VERSION"));

            let html = templates.render("report.html.tera", &context)?;
            let js = templates.render("plots.js.tera", &context)?;

            let index = format!("{}/index{}.html", index_path, page.to_string());
            let js_index = format!("{}/plot{}.js", index_path, page.to_string());
            let mut file = File::create(index)?;
            let mut js_file = File::create(js_index)?;
            file.write_all(html.as_bytes())?;
            js_file.write_all(js.as_bytes())?;
        }
    }

    // Add index1 when no variants are found
    let index1_path = output_path.to_owned() + "/indexes/index1.html";
    if !Path::new(&index1_path).exists() {
        let mut templates = Tera::default();
        templates.add_raw_template(
            "empty_report.html",
            include_str!("html/empty_report.html.tera"),
        )?;
        let mut context = Context::new();
        let local: DateTime<Local> = Local::now();
        context.insert("time", &local.format("%a %b %e %T %Y").to_string());
        context.insert("version", &env!("CARGO_PKG_VERSION"));
        let no_variants = templates.render("empty_report.html", &context)?;
        let mut file = File::create(index1_path)?;
        file.write_all(no_variants.as_bytes())?;
    }

    Ok(())
}

#[derive(new, Debug, Clone)]
struct Record {
    sample: String,
    gene: String,
    #[new(default)]
    dna_alteration: Vec<String>,
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

#[derive(new, Serialize, Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct TSVRecord {
    sample: String,
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
    variant: String,
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
    ConflictingInterpretationsOfPathogenicity,
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
            "likely_pathogenic/pathogenic" => clin_sig = ClinSig::LikelyPathogenicPathogenic,
            "likely_pathogenic" => clin_sig = ClinSig::LikelyPathogenic,
            "risk_factor" => clin_sig = ClinSig::RiskFactor,
            "drug_response" => clin_sig = ClinSig::DrugResponse,
            "affects" => clin_sig = ClinSig::Affects,
            "association" => clin_sig = ClinSig::Association,
            "uncertain_significance" => clin_sig = ClinSig::UncertainSignificance,
            "conflicting_interpretations_of_pathogenicity" => {
                clin_sig = ClinSig::ConflictingInterpretationsOfPathogenicity
            }
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

impl From<&Record> for FinalRecord {
    fn from(record: &Record) -> Self {
        FinalRecord {
            sample: record.sample.to_owned(),
            gene: record.gene.to_owned(),
            variants: record.variants.iter().sorted().unique().join("/"),
        }
    }
}

fn make_tsv_records(tsv_path: String) -> Result<HashMap<String, Vec<TSVRecord>>> {
    let mut tsv_values = HashMap::new();
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(tsv_path)?;

    let header = rdr.headers()?.clone();
    let titles: Vec<_> = header.iter().skip(1).collect();
    for res in rdr.records() {
        let row = res?;
        let sample = row[0].to_owned();
        for (i, value) in row.iter().skip(1).enumerate() {
            let rec = tsv_values
                .entry(titles[i].to_owned())
                .or_insert_with(Vec::new);
            let entry = TSVRecord {
                sample: sample.clone(),
                value: value.to_owned(),
            };
            rec.push(entry);
        }
    }
    Ok(tsv_values)
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
        v.sort_by(|(i1, a), (i2, b)| b.cmp(&a).then(i2.cmp(&i1)))
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
        v.sort_by(|(c1, a), (c2, b)| b.cmp(&a).then(c2.cmp(&c1)));
    }

    for (k, v) in order_tuples {
        let removed_count = v.into_iter().map(|(x, _)| x).collect();
        order.insert(k, removed_count);
    }
    order
}

fn make_prefixes(
    genes: Vec<&String>,
    rows_per_page: usize,
) -> HashMap<String, Vec<(&String, usize)>> {
    let mut prefix_map = HashMap::new();
    let prefix_len = 3;
    for (i, partial_table) in genes.chunks(rows_per_page).enumerate() {
        let page = i + 1;
        for gene in partial_table {
            if gene.len() >= prefix_len {
                let prefix = gene.chars().take(prefix_len).collect::<String>();
                let entry = prefix_map.entry(prefix).or_insert_with(Vec::new);
                entry.push((gene.to_owned(), page));
            }
        }
    }
    prefix_map
}

#[derive(Error, Debug)]
pub enum WriteErr {
    #[error("could not create directory at {dir_path}")]
    CantCreateDir { dir_path: String },
}
