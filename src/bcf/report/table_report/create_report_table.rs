use crate::bcf::report::table_report::fasta_reader::{get_fasta_lengths, read_fasta};
use crate::bcf::report::table_report::static_reader::{get_static_reads, Variant};
use crate::common::Region;
use anyhow::Context as AnyhowContext;
use anyhow::Result;
use chrono::{DateTime, Local};
use itertools::Itertools;
use jsonm::packer::{PackOptions, Packer};
use log::warn;
use lz_str::compress_to_utf16;
use rust_htslib::bcf::header::{HeaderView, TagType};
use rust_htslib::bcf::{HeaderRecord, Read, Record};
use rustc_serialize::json::Json;
use serde::Serialize;
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use tera::{Context, Tera};

#[derive(Serialize, Clone, Debug, PartialEq)]
pub enum VariantType {
    Deletion,
    Insertion,
    Duplicate,
    Inversion,
    Variant,
}

#[derive(Serialize, Clone, Debug, PartialEq)]
pub struct Report {
    id: String,
    name: String,
    position: i64,
    reference: String,
    var_type: VariantType,
    alternatives: Option<String>,
    ann: Option<Vec<Vec<String>>>,
    format: Option<BTreeMap<String, BTreeMap<String, Value>>>,
    info: Option<HashMap<String, Vec<Value>>>,
    json_format: Option<String>,
    json_info: Option<String>,
    vis: BTreeMap<String, String>,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn make_table_report(
    vcf_path: &Path,
    fasta_path: &Path,
    bam_sample_path: &[(String, String)],
    infos: Option<Vec<String>>,
    formats: Option<Vec<String>>,
    sample: String,
    output_path: &str,
    max_read_depth: u32,
    js_files: Vec<String>,
) -> Result<()> {
    // HashMap<gene: String, Vec<Report>>, Vec<ann_field_identifiers: String>
    // let mut reports = HashMap::new();
    let mut ann_indices = HashMap::new();
    let mut vcf = rust_htslib::bcf::Reader::from_path(&vcf_path).unwrap();
    let header = vcf.header().clone();
    let header_records = header.header_records();
    let ann_field_description: Vec<_> = get_ann_description(header_records).unwrap();
    let samples: Vec<_> = header
        .samples()
        .iter()
        .map(|s| std::str::from_utf8(s).map(|s| s.to_owned()))
        .collect::<Result<_, _>>()?;

    for (i, field) in ann_field_description.iter().enumerate() {
        ann_indices.insert(field, i);
    }

    let reference_lengths = get_fasta_lengths(fasta_path)?;

    for v in vcf.records() {
        let mut variant = v.unwrap();

        let n = header.rid2name(variant.rid().unwrap()).unwrap().to_owned();
        let i = variant.id();

        let chrom = String::from_utf8(n).unwrap();

        let id = String::from_utf8(i).unwrap();

        let pos = variant.pos();
        let end_pos = match variant.info(b"END").integer() {
            Ok(Some(end_pos)) => {
                // Subtraction of 0.5 because of the 0-based positioning in the whole plot
                let end_pos = end_pos[0] as f64 - 0.5; // -1 due to 0-basing, + 0.5 du to end pos
                Some(end_pos)
            }
            _ => None,
        };

        let (info_tags, json_info_tags) = if infos.is_some() {
            let mut info_map = HashMap::new();
            for tag in infos.clone().unwrap() {
                if tag.chars().last().unwrap().eq(&'*') {
                    let prefix_tags = header
                        .header_records()
                        .iter()
                        .filter_map(|header_record| match header_record {
                            HeaderRecord::Info { key: _, values } => {
                                if values["ID"].starts_with(&tag[..(tag.len() - 1)]) {
                                    Some(values["ID"].to_owned())
                                } else {
                                    None
                                }
                            }
                            _ => None,
                        })
                        .collect::<Vec<String>>();
                    for prefix_tag in prefix_tags {
                        read_tag_entries(&mut info_map, &mut variant, &header, &prefix_tag)?;
                    }
                } else {
                    read_tag_entries(&mut info_map, &mut variant, &header, &tag)?;
                }
            }
            (
                Some(info_map.clone()),
                Some(serde_json::to_string(&json!(info_map))?),
            )
        } else {
            (None, None)
        };

        let (format_tags, json_format_tags) = if formats.is_some() {
            let mut format_map = BTreeMap::new();
            for tag in formats.clone().unwrap() {
                let (tag_type, _) = header.format_type(tag.as_bytes())?;
                match tag_type {
                    TagType::String => {
                        let values = variant.format(tag.as_bytes()).string()?;
                        for (i, v) in values.clone().into_iter().enumerate() {
                            let value = String::from_utf8(v.to_owned())?;
                            let entry = format_map
                                .entry(tag.to_owned())
                                .or_insert_with(BTreeMap::new);
                            entry.insert(samples[i].clone(), json!(value));
                        }
                    }
                    TagType::Float => {
                        let values = variant.format(tag.as_bytes()).float()?;
                        for (i, v) in values.iter().enumerate() {
                            let value = v.to_vec();
                            let entry = format_map
                                .entry(tag.to_owned())
                                .or_insert_with(BTreeMap::new);
                            entry.insert(samples[i].clone(), json!(value));
                        }
                    }
                    TagType::Integer => {
                        let values = variant.format(tag.as_bytes()).integer()?;
                        for (i, v) in values.iter().enumerate() {
                            let value = v.to_vec();
                            let entry = format_map
                                .entry(tag.to_owned())
                                .or_insert_with(BTreeMap::new);
                            entry.insert(samples[i].clone(), json!(value));
                        }
                    }
                    _ => {}
                }
            }
            (
                Some(format_map.clone()),
                Some(serde_json::to_string(&json!(format_map))?),
            )
        } else {
            (None, None)
        };

        let alleles: Vec<_> = variant
            .alleles()
            .into_iter()
            .map(|allele| allele.to_owned())
            .collect();

        let mut annotations = Vec::new();

        let mut alterations = Vec::new();
        let mut hgvsgs = Vec::new();

        if let Some(ann) = variant.info(b"ANN").string()? {
            for entry in ann.iter() {
                let fields = entry.split(|c| *c == b'|').collect_vec();

                let get_field = |field: &str| {
                    let field: String = std::str::from_utf8(
                        fields[*ann_indices
                            .get(&field.to_owned())
                            .context({
                                format!(
                                    "No field named {} found. Please only use VEP-annotated VCF-files.",
                                    field
                                )
                            })
                            .unwrap()],
                    ).unwrap().to_string();
                    field
                };

                let alteration = if !get_field("HGVSp").is_empty() {
                    get_field("HGVSp")
                } else if !get_field("HGVSg").is_empty() {
                    get_field("HGVSg")
                } else {
                    continue;
                };

                let allele = if !get_field("Allele").is_empty() {
                    get_field("Allele")
                } else {
                    continue;
                };

                if !get_field("HGVSg").is_empty() {
                    hgvsgs.push((allele, get_field("HGVSg")));
                }

                alterations.push(alteration);

                let mut ann_strings = Vec::new();
                for f in fields {
                    let attr = String::from_utf8(f.to_owned())?;
                    ann_strings.push(attr);
                }
                annotations.push(ann_strings);
            }
        }

        if !alleles.is_empty() {
            let ref_vec = alleles[0].to_owned();
            let ref_allele = String::from_utf8(ref_vec)?;

            let len: u8 = ref_allele.len() as u8;

            for allel in alleles.iter().skip(1) {
                let alt = allel.as_slice();
                let var_string = String::from("Variant");
                let var_type: VariantType;
                let alternatives: Option<String>;
                let end_position: f64;
                let plot_start_position;

                let hgvsg: String = if alleles.len() > 2 {
                    hgvsgs.iter().find(|(a, _)| *a == std::str::from_utf8(allel).unwrap()).context(format!("Found variant {} at position {} without HGVSg field for every given allele.", &id, &pos))?.0.to_owned()
                } else {
                    let mut unique_hgsvgs = hgvsgs.iter().map(|(_, b)| b).unique().collect_vec();
                    if unique_hgsvgs.len() > 1 {
                        warn!("Found variant {} at position {} with multiple HGVSg values and only one alternative allele.", &id, &pos);
                    }
                    unique_hgsvgs
                        .pop()
                        .context(format!(
                            "Found variant {} at position {} with no HGVSg value.",
                            &id, &pos
                        ))?
                        .to_owned()
                };

                match alt {
                    b"<DEL>" => {
                        var_type = VariantType::Deletion;
                        alternatives = None;
                        end_position = end_pos.unwrap();
                        plot_start_position = pos as f64 - 0.5;
                    }
                    b"<INV>" => {
                        var_type = VariantType::Inversion;
                        let rev: String = ref_allele.chars().rev().collect();
                        alternatives = Some(rev.clone());
                        end_position = end_pos.unwrap();
                        plot_start_position = pos as f64 - 0.5;
                    }
                    b"<DUP>" => {
                        var_type = VariantType::Duplicate;
                        let dup: String = [ref_allele.clone(), ref_allele.clone()].concat();
                        alternatives = Some(dup.clone());
                        end_position = end_pos.unwrap();
                        plot_start_position = pos as f64 - 0.5;
                    }
                    _ => {
                        let mut alt_allele = String::from("");

                        for c in alt {
                            if *c as char != '<' && *c as char != '>' {
                                alt_allele.push(*c as char);
                            }
                        }

                        match alt_allele.len() {
                            a if a < ref_allele.len() => {
                                plot_start_position = pos as f64 + 0.5; // start position + 1 due to alignment with deletions from bam (example: ref: ACTT alt: A  -> deletion is just CTT)
                                end_position = pos as f64 - 0.5 + len as f64;
                                var_type = VariantType::Deletion;
                                alternatives = Some(alt_allele.clone());
                            }
                            a if a > ref_allele.len() => {
                                plot_start_position = pos as f64;
                                end_position = pos as f64 + len as f64;
                                var_type = VariantType::Insertion;
                                alternatives = Some(alt_allele.clone());
                            }
                            _ => {
                                plot_start_position = pos as f64 - 0.5;
                                end_position = pos as f64 - 0.5 + len as f64;
                                var_type = VariantType::Variant;
                                alternatives = Some(alt_allele.clone());
                            }
                        }
                    }
                }

                let var = Variant {
                    marker_type: var_string,
                    reference: ref_allele.clone(),
                    alternatives,
                    start_position: plot_start_position + 1.0,
                    end_position: end_position + 1.0,
                    row: -1,
                    var_type,
                };

                let mut visualizations = BTreeMap::new();

                for (sample, bam) in bam_sample_path {
                    let bam_path = Path::new(bam);
                    let fasta_length = *reference_lengths.get(&chrom).unwrap();
                    let visualization: String;
                    if pos < 75 {
                        let content = create_report_data(
                            fasta_path,
                            Some(var.clone()),
                            bam_path,
                            &Region {
                                target: chrom.clone(),
                                start: 0,
                                end: end_position as u64 + 75,
                            },
                            max_read_depth,
                        )?;
                        visualization = manipulate_json(content, 0, end_position as u64 + 75)?;
                    } else if end_position as i64 + 75 >= fasta_length as i64 {
                        let content = create_report_data(
                            fasta_path,
                            Some(var.clone()),
                            bam_path,
                            &Region {
                                target: chrom.clone(),
                                start: pos as u64 - 75,
                                end: fasta_length - 1,
                            },
                            max_read_depth,
                        )?;
                        visualization =
                            manipulate_json(content, pos as u64 - 75, fasta_length - 1)?;
                    } else {
                        let content = create_report_data(
                            fasta_path,
                            Some(var.clone()),
                            bam_path,
                            &Region {
                                target: chrom.clone(),
                                start: pos as u64 - 75,
                                end: end_position as u64 + 75,
                            },
                            max_read_depth,
                        )?;
                        visualization =
                            manipulate_json(content, pos as u64 - 75, end_position as u64 + 75)?;
                    }

                    visualizations.insert(sample.to_owned(), visualization.to_string());
                }

                let report_data = Report {
                    id: id.clone(),
                    name: chrom.clone(),
                    position: pos + 1,
                    reference: ref_allele.clone(),
                    var_type: var.var_type,
                    alternatives: var.alternatives,
                    ann: Some(annotations.clone()),
                    format: format_tags.clone(),
                    info: info_tags.clone(),
                    json_format: json_format_tags.clone(),
                    json_info: json_info_tags.clone(),
                    vis: visualizations,
                };

                let escaped_hgvsg = escape_hgvsg(&hgvsg);
                let detail_path = Path::new(&output_path)
                    .join(Path::new("details"))
                    .join(Path::new(sample.as_str()));
                let local: DateTime<Local> = Local::now();

                let mut templates = Tera::default();
                templates.add_raw_template(
                    "table_report.html.tera",
                    include_str!("report_table.html.tera"),
                )?;
                let mut context = Context::new();
                context.insert("variant", &report_data);
                context.insert("hgvsg", &hgvsg);
                context.insert("escaped_hgvsg", &escaped_hgvsg);
                context.insert("description", &ann_field_description);
                context.insert("sample", &sample);
                context.insert("js_imports", &js_files);
                context.insert("time", &local.format("%a %b %e %T %Y").to_string());
                context.insert("version", &env!("CARGO_PKG_VERSION"));

                let html = templates.render("table_report.html.tera", &context)?;
                let filepath = detail_path.join(Path::new(&escaped_hgvsg).with_extension("html"));
                let mut file = File::create(filepath)?;
                file.write_all(html.as_bytes())?;

                let mut templates = Tera::default();
                templates.add_raw_template("plot.js.tera", include_str!("plot.js.tera"))?;

                let plot_path = detail_path
                    .join(Path::new("plots"))
                    .join(Path::new(&escaped_hgvsg).with_extension("js"));

                let mut plot_context = Context::new();
                plot_context.insert("variant", &report_data);
                let plot_html = templates.render("plot.js.tera", &plot_context)?;
                let mut plot_file = File::create(plot_path)?;
                plot_file.write_all(plot_html.as_bytes())?;
            }
        }
    }
    Ok(())
}

fn escape_hgvsg(hgvsg: &str) -> String {
    hgvsg.replace(".", "_").replace(">", "_").replace(":", "_")
}

pub(crate) fn get_ann_description(header_records: Vec<HeaderRecord>) -> Option<Vec<String>> {
    for rec in header_records {
        if let rust_htslib::bcf::HeaderRecord::Info { key: _, values } = rec {
            if values.get("ID").unwrap() == "ANN" {
                let description = values.get("Description").unwrap();
                let fields: Vec<_> = description.split('|').collect();
                let mut owned_fields = Vec::new();
                for mut entry in fields {
                    entry = entry.trim();
                    entry = entry.trim_start_matches(
                        &"\"Consequence annotations from Ensembl VEP. Format: ",
                    );
                    entry = entry.trim_end_matches(&"\"");
                    entry = entry.trim();
                    owned_fields.push(entry.to_owned());
                }
                return Some(owned_fields);
            }
        }
    }
    None
}

pub(crate) fn read_tag_entries(
    info_map: &mut HashMap<String, Vec<Value>>,
    variant: &mut Record,
    header: &HeaderView,
    tag: &str,
) -> Result<()> {
    let (tag_type, _) = header.info_type(tag.as_bytes())?;
    match tag_type {
        TagType::String => {
            if let Some(values) = variant.info(tag.as_bytes()).string()? {
                for v in values.iter() {
                    let value = String::from_utf8(Vec::from(v.to_owned()))?;
                    let entry = info_map.entry(tag.to_owned()).or_insert_with(Vec::new);
                    entry.push(json!(value));
                }
            }
        }
        TagType::Float => {
            if let Some(values) = variant.info(tag.as_bytes()).float()? {
                for v in values.iter() {
                    let entry = info_map.entry(tag.to_owned()).or_insert_with(Vec::new);
                    entry.push(json!(v));
                }
            }
        }
        TagType::Integer => {
            if let Some(values) = variant.info(tag.as_bytes()).integer()? {
                for v in values.iter() {
                    let entry = info_map.entry(tag.to_owned()).or_insert_with(Vec::new);
                    entry.push(json!(v));
                }
            }
        }
        _ => {}
    }
    Ok(())
}

pub(crate) fn create_report_data<P: AsRef<Path> + std::fmt::Debug>(
    fasta_path: P,
    variant: Option<Variant>,
    bam_path: P,
    region: &Region,
    max_read_depth: u32,
) -> Result<Json> {
    let mut data = Vec::new();

    for f in read_fasta(&fasta_path, region, true)? {
        let nucleobase = json!(f);
        data.push(nucleobase);
    }

    let (bases, matches) = get_static_reads(
        bam_path,
        fasta_path,
        region,
        max_read_depth,
        variant.as_ref(),
    )?;

    for b in bases {
        let base = json!(b);
        data.push(base);
    }

    for m in matches {
        let mat = json!(m);
        data.push(mat);
    }

    if variant.is_some() {
        data.push(json!(variant));
    }

    Ok(Json::from_str(&json!(data).to_string()).unwrap())
}

/// Inserts the json containing the genome data into the vega specs.
/// It also changes keys and values of the json data for the vega plot to look better and compresses the json with jsonm.
pub(crate) fn manipulate_json(data: Json, from: u64, to: u64) -> Result<String> {
    let json_string = include_str!("vegaSpecs.json");

    let mut vega_specs: Value = serde_json::from_str(json_string)?;
    let values: Value = serde_json::from_str(&data.to_string())?;
    let mut values = json!({"values": values, "name": "fasta"});

    let v = values["values"].as_array().unwrap().clone();

    let mut row = 0;

    for (i, _) in v.iter().enumerate() {
        let k = v[i]["marker_type"].clone().as_str().unwrap().to_owned();
        let r = v[i]["row"].clone().as_i64().unwrap();
        if r > row {
            row = r;
        }

        if k == "A" || k == "T" || k == "G" || k == "C" || k == "U" || k == "N" {
            values["values"][i]["base"] = values["values"][i]["marker_type"].clone();
        } else if k == "Deletion"
            || k == "Match"
            || k == "Pairing"
            || k == "Duplicate"
            || k == "Inversion"
        {
            values["values"][i]["typ"] = values["values"][i]["marker_type"].clone();
        } else if k == "Insertion" {
            values["values"][i]["typ"] = values["values"][i]["marker_type"].clone();
            values["values"][i]["inserts"] = values["values"][i]["bases"].clone();
        }
    }

    vega_specs["width"] = json!(700);
    vega_specs["height"] = json!(core::cmp::max(10 * row + 60, 203));
    let domain = json!([from, to]);

    vega_specs["scales"][0]["domain"] = domain;
    vega_specs["data"][1] = values;

    let mut packer = Packer::new();
    packer.set_max_dict_size(100000);
    let options = PackOptions::new();
    let packed_specs = packer.pack(&vega_specs, &options)?;

    Ok(json!(compress_to_utf16(&packed_specs.to_string())).to_string())
}
