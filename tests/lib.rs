use bio::io::fastq;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::fs;
use std::process::Command;

/// Compare an output file to the expected output and delete the output file.
fn test_output(result: &str, expected: &str) {
    assert!(Command::new("cmp")
        .arg(result)
        .arg(expected)
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());
    fs::remove_file(result).unwrap();
}

/// Compare two fastq files, ignoring the name lines
/// Reads are sorted by their sequence, which is not 100% robust
/// if mutations/ sequencing errors are considered.
fn compare_fastq(result: &str, expected: &str) {
    let result_reader = fastq::Reader::from_file(result).unwrap();
    let mut result_recs: Vec<fastq::Record> =
        result_reader.records().filter_map(Result::ok).collect();
    result_recs.sort_by_key(|x| x.seq().to_owned());
    let expected_reader = fastq::Reader::from_file(expected).unwrap();
    let mut expected_recs: Vec<fastq::Record> =
        expected_reader.records().filter_map(Result::ok).collect();
    expected_recs.sort_by_key(|x| x.seq().to_owned());

    for (result, expected) in result_recs.iter().zip(expected_recs.iter()) {
        assert_eq!(result.seq(), expected.seq());
        assert_eq!(result.qual(), expected.qual());
    }
}

fn compare_bam(result: &str, expected: &str) {
    let mut result_reader = bam::Reader::from_path(result).unwrap();
    let mut result_recs: Vec<bam::Record> =
        result_reader.records().filter_map(Result::ok).collect();
    result_recs.sort_by_key(|x| x.seq().as_bytes().to_owned());
    let mut expected_reader = bam::Reader::from_path(expected).unwrap();
    let mut expected_recs: Vec<bam::Record> =
        expected_reader.records().filter_map(Result::ok).collect();
    expected_recs.sort_by_key(|x| x.seq().as_bytes().to_owned());
    for (result, expected) in result_recs.iter().zip(expected_recs.iter()) {
        assert_eq!(result.seq().as_bytes(), expected.seq().as_bytes());
        assert_eq!(result.qual(), expected.qual());
    }
}

#[test]
fn fastq_split() {
    assert!(Command::new("bash")
        .arg("-c")
        .arg("target/debug/rbt fastq-split tests/A.fastq tests/B.fastq < tests/test.fastq")
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());
    test_output("tests/A.fastq", "tests/expected/A.fastq");
    test_output("tests/B.fastq", "tests/expected/B.fastq");
}

#[test]
fn fastq_filter() {
    assert!(Command::new("bash")
        .arg("-c")
        .arg(
            "target/debug/rbt fastq-filter tests/ids.txt < tests/test.fastq > tests/filtered.fastq"
        )
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());
    test_output("tests/filtered.fastq", "tests/expected/B.fastq");
}

#[test]
fn bam_depth() {
    assert!(Command::new("bash")
        .arg("-c")
        .arg("target/debug/rbt bam-depth tests/test.bam < tests/pos.txt > tests/depth.txt")
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());
    test_output("tests/depth.txt", "tests/expected/depth.txt");
}

#[test]
fn vcf_to_txt() {
    assert!(Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt vcf-to-txt --genotypes --fmt S --info T X SOMATIC < tests/test.vcf > tests/variant-table.txt")
            .spawn().unwrap().wait().unwrap().success());
    test_output(
        "tests/variant-table.txt",
        "tests/expected/variant-table.txt",
    );
}

#[test]
fn vcf_match() {
    assert!(Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt vcf-match -d 50 -l 20 tests/test3.vcf < tests/test2.vcf > tests/matching.bcf")
            .spawn().unwrap().wait().unwrap().success());
    test_output("tests/matching.bcf", "tests/expected/matching.bcf");
}

#[test]
fn vcf_match_same() {
    assert!(Command::new("bash").arg("-c")
                                .arg("target/debug/rbt vcf-match -d 50 -l 20 tests/test4.vcf < tests/test4.vcf > tests/matching-same.bcf")
                                .spawn().unwrap().wait().unwrap().success());
    test_output(
        "tests/matching-same.bcf",
        "tests/expected/matching-same.bcf",
    );
}

#[test]
fn vcf_fix_iupac_alleles() {
    assert!(Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt vcf-fix-iupac-alleles < tests/test-iupac.vcf > tests/iupac-fixed.bcf")
            .spawn().unwrap().wait().unwrap().success());
    test_output("tests/iupac-fixed.bcf", "tests/expected/iupac-fixed.bcf");
}

#[test]
fn vcf_baf() {
    assert!(Command::new("bash")
        .arg("-c")
        .arg("target/debug/rbt vcf-baf < tests/test-freebayes.vcf > tests/baf.bcf")
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success());
    test_output("tests/baf.bcf", "tests/expected/baf.bcf");
}

#[test]
fn test_call_consensus_reads_two_cluster() {
    assert!(
        Command::new("bash")
                .arg("-c")
                .arg("target/debug/rbt call-consensus-reads fastq --umi-len 3 -u --max-umi-dist 0 --max-seq-dist 2 tests/test-consensus.fastq tests/test-consensus.fastq /tmp/test-consensus.1.fastq /tmp/test-consensus.2.fastq")
            .spawn().unwrap().wait().unwrap().success());
    compare_fastq(
        "/tmp/test-consensus.1.fastq",
        "tests/expected/test-consensus.1.fastq",
    );
    compare_fastq(
        "/tmp/test-consensus.2.fastq",
        "tests/expected/test-consensus.2.fastq",
    );
}

#[test]
fn test_call_consensus_reads_single_cluster() {
    assert!(
        Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt call-consensus-reads fastq --umi-len 3 -u --max-umi-dist 2 --max-seq-dist 2 tests/test-consensus.fastq tests/test-consensus.fastq /tmp/test-consensus_single.1.fastq /tmp/test-consensus_single.2.fastq")
            .spawn().unwrap().wait().unwrap().success());
    compare_fastq(
        "/tmp/test-consensus_single.1.fastq",
        "tests/expected/test-consensus_single.1.fastq",
    );
    compare_fastq(
        "/tmp/test-consensus_single.2.fastq",
        "tests/expected/test-consensus_single.2.fastq",
    );
}

#[test]
fn test_call_overlapping_consensus_reads() {
    assert!(
        Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt call-consensus-reads fastq --umi-len 10 --max-umi-dist 0 --max-seq-dist 8 --insert-size 450 --std-dev 50  tests/overlapping-consensus.1.fastq tests/overlapping-consensus.2.fastq /tmp/test_overlapping-consensus.1.fastq /tmp/test_overlapping-consensus.2.fastq /tmp/test_overlapping-consensus.3.fastq")
            .spawn().unwrap().wait().unwrap().success());
    compare_fastq(
        "/tmp/test_overlapping-consensus.1.fastq",
        "tests/expected/test_overlapping-consensus.1.fastq",
    );
    compare_fastq(
        "/tmp/test_overlapping-consensus.2.fastq",
        "tests/expected/test_overlapping-consensus.2.fastq",
    );
    compare_fastq(
        "/tmp/test_overlapping-consensus.3.fastq",
        "tests/expected/test_overlapping-consensus.3.fastq",
    );
}

#[test]
fn test_call_consensus_from_bam() {
    assert!(
    Command::new("bash")
        .arg("-c")
        .arg("target/debug/rbt call-consensus-reads bam --max-seq-dist 8 tests/overlapping_consensus_marked.bam /tmp/overlapping_consensus_marked.bam")
        .spawn().unwrap().wait().unwrap().success());
    compare_bam(
        "/tmp/overlapping_consensus_marked.bam",
        "tests/expected/overlapping_consensus_marked.bam",
    );
}

#[test]
fn test_vcf_annotate_dgidb() {
    assert!(
        Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt vcf-annotate-dgidb tests/annotate_dgidb_test.vcf > /tmp/annotate_dgidb_test.vcf")
            .spawn().unwrap().wait().unwrap().success());
    test_output(
        "/tmp/annotate_dgidb_test.vcf",
        "tests/expected/annotate_dgidb_test.vcf",
    );
}
