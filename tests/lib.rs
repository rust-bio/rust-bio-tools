use std::process::Command;
use std::fs;


fn test_output(result: &str, expected: &str) {
    assert!(Command::new("cmp")
            .arg(result)
            .arg(expected)
            .spawn().unwrap().wait().unwrap().success());
    fs::remove_file(result).unwrap();
}


#[test]
fn fastq_split() {
    assert!(Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt fastq-split tests/A.fastq tests/B.fastq < tests/test.fastq")
            .spawn().unwrap().wait().unwrap().success());
    test_output("tests/A.fastq", "tests/expected/A.fastq");
    test_output("tests/B.fastq", "tests/expected/B.fastq");
}

#[test]
fn fastq_filter() {
    assert!(Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt fastq-filter tests/ids.txt < tests/test.fastq > tests/filtered.fastq")
            .spawn().unwrap().wait().unwrap().success());
    test_output("tests/filtered.fastq", "tests/expected/B.fastq");
}

#[test]
fn bam_depth() {
    assert!(Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt bam-depth tests/test.bam < tests/pos.txt > tests/depth.txt")
            .spawn().unwrap().wait().unwrap().success());
    test_output("tests/depth.txt", "tests/expected/depth.txt");
}


#[test]
fn vcf_to_txt() {
    assert!(Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt vcf-to-txt --genotypes --fmt S --info T X SOMATIC < tests/test.vcf > tests/variant-table.txt")
            .spawn().unwrap().wait().unwrap().success());
    test_output("tests/variant-table.txt", "tests/expected/variant-table.txt");
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
    test_output("tests/matching-same.bcf", "tests/expected/matching-same.bcf");
}

#[test]
fn vcf_baf() {
    assert!(Command::new("bash")
            .arg("-c")
            .arg("target/debug/rbt vcf-baf < tests/test-freebayes.vcf > tests/baf.bcf")
            .spawn().unwrap().wait().unwrap().success());
    test_output("tests/baf.bcf", "tests/expected/baf.bcf");
}

#[test]
fn test_group_by_umi() {
    assert!(
        Command::new("bash")
                .arg("-c")
                .arg("cargo build --release && target/release/rbt group-by-umi tests/test-group-umi.bam > tests/test-grouped-umi.bam")
                .spawn().unwrap().wait().unwrap().success()
    );
}
