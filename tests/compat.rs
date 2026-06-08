use std::io::Write;
use std::path::Path;
use std::process::Command;

use rsomics_bed_maskfasta::{MaskMode, maskfasta};
use tempfile::NamedTempFile;

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn bedtools_present() -> bool {
    Command::new("bedtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
}

// bedtools maskfasta writes to a -fo file; ours writes to stdout. Compare contents.
fn run_bedtools(fasta: &Path, bed: &Path, extra: &[&str]) -> String {
    let out = NamedTempFile::new().unwrap();
    let status = Command::new("bedtools")
        .args(["maskfasta", "-fi"])
        .arg(fasta)
        .arg("-bed")
        .arg(bed)
        .arg("-fo")
        .arg(out.path())
        .args(extra)
        .status()
        .expect("run bedtools");
    assert!(status.success());
    std::fs::read_to_string(out.path()).unwrap()
}

fn run_ours(fasta: &Path, bed: &Path, mode: &MaskMode) -> String {
    let mut buf = Vec::new();
    maskfasta(fasta, bed, mode, &mut buf).unwrap();
    String::from_utf8(buf).unwrap()
}

// Byte-for-byte against output frozen from bedtools v2.31.1
// (`bedtools maskfasta` hard / -soft / -mc X). Always runs so CI guards masking
// and FASTA line-wrap drift even where bedtools is absent.
#[test]
fn matches_bedtools_golden() {
    let (fa, bed) = (golden("genome.fa"), golden("regions.bed"));
    for (mode, expected) in [
        (MaskMode::Hard(b'N'), "mask_hard.expected"),
        (MaskMode::Soft, "mask_soft.expected"),
        (MaskMode::Hard(b'X'), "mask_mc.expected"),
    ] {
        let want = std::fs::read_to_string(golden(expected)).unwrap();
        assert_eq!(run_ours(&fa, &bed, &mode), want, "{expected}");
    }
}

#[test]
fn matches_bedtools_hard_mask() {
    if !bedtools_present() {
        eprintln!("SKIP: bedtools not on PATH");
        return;
    }
    let (fa, bed) = (golden("genome.fa"), golden("regions.bed"));
    assert_eq!(
        run_bedtools(&fa, &bed, &[]),
        run_ours(&fa, &bed, &MaskMode::Hard(b'N'))
    );
}

#[test]
fn matches_bedtools_soft_mask() {
    if !bedtools_present() {
        eprintln!("SKIP: bedtools not on PATH");
        return;
    }
    let (fa, bed) = (golden("genome.fa"), golden("regions.bed"));
    assert_eq!(
        run_bedtools(&fa, &bed, &["-soft"]),
        run_ours(&fa, &bed, &MaskMode::Soft)
    );
}

#[test]
fn matches_bedtools_mask_char() {
    if !bedtools_present() {
        eprintln!("SKIP: bedtools not on PATH");
        return;
    }
    let (fa, bed) = (golden("genome.fa"), golden("regions.bed"));
    assert_eq!(
        run_bedtools(&fa, &bed, &["-mc", "X"]),
        run_ours(&fa, &bed, &MaskMode::Hard(b'X'))
    );
}

#[test]
fn hard_mask_replaces_with_n() {
    let mut fa = NamedTempFile::new().unwrap();
    writeln!(fa, ">c\nACGTACGT").unwrap();
    let mut bed = NamedTempFile::new().unwrap();
    writeln!(bed, "c\t2\t5").unwrap();
    let got = run_ours(fa.path(), bed.path(), &MaskMode::Hard(b'N'));
    assert_eq!(got, ">c\nACNNNCGT\n");
}
