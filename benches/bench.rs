use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bed_maskfasta(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bed-maskfasta");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fasta = manifest.join("tests/golden/genome.fa");
    let bed = manifest.join("tests/golden/regions.bed");
    c.bench_function("rsomics-bed-maskfasta golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args(["-f", fasta.to_str().unwrap(), "-b", bed.to_str().unwrap()])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_bed_maskfasta);
criterion_main!(benches);
