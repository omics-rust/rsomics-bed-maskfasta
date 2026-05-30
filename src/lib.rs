//! Mask FASTA bases overlapping BED intervals — `bedtools maskfasta`.
//!
//! Default replaces masked bases with `N`; `Soft` lowercases them; `Hard(c)`
//! uses a chosen character. The original per-sequence line width (from each
//! record's first data line) is preserved, so output is byte-identical to
//! bedtools. BED intervals are merged per chromosome and the FASTA is streamed
//! once against them — O(N + M log M) for N bases, M intervals.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub enum MaskMode {
    /// Replace with `N` or the given character.
    Hard(u8),
    /// Replace with the lowercase equivalent.
    Soft,
}

pub fn maskfasta(
    fasta_path: &Path,
    bed_path: &Path,
    mode: &MaskMode,
    output: &mut dyn Write,
) -> Result<()> {
    let intervals = load_intervals(bed_path)?;
    let out = &mut BufWriter::with_capacity(256 * 1024, output);
    stream_and_mask(fasta_path, &intervals, mode, out)
}

/// Sorted, merged half-open intervals [start, end) for one chromosome.
type Iv = (u64, u64);

fn load_intervals(bed_path: &Path) -> Result<HashMap<String, Vec<Iv>>> {
    let file = File::open(bed_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bed_path.display())))?;
    let mut by_chrom: HashMap<String, Vec<Iv>> = HashMap::new();

    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.splitn(4, '\t');
        let chrom = match fields.next() {
            Some(c) => c.to_string(),
            None => continue,
        };
        let start: u64 = match fields.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let end: u64 = match fields.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        by_chrom.entry(chrom).or_default().push((start, end));
    }

    for ivs in by_chrom.values_mut() {
        ivs.sort_unstable_by_key(|&(s, e)| (s, e));
        let mut merged: Vec<Iv> = Vec::with_capacity(ivs.len());
        for &(s, e) in ivs.iter() {
            if let Some(last) = merged.last_mut()
                && s < last.1
            {
                last.1 = last.1.max(e);
                continue;
            }
            merged.push((s, e));
        }
        *ivs = merged;
    }

    Ok(by_chrom)
}

fn stream_and_mask(
    fasta_path: &Path,
    intervals: &HashMap<String, Vec<Iv>>,
    mode: &MaskMode,
    out: &mut impl Write,
) -> Result<()> {
    let file = File::open(fasta_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", fasta_path.display())))?;

    let mut cur_name: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::with_capacity(1 << 20);
    let mut cur_line_width: usize = 0;

    for line_res in BufReader::new(file).lines() {
        let line = line_res.map_err(RsomicsError::Io)?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = cur_name.take() {
                flush_seq(out, &name, &cur_seq, cur_line_width, intervals, mode)?;
                cur_seq.clear();
                cur_line_width = 0;
            }
            let name = rest.split_whitespace().next().unwrap_or("").to_string();
            cur_name = Some(name);
        } else if cur_name.is_some() {
            let data = line.as_bytes();
            if cur_line_width == 0 && !data.is_empty() {
                cur_line_width = data.len();
            }
            cur_seq.extend_from_slice(data);
        }
    }

    if let Some(name) = cur_name {
        flush_seq(out, &name, &cur_seq, cur_line_width, intervals, mode)?;
    }

    out.flush().map_err(RsomicsError::Io)
}

fn flush_seq(
    out: &mut impl Write,
    name: &str,
    seq: &[u8],
    line_width: usize,
    intervals: &HashMap<String, Vec<Iv>>,
    mode: &MaskMode,
) -> Result<()> {
    writeln!(out, ">{name}").map_err(RsomicsError::Io)?;

    if seq.is_empty() {
        return Ok(());
    }

    let ivs = intervals.get(name).map(|v| v.as_slice()).unwrap_or(&[]);
    let wrap = if line_width == 0 {
        seq.len()
    } else {
        line_width
    };

    let mut flat_pos: u64 = 0;
    let mut iv_idx: usize = 0;

    for chunk in seq.chunks(wrap) {
        for &b in chunk {
            while iv_idx < ivs.len() && ivs[iv_idx].1 <= flat_pos {
                iv_idx += 1;
            }
            let masked = iv_idx < ivs.len() && ivs[iv_idx].0 <= flat_pos;
            let out_byte = if masked {
                match mode {
                    MaskMode::Hard(ch) => *ch,
                    MaskMode::Soft => b.to_ascii_lowercase(),
                }
            } else {
                b
            };
            out.write_all(&[out_byte]).map_err(RsomicsError::Io)?;
            flat_pos += 1;
        }
        writeln!(out).map_err(RsomicsError::Io)?;
    }

    Ok(())
}
