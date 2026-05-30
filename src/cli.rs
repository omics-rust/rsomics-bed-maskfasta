use std::io;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_maskfasta::{MaskMode, maskfasta};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-maskfasta", disable_help_flag = true)]
pub struct Cli {
    /// Reference FASTA to mask (bedtools -fi)
    #[arg(short = 'f', long, value_name = "FASTA")]
    fasta: PathBuf,
    /// BED intervals to mask (bedtools -bed)
    #[arg(short = 'b', long)]
    bed: PathBuf,
    /// Output FASTA (default: stdout; bedtools -fo)
    #[arg(short = 'o', long)]
    output: Option<PathBuf>,
    /// Soft-mask: lowercase instead of replacing with N (bedtools -soft)
    #[arg(long = "soft")]
    soft: bool,
    /// Mask with this character instead of N (bedtools -mc)
    #[arg(long = "mc", value_name = "CHAR")]
    mask_char: Option<char>,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mode = if self.soft {
            MaskMode::Soft
        } else {
            MaskMode::Hard(self.mask_char.map_or(b'N', |c| c as u8))
        };

        let mut stdout_lock;
        let mut file_out;
        let out: &mut dyn io::Write = if let Some(ref p) = self.output {
            file_out = std::fs::File::create(p).map_err(RsomicsError::Io)?;
            &mut file_out
        } else {
            stdout_lock = io::stdout().lock();
            &mut stdout_lock
        };
        maskfasta(&self.fasta, &self.bed, &mode, out)
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Mask FASTA bases overlapping BED intervals (bedtools maskfasta equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["-f <FASTA> -b <BED> [OPTIONS]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('f'),
                long: "fasta",
                aliases: &[],
                value: Some("<FASTA>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Reference FASTA to mask (bedtools -fi)",
                why_default: None,
            },
            FlagSpec {
                short: Some('b'),
                long: "bed",
                aliases: &[],
                value: Some("<BED>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "BED intervals to mask (bedtools -bed)",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdout"),
                description: "Output FASTA path (bedtools -fo)",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "soft",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Soft-mask: lowercase instead of N (bedtools -soft)",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "mc",
                aliases: &[],
                value: Some("<CHAR>"),
                type_hint: Some("char"),
                required: false,
                default: None,
                description: "Mask with this character instead of N (bedtools -mc)",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Hard-mask repeat regions to N",
        command: "rsomics-bed-maskfasta -f genome.fa -b repeats.bed -o masked.fa",
    }],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
