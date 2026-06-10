use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use std::{fmt, path::PathBuf};
use svcf::breakend::vcf_to_structural_variants;

#[derive(Debug, Clone, Copy, ValueEnum)]
enum VcfTypes {
    Purple,
}

impl fmt::Display for VcfTypes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VcfTypes::Purple => write!(f, "Purple"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, ValueEnum)]
enum OutputTypes {
    Bedpe,
    BreakendTsv,
}

#[derive(Parser)]
#[command(
    version,
    about = "Transform structural variant VCFs to other formats",
    long_about = None
)]
struct Cli {
    /// Path to a structural variant vcf
    #[arg(short = 'i', long = "input", value_name = "vcf")]
    vcf: PathBuf,

    /// Input sv vcf filetype
    #[arg(long, value_enum, default_value_t = VcfTypes::Purple)]
    from: VcfTypes,

    /// Output filetype
    #[arg(long, value_enum, value_name = "filetype")]
    to: OutputTypes,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if !cli.vcf.exists() {
        anyhow::bail!("Failed to find VCF file [{}]", cli.vcf.display())
    }

    // CONFIG
    // Read svVCF
    let vcf = cli
        .vcf
        .to_str()
        .context("Failed to convert vcf path to str")?;

    // Define VAF field based on input tool
    let vaf_field = match cli.from {
        VcfTypes::Purple => "PURPLE_AF",
    };

    match cli.to {
        OutputTypes::Bedpe => vcf_to_bedpe(vcf, vaf_field)?,
        OutputTypes::BreakendTsv => todo!("BreakendTsv conversion is not yet implemented"),
    };

    Ok(())
}

fn vcf_to_bedpe(vcf: &str, vaf_field: &str) -> Result<()> {
    // Get serialised version fo sv VCF
    let structural_variants = vcf_to_structural_variants(vcf, vaf_field)?;

    // Print to STDOUT
    let stdout = std::io::stdout();
    structural_variants.write_bedpe_tsv(&stdout)?;

    Ok(())
}
