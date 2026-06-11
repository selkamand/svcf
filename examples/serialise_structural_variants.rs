//! Read structural variant VCF file into StructuralVariants object and convert to bedpe
//!
//! Output matches `svcf --input <vcf> --from purple --to bedpe`
//!
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let vcf = std::env::args().nth(1).expect("missing vcf");
    let structural_variants = svcf::breakend::vcf_to_structural_variants(&vcf, "PURPLE_AF")?;

    let stdout = std::io::stdout();
    structural_variants.write_bedpe_tsv(&stdout)?;

    Ok(())
}
