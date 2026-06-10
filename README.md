# svcf

> [!WARNING]
> This package is in early development and API is not yet stable 

Converts a single sample somatic structural variant VCF file (from gridss/purple) to other file formats. 

## Installation

Binaries and installation scripts are available in [releases](https://github.com/selkamand/svcf/releases/)

Or to compile development version from source run:

```
cargo install --git https://github.com/selkamand/svcf
```

## Target Formats 

### BEDPE

> [!NOTE]
> `svcf bedpe` Outputs ONLY breakpoints where both breakends are described (and have FILTER=PASS). 
> Single breakends are excluded.

```
svcf -i <vcf> --from purple --to bedpe
```

Outputs BEDPE format data with one row per breakpoint including the following columns: 

1. chrom1: Chromosome of one side of first breakend in pair.
2. start1: Zero-based starting position of the lower confidence interval of first breakend.
3. end1: One-based end position of the upper confidence interval of first breakend.
4. chrom2: Chromosome of second breakend in pair.
5. start2: Zero-based starting position of the lower confidence interval of second breakend in pair.
6. end2: One-based end position of the upper confidence interval of second breakend in pair.
7. name: Breakpoint identifier.
8. score: quality score (from first breakpoint in svcf)
9. strand1: strand of the first breakend in pair
10. strand2: strand for the second breakend in pair

Plus additional columns: (downstream tools like bedtools allow any number of additional columns - these will just be passed-through)

11. vaf1: purity adjusted VAF of first breakend in pair (e.g. from PURPLE_VAF info field if `--from purple`)
12. vaf2: purity adjusted VAF of second breakend in pair (e.g. from PURPLE_VAF info field if `--from purple`)

> [!Warning]
> A header row is included, in contrast to the official BEDPE specification. We'll leave up to user to drop before plugging into bedtools / other tools

#### Conversion Logic 

The original SV vcf will include 1 row per breakend with classic chrom:pos ref > alt 

They also have an ID col (std vcf) AND a MATEID info column that tells you the identifier of the mate.
Single breakends have MATEID=. and an '.' in their ALT allele (e.g. ACATCG.) instead of '\[' and '\]'
Each SV also has a "CIPOS" INFO representing a Confidence interval around POS for imprecise variants. 
Imagine the SV hits a region of microhomology - it might be very hard to tell exactly which position breaks. 
This CIPOS (confidence interval) is the reason each breakend has its own start & end coord. 

Note gridds assigns a QUAL score per breakpoint (so each side of SV should have the same QUAL score - we just use the first one encountered)

> [!NOTE]
> Which breakend is first / second is determined entirely by the order of the VCF. Presort by coordinate if you want to guarentee lower coord breakends are chrom1/start1/end1 etc. 


### Breakend TSV

> [!WARNING]
> Breakend TSV conversion not implemented yet

```
svcf -i <vcf> --from purple --to breakend-tsv
```

Outputs a TSV with one row per breakend (each side paired breakpoints will have their own row). 

Columns include: 

1. chrom: Chromosome of breakend.
2. position: 1-based position of breakend as described by POS column in vcf.
3. vaf: Purity adjusted variant allele frequency supporting breakend (e.g. from PURPLE_VAF info field if `--from purple`).
4. id: id of breakend.
5. mateid: id of mate (set to `.` if single breakend)
6. qual: quality of breakend.
