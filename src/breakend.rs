use anyhow::Context;
use anyhow::Result;
use anyhow::bail;
use noodles::vcf;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record::Filters;
use noodles::vcf::variant::record::Ids;
use noodles::vcf::variant::record::info::field;
use noodles::vcf::variant::record::info::field::value::Array;
use std::fs::File;
use std::io::BufReader;
use std::io::Write;
use std::{
    collections::{HashMap, HashSet},
    fmt,
};

// Note start and pos are 0-based & inclusive while end position is non-inclusive (e.g. a region 1000)
pub struct Breakend {
    chrom: String,
    start: u64,
    end: u64,
    pos: u64,
    id: String,
    mateid: Option<String>,
    strand: Strand,
    qual: f32,
    vaf: f32,
}
impl std::fmt::Display for Breakend {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}-{} (pos={}) id={} mateid={:?} strand={} qual={} vaf={}",
            self.chrom,
            self.start,
            self.end,
            self.pos,
            self.id,
            self.mateid,
            self.strand,
            self.qual,
            self.vaf,
        )
    }
}
pub struct Breakpoint {
    first: Breakend,
    second: Breakend,
}

impl Breakpoint {
    fn write_bedpe_record<W: Write>(&self, writer: &mut csv::Writer<W>) -> Result<()> {
        writer.write_record([
            self.first.chrom.as_str(),
            &self.first.start.to_string(),
            &self.first.end.to_string(),
            self.second.chrom.as_str(),
            &self.second.start.to_string(),
            &self.second.end.to_string(),
            &self.id(),
            &self.qual().to_string(),
            &self.first.strand.to_string(),
            &self.second.strand.to_string(),
            &self.first.vaf.to_string(),
            &self.second.vaf.to_string(),
        ])?;

        Ok(())
    }

    fn id(&self) -> String {
        format!("{}.{}", self.first.id, self.second.id)
    }

    fn qual(&self) -> f32 {
        self.first.qual
    }
}

enum Strand {
    Plus,
    Minus,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
        }
    }
}

impl fmt::Debug for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
        }
    }
}

#[derive(Default)]
pub struct BreakpointPairer {
    pub pending_by_id: HashMap<String, Breakend>,
    pub seen_ids: HashSet<String>,
    pub breakpoints: Vec<Breakpoint>,
    pub single_breakends: Vec<Breakend>,
}

// All structural variants
// Breakpoints = Paired
// Single Breakends = No MateID
// Unmatched = had a MateID specified but mate not found in VCF
pub struct StructuralVariants {
    pub breakpoints: Vec<Breakpoint>,
    pub single_breakends: Vec<Breakend>,
    pub unmatched_breakends: Vec<Breakend>,
}

impl BreakpointPairer {
    pub fn push(&mut self, breakend: Breakend) -> Result<()> {
        let id = breakend.id.clone();

        if !self.seen_ids.insert(id.clone()) {
            bail!("duplicate breakend ID found in VCF: {id}");
        }

        let Some(mate_id) = breakend.mateid.clone() else {
            self.single_breakends.push(breakend);
            return Ok(());
        };

        if mate_id == id {
            bail!("breakend {id} has itself as MATEID");
        }

        if let Some(mate) = self.pending_by_id.remove(&mate_id) {
            validate_reciprocal_mates(&mate, &breakend)?;

            self.breakpoints.push(Breakpoint {
                first: mate,
                second: breakend,
            });
        } else {
            self.pending_by_id.insert(id, breakend);
        }

        Ok(())
    }

    pub fn finish(self) -> StructuralVariants {
        StructuralVariants {
            breakpoints: self.breakpoints,
            single_breakends: self.single_breakends,
            unmatched_breakends: self.pending_by_id.into_values().collect(),
        }
    }
}

fn validate_reciprocal_mates(first: &Breakend, second: &Breakend) -> Result<()> {
    if first.mateid.as_deref() != Some(second.id.as_str()) {
        bail!(
            "non-reciprocal MATEID: breakend {} has mateid {:?}, but expected {}",
            first.id,
            first.mateid,
            second.id
        );
    }

    if second.mateid.as_deref() != Some(first.id.as_str()) {
        bail!(
            "non-reciprocal MATEID: breakend {} has mateid {:?}, but expected {}",
            second.id,
            second.mateid,
            first.id
        );
    }

    Ok(())
}

impl StructuralVariants {
    pub fn write_bedpe_tsv<W: Write>(&self, writer: W) -> Result<()> {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(writer);

        writer.write_record([
            "chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1",
            "strand2", "vaf1", "vaf2",
        ])?;

        for breakpoint in &self.breakpoints {
            breakpoint.write_bedpe_record(&mut writer)?;
        }

        writer.flush()?;

        Ok(())
    }

    pub fn n_breakpoints(&self) -> usize {
        self.breakpoints.len()
    }

    pub fn n_single_breakends(&self) -> usize {
        self.single_breakends.len()
    }

    pub fn n_unmatched_breakends(&self) -> usize {
        self.unmatched_breakends.len()
    }
}

pub fn vcf_to_structural_variants(vcf: &str, vaf_field: &str) -> Result<StructuralVariants> {
    let mut reader = File::open(vcf)
        .map(BufReader::new)
        .map(noodles::vcf::io::Reader::new)?;
    let header = reader.read_header()?;

    //let mut breakpoints = Vec::<Breakpoint>::new();
    let mut pairer = BreakpointPairer::default();

    // Iterate through VCF record
    for result in reader.records() {
        let record = result?; // Skip non-pass records

        // Skip non-pass variants
        if !record.filters().is_pass(&header)? {
            continue;
        }

        let current_breakend = record_to_breakend(&record, &header, vaf_field)?;

        // Add to 'BreakpointPairer'
        pairer.push(current_breakend)?;
    }

    let pairing = pairer.finish();
    Ok(pairing)
}

pub fn record_to_breakend(
    record: &vcf::Record,
    header: &vcf::Header,
    vaf_field: &str,
) -> Result<Breakend> {
    // Grab ID
    let breakend_id = parse_id(record)?;

    // Fetch MATEID (if none: Single breakend
    let mateid = parse_mate_id(record, header).context(format!("Variant ID: {breakend_id}"))?;

    // Grab Position
    let Some(pos_res) = record.variant_start() else {
        bail! {
            "Failed to get position in row for variant: {breakend_id:#?}",
        }
    };
    let pos_usize = pos_res?.get(); // 1-based position
    let pos_1based: u64 = pos_usize.try_into()?; // Convert to u64
    let pos: u64 = pos_1based.saturating_sub(1); // Convert to 0-based position

    // Grab position confidence interval CIPOS or if anything goes wrong set as 0,0
    let cipos = parse_cipos(record, header).unwrap_or_default(); // default will
    // be (0_i32, 0_i32). NOTE in PURPLE VCFs CIPOS lower element will be negative

    // Get the absolute value of lower and higher end
    let cipos_low: u64 = cipos.0.unsigned_abs().into();
    let cipos_high: u64 = cipos.1.unsigned_abs().into();

    // Get Start Position by subtracting absolute value of lower CIPOS to 0-based position
    let start = pos.saturating_sub(cipos_low);

    // Get End Position by adding absolute value of upper CIPOS to 1-based position since bedpe
    // end positions are non-inclusive
    let end = pos_1based.saturating_add(cipos_high);

    // Grab VAF
    let vaf = parse_vaf(record, header, vaf_field).context(format!("Variant: {breakend_id}"))?;

    // Grab QUAL  (or NAN if anything went wrong)
    let qual = record
        .quality_score()
        .unwrap_or(Ok(f32::NAN))
        .unwrap_or(f32::NAN);

    // Grab CHROM
    let chrom = record.reference_sequence_name();

    // Grab ALT (used to infer strand)
    let altbases = record.alternate_bases();
    if altbases.len() != 1 {
        bail!(
            "Expected a single alternative sequece for variant {breakend_id} but found {}",
            altbases.len()
        )
    }
    let alt = get_first_alt_as_string(record)
        .context("Failed to pull a valid alternative sequence for variant {breakend_id}")?;

    // Infer strand from alt sequence Grab strand
    let strand = alt_to_strand(alt)?;

    // Create Breakend
    Ok(Breakend {
        chrom: chrom.to_string(),
        start,
        end,
        pos,
        id: breakend_id,
        mateid,
        strand,
        qual,
        vaf,
    })
}

// Extract MATEID info field
fn parse_mate_id(record: &vcf::Record, header: &vcf::Header) -> Result<Option<String>> {
    let info = record.info();
    let mateid_key: &str = "MATEID";
    let Some(value_result) = info.get(header, "MATEID") else {
        return Ok(None);
    };

    let Some(value) = value_result? else {
        return Ok(None);
    };

    let mate_id = match value {
        field::Value::Integer(n) => n.to_string(),
        field::Value::Float(n) => n.to_string(),
        field::Value::Character(c) => c.to_string(),
        field::Value::String(s) => s.to_string(),
        field::Value::Flag => {
            anyhow::bail!("MATEID field is a flag type, which can not be coerced to a string")
        }
        field::Value::Array(arr) => {
            parse_one_string_from_array(arr, mateid_key).context("Failed to extract MATEID")?
        }
    };

    if mate_id == "." {
        Ok(None)
    } else {
        Ok(Some(mate_id))
    }
}

// Get stock standard ID column from VCF
fn parse_id(record: &vcf::Record) -> Result<String> {
    let ids = record.ids();
    if ids.len() > 1 {
        bail!(
            "Multiple IDs found in a single ID column of SV vcf. This is unexpected and should be resolved before converting to BEDPE"
        )
    }
    if ids.len() == 0 {
        bail!(
            "Some SV entries lack any value in the ID column. This is unexpected and should be resolved before converting to BEDPE"
        )
    }
    // Get first ID
    let Some(id_str) = ids.iter().next() else {
        bail!(
            "At least one SV entry lacks a value in the ID column. This is unexpected and should be resolved before converting to BEDPE"
        )
    };

    let id = id_str.to_owned();

    if id == "." {
        bail!(
            "At least one SV entry lacks a value in the ID column. This is unexpected and should be resolved before converting to BEDPE"
        )
    } else {
        Ok(id)
    }
}

fn parse_cipos(record: &vcf::Record, header: &vcf::Header) -> Result<(i32, i32)> {
    let info = record.info();

    let Some(value_result) = info.get(header, "CIPOS") else {
        bail!("Can not find CIPOS field")
    };

    let Some(value) = value_result? else {
        bail!("Can not find CIPOS field")
    };

    let array = match value {
        field::Value::Array(array) => array,
        _ => bail!("INFO/CIPOS expected an array, got {value:#?}"),
    };

    let array_int = match array {
        field::value::Array::Integer(values) => values,
        _ => bail!("INFO/CIPOS expected to be an integer array, got {array:#?}"),
    };

    if array_int.len() != 2 {
        bail!(
            "INFO/CIPOS should contain 2 numbers per entry, found {}",
            array_int.len()
        )
    }

    // Grab the first two elements of the array
    let mut iter = array_int.iter();

    let lo_option = iter
        .next()
        .transpose()?
        .context("INFO/CIPOS expected first integer value")?;

    let hi_option = iter
        .next()
        .transpose()?
        .context("INFO/CIPOS expected second integer value")?;

    let Some(lo) = lo_option else {
        bail!("INFO/CIPOS expected first integer value");
    };
    let Some(hi) = hi_option else {
        bail!("INFO/CIPOS expected second integer value");
    };

    Ok((lo, hi))
}

fn parse_vaf(record: &vcf::Record, header: &vcf::Header, vaf_field: &str) -> Result<f32> {
    // Keep this binding, otherwise record.info() may be a temporary that gets dropped.
    let info = record.info();

    let Some(value_result) = info.get(header, vaf_field) else {
        bail!("cannot find INFO/{vaf_field} field");
    };

    let Some(value) = value_result? else {
        bail!("INFO/{vaf_field} is present but has no value");
    };

    match value {
        // This may happen if the header or parser treats it as a scalar.
        field::Value::Float(vaf) => Ok(vaf),

        // This is the likely case for Number=.
        field::Value::Array(array) => parse_first_float_from_array(array, vaf_field),

        field::Value::Integer(_) => {
            bail!("INFO/{vaf_field} expected Float, got Integer")
        }

        field::Value::Flag => {
            bail!("INFO/{vaf_field} expected Float, got Flag")
        }

        field::Value::Character(_) => {
            bail!("INFO/{vaf_field} expected Float, got Character")
        }

        field::Value::String(_) => {
            bail!("INFO/{vaf_field} expected Float, got String")
        }
    }
}

fn parse_first_float_from_array(array: Array<'_>, field_name: &str) -> Result<f32> {
    let Array::Float(values) = array else {
        bail!("INFO/{field_name} expected a Float array");
    };

    // if values.len() != 1 {
    //     bail!(
    //         "INFO/{field_name} expected exactly one float value, got {}",
    //         values.len()
    //     );
    // }

    let mut iter = values.iter();

    let vaf = iter
        .next()
        .context(format!("INFO/{field_name} expected one float value"))??
        .context(format!("INFO/{field_name} contains a missing float value"))?;

    Ok(vaf)
}

fn parse_one_string_from_array(array: Array<'_>, field_name: &str) -> Result<String> {
    let Array::String(values) = array else {
        bail!("INFO/{field_name} expected a String array");
    };

    if values.len() != 1 {
        bail!(
            "INFO/{field_name} expected exactly one string value, got {}",
            values.len()
        );
    }

    let mut iter = values.iter();

    let vaf = iter
        .next()
        .context(format!("INFO/{field_name} expected one string value"))??
        .context(format!("INFO/{field_name} contains a missing string value"))?;

    Ok(vaf.to_string())
}

// ALT pattern             local breakend strand
// ------------------------------------------------
// s[chr:pos[              +
// s]chr:pos]              +
// ]chr:pos]s              -
// [chr:pos[s              -
// s.                      +
// .s                      -
fn alt_to_strand(alt: String) -> Result<Strand> {
    if alt.ends_with('[') {
        return Ok(Strand::Plus);
    }

    if alt.ends_with(']') {
        return Ok(Strand::Plus);
    }

    if alt.ends_with('.') {
        return Ok(Strand::Plus);
    }

    if alt.starts_with(']') {
        return Ok(Strand::Minus);
    }

    if alt.starts_with('[') {
        return Ok(Strand::Minus);
    }

    if alt.starts_with('.') {
        return Ok(Strand::Minus);
    }

    bail!("Failed to infer strand from ALT sequence: {alt}")
}

// Get first alternate base as string. If anything goes wrong or alt is empty return None
fn get_first_alt_as_string(record: &vcf::Record) -> Option<String> {
    let alternative_bases = record.alternate_bases();
    let first = alternative_bases.iter().next()?;

    let unwrapped = match first {
        Ok(b) => b,
        Err(_) => return None,
    };

    Some(unwrapped.to_string())
}

pub fn write_bedpe_tsv<W: Write>(writer: W, breakpoints: &[Breakpoint]) -> Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);

    writer.write_record([
        "chrom1", "start1", "end1", "chrom2", "start2", "end2", "id", "qual", "strand1", "strand2",
    ])?;

    for breakpoint in breakpoints {
        //let svclass = infer_svclass(breakpoint);

        writer.write_record([
            breakpoint.first.chrom.as_str(),
            &breakpoint.first.start.to_string(),
            &breakpoint.first.end.to_string(),
            breakpoint.second.chrom.as_str(),
            &breakpoint.second.start.to_string(),
            &breakpoint.second.end.to_string(),
            &format!(
                "{}.{}",
                breakpoint.first.id.as_str(),
                breakpoint.second.id.as_str()
            ),
            &breakpoint.first.qual.to_string(),
            &breakpoint.first.strand.to_string(),
            &breakpoint.second.strand.to_string(),
            // svclass,
        ])?;
    }

    writer.flush()?;

    Ok(())
}
