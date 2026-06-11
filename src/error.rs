use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub enum Error {
    #[error("failed to open or read VCF: {path}")]
    ReadVcf {
        path: String,
        #[source]
        source: std::io::Error,
    },

    #[error("failed to parse VCF header: {path}")]
    ParseVcfHeader {
        path: String,
        #[source]
        source: std::io::Error,
    },

    #[error("failed to parse VCF record: {path}")]
    ParseVcfRecord {
        path: String,
        #[source]
        source: std::io::Error,
    },

    #[error("failed to evaluate FILTER/PASS status for VCF record {record_id}")]
    FilterStatus {
        record_id: String,
        #[source]
        source: std::io::Error,
    },

    #[error("failed to write BEDPE TSV")]
    WriteBedpe(#[source] csv::Error),

    #[error("failed to flush BEDPE TSV writer")]
    FlushBedpe(#[source] std::io::Error),

    #[error("invalid breakend pairing: {0}")]
    InvalidPairing(String),

    #[error("invalid VCF record: {message}")]
    InvalidRecord { message: String },

    #[error("invalid VCF record {variant}: {message}")]
    InvalidVariantRecord { variant: String, message: String },

    #[error("invalid INFO/{field}: {message}")]
    InvalidInfo { field: String, message: String },

    #[error("invalid ALT allele {alt:?}: failed to infer breakend strand")]
    InvalidAlt { alt: String },

    #[error("integer conversion failed")]
    IntConversion(#[from] std::num::TryFromIntError),
}

impl Error {
    pub(crate) fn invalid_record(message: impl Into<String>) -> Self {
        Self::InvalidRecord {
            message: message.into(),
        }
    }

    pub(crate) fn invalid_variant_record(
        variant: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self::InvalidVariantRecord {
            variant: variant.into(),
            message: message.into(),
        }
    }

    pub(crate) fn invalid_info(field: impl Into<String>, message: impl Into<String>) -> Self {
        Self::InvalidInfo {
            field: field.into(),
            message: message.into(),
        }
    }
}
