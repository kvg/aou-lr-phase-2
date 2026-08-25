/// Source-specific FORMAT/INFO field propagation rules.
///
/// Destination IDs keep the source names (no prefixes). Genotypes are not
/// propagated. Sparse per-site FORMAT emits only fields from matched sources;
/// SNV/indel and SV share `GQ`/`DP`/`AD` but annotate disjoint site sets.
#[derive(Debug, Clone, Copy)]
pub struct SourceSpec {
    pub label: &'static str,
    /// (source FORMAT tag, output FORMAT tag) — usually identical.
    pub format_fields: &'static [(&'static str, &'static str)],
    pub info_fields: &'static [(&'static str, &'static str)],
}

pub const SNV_INDEL: SourceSpec = SourceSpec {
    label: "snv_indel",
    format_fields: &[("GQ", "GQ"), ("DP", "DP"), ("AD", "AD")],
    info_fields: &[("AF", "AF"), ("AQ", "AQ")],
};

pub const SV: SourceSpec = SourceSpec {
    label: "sv",
    format_fields: &[("GQ", "GQ"), ("DP", "DP"), ("AD", "AD")],
    info_fields: &[
        ("SVTYPE", "SVTYPE"),
        ("SVLEN", "SVLEN"),
        ("SCORE", "SCORE"),
        ("SQ", "SQ"),
        ("GQ", "GQ"),
        ("SUPP_PAV", "SUPP_PAV"),
        ("SUPP_PBSV", "SUPP_PBSV"),
        ("SUPP_SNIFFLES", "SUPP_SNIFFLES"),
    ],
};

pub const FLARE: SourceSpec = SourceSpec {
    label: "flare",
    format_fields: &[
        ("AN1", "AN1"),
        ("AN2", "AN2"),
        ("ANP1", "ANP1"),
        ("ANP2", "ANP2"),
    ],
    info_fields: &[],
};

pub fn all_output_format_ids() -> Vec<&'static str> {
    unique_format_ids([SNV_INDEL, SV, FLARE].into_iter())
}

/// FORMAT IDs contributed by the given source labels (order preserved, deduped).
pub fn format_ids_for_labels(labels: &[&str]) -> Vec<&'static str> {
    unique_format_ids(
        [SNV_INDEL, SV, FLARE]
            .into_iter()
            .filter(|spec| labels.iter().any(|l| *l == spec.label)),
    )
}

fn unique_format_ids<'a>(
    specs: impl Iterator<Item = SourceSpec>,
) -> Vec<&'static str> {
    let mut out = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for spec in specs {
        for &(_, dst) in spec.format_fields {
            if seen.insert(dst) {
                out.push(dst);
            }
        }
    }
    out
}

pub fn format_description(id: &str) -> &'static str {
    match id {
        "GQ" => "Genotype quality",
        "DP" => "Read depth",
        "AD" => "Allele depths",
        "AN1" => "FLARE ancestry haplotype 1",
        "AN2" => "FLARE ancestry haplotype 2",
        "ANP1" => "FLARE ancestry probabilities haplotype 1",
        "ANP2" => "FLARE ancestry probabilities haplotype 2",
        _ => "Propagated annotation",
    }
}

pub fn format_number(id: &str) -> &'static str {
    match id {
        "AD" | "ANP1" | "ANP2" => ".",
        _ => "1",
    }
}

pub fn format_type(id: &str) -> &'static str {
    match id {
        "AD" | "ANP1" | "ANP2" => "String",
        _ => "Integer",
    }
}
