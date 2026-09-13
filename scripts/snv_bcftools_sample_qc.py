#!/usr/bin/env python3
"""Per-chromosome ``bcftools stats`` for the DeepVariant + GLnexus joint callset.

One sequential scan per shard. Parallelize across chromosomes. Resume skips
shards whose ``*.stats.txt`` already looks complete.

Usage:
  python scripts/snv_bcftools_sample_qc.py --from-firecloud --pull-stats --out-dir summaries
  python scripts/snv_bcftools_sample_qc.py --merge-dir work/shards --out-dir summaries
  python scripts/snv_bcftools_sample_qc.py --manifest-tsv GL_INTERVAL_set.tsv --run --jobs 8
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
NUCLEAR = AUTOSOMES | {"chrX", "chrY"}
DROP = {"chrM", "chrMT", "MT"}
STATS_COLUMN = "stats"

PSC_INT_COLS = (
    "n_ref_hom",
    "n_hom_var_snp",
    "n_het_snp",
    "n_ti",
    "n_tv",
    "n_indel",
    "n_singletons",
    "n_hap_ref",
    "n_hap_alt",
    "n_missing",
)


def env_flag(name: str, *, default: bool = False) -> bool:
    value = os.environ.get(name, "")
    if not value:
        return default
    return value.lower() in {"1", "true", "yes", "on"}


def gsutil_prefix() -> list[str]:
    cmd = ["gsutil"]
    project = os.environ.get("GOOGLE_PROJECT", "").strip()
    if project:
        cmd.extend(["-u", project])
    return cmd


def stats_uri_column(df: pd.DataFrame, preferred: str = STATS_COLUMN) -> str:
    if preferred in df.columns:
        return preferred
    for name in ("stats", "stats_txt", "bcftools_stats"):
        if name in df.columns:
            return name
    raise ValueError(
        f"No stats URI column (looked for {preferred!r}); columns={list(df.columns)}"
    )


def select_merge_intervals(
    intervals: list[str],
    *,
    autosomes_only: bool = False,
) -> list[str]:
    from resolve_gl_interval_manifest import natural_chrom_key

    wanted = AUTOSOMES if autosomes_only else NUCLEAR
    dropped = sorted(i for i in intervals if i in DROP)
    keep = [i for i in intervals if i in wanted]
    missing = sorted(wanted - set(keep), key=natural_chrom_key)
    if missing:
        raise ValueError(f"missing required shards: {missing}")
    extra = sorted(set(intervals) - wanted - DROP, key=natural_chrom_key)
    if extra:
        print(f"ignoring extra intervals: {extra}", flush=True)
    if dropped:
        print(f"skipping {dropped} (not in nuclear catalog)", flush=True)
    return sorted(keep, key=natural_chrom_key)


def pull_stats_from_table(
    df: pd.DataFrame,
    stats_dir: Path,
    *,
    id_column: str = "interval_id",
    stats_column: str | None = None,
    jobs: int = 8,
    force: bool = False,
) -> dict[str, Any]:
    """Download table `stats` URIs to `{interval_id}.stats.txt`."""
    stats_dir.mkdir(parents=True, exist_ok=True)
    col = stats_column or stats_uri_column(df)
    if id_column not in df.columns:
        raise ValueError(f"Missing {id_column}; columns={list(df.columns)}")

    jobs_spec: list[tuple[str, str, Path]] = []
    for _, row in df.iterrows():
        interval_id = str(row[id_column]).strip()
        uri = str(row[col]).strip()
        if not interval_id or not uri.startswith("gs://"):
            continue
        dest = stats_dir / f"{interval_id}.stats.txt"
        jobs_spec.append((interval_id, uri, dest))
    if not jobs_spec:
        raise ValueError(f"Column {col} has no gs:// stats URIs ({len(df)} table rows)")

    def one(interval_id: str, uri: str, dest: Path) -> str:
        if dest.is_file() and not force and stats_looks_complete(dest):
            return "skip"
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_suffix(dest.suffix + ".partial")
        subprocess.check_call(gsutil_prefix() + ["cp", uri, str(tmp)])
        tmp.replace(dest)
        if not stats_looks_complete(dest):
            dest.unlink(missing_ok=True)
            raise RuntimeError(f"incomplete stats after download: {interval_id} {uri}")
        return "ok"

    workers = max(1, int(jobs))
    if workers == 1:
        results = [one(*job) for job in jobs_spec]
    else:
        results = []
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(one, *job): job[0] for job in jobs_spec}
            for fut in as_completed(futs):
                interval_id = futs[fut]
                try:
                    results.append(fut.result())
                except Exception as exc:  # noqa: BLE001
                    raise RuntimeError(f"{interval_id} download failed: {exc}") from exc
    n_ok = sum(1 for status in results if status != "skip")
    n_skip = sum(1 for status in results if status == "skip")
    return {
        "stats_uri_column": col,
        "n_table_rows": int(len(df)),
        "n_with_uri": len(jobs_spec),
        "n_downloaded": n_ok,
        "n_skipped": n_skip,
        "stats_dir": str(stats_dir),
    }


def stats_looks_complete(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size < 200:
        return False
    text = path.read_text(errors="replace")
    if "# PSC" not in text or "\nPSC\t" not in text:
        return False
    # A killed pipe often leaves a truncated last line without a trailing newline.
    return text.endswith("\n") and "Broken pipe" not in text


def parse_stats_sections(path: Path) -> dict[str, pd.DataFrame]:
    sn_rows: list[dict[str, Any]] = []
    tstv_rows: list[dict[str, Any]] = []
    idd_rows: list[dict[str, Any]] = []
    psc_rows: list[dict[str, Any]] = []
    with path.open() as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            kind = parts[0]
            if kind == "SN" and len(parts) >= 4:
                sn_rows.append({"key": parts[2].rstrip(":").strip(), "value": int(float(parts[3]))})
            elif kind == "TSTV" and len(parts) >= 5:
                tstv_rows.append({"n_ti": int(parts[2]), "n_tv": int(parts[3])})
            elif kind == "IDD" and len(parts) >= 4:
                idd_rows.append({"length": int(parts[2]), "count": int(parts[3])})
            elif kind == "PSC" and len(parts) >= 9:
                psc_rows.append(
                    {
                        "research_id": parts[2],
                        "n_ref_hom": int(parts[3]),
                        "n_hom_var_snp": int(parts[4]),
                        "n_het_snp": int(parts[5]),
                        "n_ti": int(parts[6]),
                        "n_tv": int(parts[7]),
                        "n_indel": int(parts[8]),
                        "n_singletons": int(parts[10]) if len(parts) > 10 else 0,
                        "n_hap_ref": int(parts[11]) if len(parts) > 11 else 0,
                        "n_hap_alt": int(parts[12]) if len(parts) > 12 else 0,
                        "n_missing": int(parts[13]) if len(parts) > 13 else 0,
                    }
                )
    return {
        "sn": pd.DataFrame(sn_rows),
        "tstv": pd.DataFrame(tstv_rows),
        "idd": pd.DataFrame(idd_rows),
        "psc": pd.DataFrame(psc_rows),
    }


def indel_bin_counts(idd: pd.DataFrame) -> dict[str, int]:
    if idd.empty:
        return {
            "n_ins": 0,
            "n_del": 0,
            "n_ins_lt20": 0,
            "n_del_lt20": 0,
            "n_ins_lt50": 0,
            "n_del_lt50": 0,
        }
    length = idd["length"].astype(int)
    count = idd["count"].astype(int)
    ins = length > 0
    dele = length < 0
    abs_len = length.abs()
    return {
        "n_ins": int(count.loc[ins].sum()),
        "n_del": int(count.loc[dele].sum()),
        "n_ins_lt20": int(count.loc[ins & (abs_len < 20)].sum()),
        "n_del_lt20": int(count.loc[dele & (abs_len < 20)].sum()),
        "n_ins_lt50": int(count.loc[ins & (abs_len < 50)].sum()),
        "n_del_lt50": int(count.loc[dele & (abs_len < 50)].sum()),
    }


def sn_value(sn: pd.DataFrame, *keys: str) -> int:
    if sn.empty:
        return 0
    for key in keys:
        hit = sn.loc[sn["key"] == key, "value"]
        if len(hit):
            return int(hit.iloc[0])
    return 0


def shard_site_row(interval_id: str, parsed: dict[str, pd.DataFrame]) -> dict[str, Any]:
    sn = parsed["sn"]
    tstv = parsed["tstv"]
    n_ti = int(tstv["n_ti"].sum()) if len(tstv) else 0
    n_tv = int(tstv["n_tv"].sum()) if len(tstv) else 0
    row = {
        "interval_id": interval_id,
        "n_records": sn_value(sn, "number of records"),
        "n_snv": sn_value(sn, "number of SNPs"),
        "n_mnp": sn_value(sn, "number of MNPs"),
        "n_indel_sites": sn_value(sn, "number of indels"),
        "n_other": sn_value(sn, "number of others"),
        "n_multiallelic": sn_value(sn, "number of multiallelic sites"),
        "n_multiallelic_snp": sn_value(sn, "number of multiallelic SNP sites"),
        "n_ti": n_ti,
        "n_tv": n_tv,
        **indel_bin_counts(parsed["idd"]),
    }
    return row


def add_ratios(df: pd.DataFrame, prefix: str = "") -> pd.DataFrame:
    p = f"{prefix}_" if prefix else ""
    ti = df[f"{p}n_ti"]
    tv = df[f"{p}n_tv"]
    het = df[f"{p}n_het_snp"]
    hom = df[f"{p}n_hom_var_snp"]
    df[f"{p}n_snp"] = het + hom
    df[f"{p}ti_tv"] = ti / tv.where(tv > 0)
    df[f"{p}het_hom_snp"] = het / hom.where(hom > 0)
    return df


def run_one_shard(
    *,
    interval_id: str,
    vcf_uri: str,
    out_txt: str,
    threads: int,
    pass_only: bool,
    force: bool,
) -> dict[str, Any]:
    dest = Path(out_txt)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file() and not force and stats_looks_complete(dest):
        return {"interval_id": interval_id, "status": "skip", "path": str(dest)}

    bcftools = shutil.which("bcftools")
    if not bcftools:
        raise FileNotFoundError("bcftools is not on PATH")

    cmd = [bcftools, "stats", "-s", "-", "--threads", str(threads)]
    if pass_only:
        # GLnexus writes FILTER="." (unfiltered), not PASS.
        cmd.extend(["-f", "PASS,."])

    tmp = dest.with_suffix(dest.suffix + ".partial")
    if tmp.exists():
        tmp.unlink()

    if vcf_uri.startswith("gs://"):
        cat = gsutil_prefix() + ["cat", vcf_uri]
        with tmp.open("w") as handle:
            producer = subprocess.Popen(cat, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            consumer = subprocess.Popen(cmd, stdin=producer.stdout, stdout=handle, stderr=subprocess.PIPE)
            if producer.stdout is not None:
                producer.stdout.close()
            _out, consumer_err = consumer.communicate()
            _ignored, producer_err = producer.communicate()
        # 141 = SIGPIPE if bcftools exits first; treat as success when stats look complete.
        if producer.returncode not in (0, None, 141) and producer.returncode != 0:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(
                f"gsutil cat failed for {interval_id} ({producer.returncode}): "
                f"{producer_err.decode(errors='replace')[:800]}"
            )
        if consumer.returncode != 0:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(
                f"bcftools stats failed for {interval_id} ({consumer.returncode}): "
                f"{consumer_err.decode(errors='replace')[:800]}"
            )
    else:
        local = Path(vcf_uri)
        if not local.is_file():
            raise FileNotFoundError(vcf_uri)
        with tmp.open("w") as handle:
            proc = subprocess.run(
                [*cmd, str(local)],
                stdout=handle,
                stderr=subprocess.PIPE,
                check=False,
            )
        if proc.returncode != 0:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(
                f"bcftools stats failed for {interval_id} ({proc.returncode}): "
                f"{proc.stderr.decode(errors='replace')[:800]}"
            )

    if not stats_looks_complete(tmp):
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"incomplete bcftools stats output for {interval_id}: {tmp}")
    tmp.replace(dest)
    return {"interval_id": interval_id, "status": "ok", "path": str(dest)}


def _run_one_shard_star(payload: dict[str, Any]) -> dict[str, Any]:
    return run_one_shard(**payload)


def merge_shards(shard_dir: Path, intervals: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    site_rows: list[dict[str, Any]] = []
    psc_parts: list[pd.DataFrame] = []
    for interval_id in intervals:
        path = shard_dir / f"{interval_id}.stats.txt"
        if not stats_looks_complete(path):
            raise FileNotFoundError(f"missing or incomplete stats for {interval_id}: {path}")
        parsed = parse_stats_sections(path)
        site_row = shard_site_row(interval_id, parsed)
        if site_row["n_records"] <= 0:
            raise ValueError(
                f"{interval_id}: 0 records in {path}. GLnexus FILTER is usually '.' not "
                "PASS; rerun BcftoolsGlnexusStats with apply_filters=PASS,."
            )
        site_rows.append(site_row)
        psc = parsed["psc"]
        if psc.empty:
            raise ValueError(f"no PSC rows in {path}; rerun with bcftools stats -s -")
        psc = psc.copy()
        psc["interval_id"] = interval_id
        psc_parts.append(psc)

    sites = pd.DataFrame(site_rows)
    psc = pd.concat(psc_parts, ignore_index=True)

    def sum_sites(mask: pd.Series) -> dict[str, Any]:
        sub = sites.loc[mask]
        out = {
            "n_sites": int(sub["n_records"].sum()),
            "n_snv": int(sub["n_snv"].sum()),
            "n_ins": int(sub["n_ins"].sum()),
            "n_del": int(sub["n_del"].sum()),
            "n_ins_lt20": int(sub["n_ins_lt20"].sum()),
            "n_del_lt20": int(sub["n_del_lt20"].sum()),
            "n_ins_lt50": int(sub["n_ins_lt50"].sum()),
            "n_del_lt50": int(sub["n_del_lt50"].sum()),
            "n_ti": int(sub["n_ti"].sum()),
            "n_tv": int(sub["n_tv"].sum()),
            "n_multiallelic": int(sub["n_multiallelic"].sum()),
        }
        out["ti_tv"] = (out["n_ti"] / out["n_tv"]) if out["n_tv"] else float("nan")
        return out

    auto_mask = sites["interval_id"].isin(AUTOSOMES)
    nuc_mask = sites["interval_id"].isin(NUCLEAR)
    catalog = {"n_intervals": len(sites)}
    for prefix, mask in (("auto", auto_mask), ("nuclear", nuc_mask)):
        for key, value in sum_sites(mask).items():
            catalog[f"{prefix}_{key}"] = value
    catalog_df = pd.DataFrame([catalog])

    def collapse(mask: pd.Series, prefix: str) -> pd.DataFrame:
        sub = psc.loc[psc["interval_id"].isin(set(sites.loc[mask, "interval_id"]))]
        grouped = sub.groupby("research_id", as_index=False)[list(PSC_INT_COLS)].sum()
        grouped = add_ratios(grouped)
        grouped = grouped.rename(columns={c: f"{prefix}_{c}" for c in grouped.columns if c != "research_id"})
        return grouped

    auto = collapse(auto_mask, "auto")
    nuclear = collapse(nuc_mask, "nuclear")
    sample_df = auto.merge(nuclear, on="research_id", how="outer")
    return catalog_df, sample_df.sort_values("research_id").reset_index(drop=True)


def hg_na_mask(ids: pd.Series) -> pd.Series:
    s = ids.astype(str).str.strip()
    return s.str.startswith("HG") | s.str.startswith("NA")


def write_outputs(
    catalog: pd.DataFrame,
    sample_df: pd.DataFrame,
    out_dir: Path,
    *,
    covariates: Path | None = None,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    site_tsv = out_dir / "snv_indel_site_counts.tsv"
    sample_tsv = out_dir / "snv_indel_sample_qc.tsv"
    summary_tsv = out_dir / "snv_indel_sample_qc_summary.tsv"
    site_md = out_dir / "snv_indel_site_counts.md"
    summary_md = out_dir / "snv_indel_sample_qc_summary.md"
    catalog.to_csv(site_tsv, sep="\t", index=False)
    sample_df.to_csv(sample_tsv, sep="\t", index=False)
    summary = sample_qc_strata(sample_df, covariates)
    summary.to_csv(summary_tsv, sep="\t", index=False)
    write_markdown(catalog, sample_df, site_md, summary_md, strata=summary)
    print(f"wrote {site_tsv}", flush=True)
    print(f"wrote {sample_tsv}", flush=True)
    print(f"wrote {summary_tsv}", flush=True)
    print(site_md.read_text())
    print(summary_md.read_text())
    sentence = manuscript_per_participant_sentence(summary)
    if sentence:
        print(sentence, flush=True)
    return {
        "site_tsv": site_tsv,
        "sample_tsv": sample_tsv,
        "summary_tsv": summary_tsv,
        "site_md": site_md,
        "summary_md": summary_md,
    }


def sample_qc_strata(sample_df: pd.DataFrame, covariates: Path | None = None) -> pd.DataFrame:
    def row(name: str, df: pd.DataFrame) -> dict[str, Any]:
        empty = len(df) == 0

        def mean_of(col: str) -> float:
            return float("nan") if empty else float(df[col].mean())

        def sd_of(col: str) -> float:
            return float("nan") if empty else float(df[col].std())

        return {
            "stratum": name,
            "n": int(len(df)),
            "mean_n_snp_auto": mean_of("auto_n_snp"),
            "sd_n_snp_auto": sd_of("auto_n_snp"),
            "mean_n_indel_auto": mean_of("auto_n_indel"),
            "sd_n_indel_auto": sd_of("auto_n_indel"),
            "mean_n_snp_nuclear": mean_of("nuclear_n_snp"),
            "sd_n_snp_nuclear": sd_of("nuclear_n_snp"),
            "mean_n_indel_nuclear": mean_of("nuclear_n_indel"),
            "sd_n_indel_nuclear": sd_of("nuclear_n_indel"),
            "mean_ti_tv_auto": mean_of("auto_ti_tv"),
            "mean_het_hom_snp_auto": mean_of("auto_het_hom_snp"),
        }

    control = hg_na_mask(sample_df["research_id"])
    rows = [
        row("all VCF samples", sample_df),
        row("non-control (not HG*/NA*)", sample_df.loc[~control]),
    ]
    if covariates is None or not Path(covariates).is_file():
        return pd.DataFrame(rows)
    cov = pd.read_csv(covariates, dtype={"research_id": str}, low_memory=False)
    if "research_id" not in cov.columns:
        return pd.DataFrame(rows)
    merged = sample_df.merge(cov, on="research_id", how="inner")
    rows.append(row("in covariates", merged))
    if "final_releasable_v9" in merged.columns and "technology" in merged.columns:
        rel = merged["final_releasable_v9"]
        if rel.dtype == bool:
            rel_mask = rel
        else:
            rel_mask = rel.astype(str).str.lower().isin({"true", "1", "yes"})
        discovery = merged.loc[
            rel_mask & merged["technology"].astype(str).str.contains("PacBio", case=False, na=False)
        ]
        rows.append(row("Phase 2 PacBio discovery", discovery))
    return pd.DataFrame(rows)


def manuscript_per_participant_sentence(strata: pd.DataFrame) -> str:
    """Nuclear SNV/indel means for the manuscript 'per participant' sentence."""
    preferred = ("Phase 2 PacBio discovery", "non-control (not HG*/NA*)", "all VCF samples")
    hit = None
    for name in preferred:
        sub = strata.loc[strata["stratum"] == name]
        if len(sub) and int(sub.iloc[0]["n"]) > 0:
            hit = sub.iloc[0]
            break
    if hit is None:
        return ""
    n_snp = hit["mean_n_snp_nuclear"]
    n_indel = hit["mean_n_indel_nuclear"]
    if pd.isna(n_snp) or pd.isna(n_indel):
        return ""
    return (
        f"Manuscript sentence ({hit['stratum']}, n={int(hit['n']):,}; nuclear means):\n"
        f"yielding approximately {n_snp:,.0f} SNVs and {n_indel:,.0f} indels per participant."
    )


def write_markdown(
    catalog: pd.DataFrame,
    sample_df: pd.DataFrame,
    site_md: Path,
    summary_md: Path,
    *,
    strata: pd.DataFrame,
) -> None:
    row = catalog.iloc[0]

    def fmt_int(value: Any) -> str:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return "NA"
        return f"{int(value):,}"

    def mean_sd_val(value: Any, *, digits: int) -> str:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return "NA"
        return f"{float(value):,.{digits}f}"

    site_md.write_text(
        "\n".join(
            [
                "# DeepVariant + GLnexus SNV/indel site counts (bcftools stats)",
                "",
                "All samples in the joint VCF. Sites with FILTER PASS or unfiltered (`.`);",
                "GLnexus typically does not set PASS. `bcftools stats -s - -f PASS,.`.",
                "One scan per chromosome shard. Indel length bins (<20 / <50 bp) come from",
                "the IDD section (site catalog only; per-sample PSC n_indel is unbinned).",
                "",
                f"- Samples: {len(sample_df):,}",
                f"- Intervals: {int(row['n_intervals']):,}",
                "",
                "| Metric | Autosomes | Nuclear (chr1–22, X, Y) |",
                "|---|---:|---:|",
                f"| SNVs | {fmt_int(row['auto_n_snv'])} | {fmt_int(row['nuclear_n_snv'])} |",
                f"| Insertions < 20 bp | {fmt_int(row['auto_n_ins_lt20'])} | {fmt_int(row['nuclear_n_ins_lt20'])} |",
                f"| Deletions < 20 bp | {fmt_int(row['auto_n_del_lt20'])} | {fmt_int(row['nuclear_n_del_lt20'])} |",
                f"| Insertions < 50 bp | {fmt_int(row['auto_n_ins_lt50'])} | {fmt_int(row['nuclear_n_ins_lt50'])} |",
                f"| Deletions < 50 bp | {fmt_int(row['auto_n_del_lt50'])} | {fmt_int(row['nuclear_n_del_lt50'])} |",
                f"| All insertions | {fmt_int(row['auto_n_ins'])} | {fmt_int(row['nuclear_n_ins'])} |",
                f"| All deletions | {fmt_int(row['auto_n_del'])} | {fmt_int(row['nuclear_n_del'])} |",
                f"| Records | {fmt_int(row['auto_n_sites'])} | {fmt_int(row['nuclear_n_sites'])} |",
                f"| Ti/Tv (SNV sites) | {row['auto_ti_tv']:.3f} | {row['nuclear_ti_tv']:.3f} |",
                "",
            ]
        )
        + "\n"
    )
    def mean_pm_sd(mean_value: Any, sd_value: Any, *, digits: int) -> str:
        if mean_value is None or (isinstance(mean_value, float) and pd.isna(mean_value)):
            return "NA"
        if sd_value is None or (isinstance(sd_value, float) and pd.isna(sd_value)):
            return f"{float(mean_value):,.{digits}f}"
        return f"{float(mean_value):,.{digits}f} ± {float(sd_value):,.{digits}f}"

    table_lines = [
        "| Group | n | SNVs (nuclear) | indels (nuclear) | SNVs (auto) | indels (auto) | Ti/Tv (auto) | het/hom (auto) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, srow in strata.iterrows():
        table_lines.append(
            "| {stratum} | {n:,} | {snp_n} | {indel_n} | {snp_a} | {indel_a} | {ti} | {het} |".format(
                stratum=srow["stratum"],
                n=int(srow["n"]),
                snp_n=mean_pm_sd(srow.get("mean_n_snp_nuclear"), srow.get("sd_n_snp_nuclear"), digits=0),
                indel_n=mean_pm_sd(srow.get("mean_n_indel_nuclear"), srow.get("sd_n_indel_nuclear"), digits=0),
                snp_a=mean_pm_sd(srow.get("mean_n_snp_auto"), srow.get("sd_n_snp_auto"), digits=0),
                indel_a=mean_pm_sd(srow.get("mean_n_indel_auto"), srow.get("sd_n_indel_auto"), digits=0),
                ti=mean_sd_val(srow.get("mean_ti_tv_auto"), digits=3),
                het=mean_sd_val(srow.get("mean_het_hom_snp_auto"), digits=3),
            )
        )
    sentence = manuscript_per_participant_sentence(strata)
    summary_md.write_text(
        "\n".join(
            [
                "# Per-sample DeepVariant + GLnexus QC (bcftools stats)",
                "",
                "PSC columns: het / hom-alt SNVs, transitions, transversions, indel genotypes.",
                "Non-control drops VCF sample IDs starting with `HG` or `NA`.",
                "Nuclear means are chr1–22, X, Y (catalog-matched). Autosomal columns exclude X/Y.",
                "Per-sample `n_indel` is unbinned.",
                "",
                *table_lines,
                "",
                sentence,
                "",
            ]
        )
        + "\n"
    )


def load_manifest(args: argparse.Namespace) -> pd.DataFrame:
    from resolve_gl_interval_manifest import (
        DEFAULT_ENTITY_TYPE,
        DEFAULT_NAMESPACE,
        DEFAULT_WORKSPACE,
        load_gl_interval_manifest_tsv,
        normalize_gl_interval_manifest,
        resolve_gl_interval_manifest,
    )

    if args.manifest_tsv:
        manifest = load_gl_interval_manifest_tsv(args.manifest_tsv)
        source = str(args.manifest_tsv)
        out = normalize_gl_interval_manifest(manifest, autosomes_only=False, source_label=source)
    elif args.from_firecloud:
        out = resolve_gl_interval_manifest(
            namespace=args.namespace or DEFAULT_NAMESPACE,
            workspace=args.workspace or DEFAULT_WORKSPACE,
            entity_type=args.entity_type or DEFAULT_ENTITY_TYPE,
            from_firecloud=True,
            autosomes_only=False,
        )
    else:
        raise ValueError("Provide --manifest-tsv, --from-firecloud, or --merge-dir only")
    keep = AUTOSOMES if args.autosomes_only else NUCLEAR
    out = out.loc[~out["interval_id"].isin(DROP)].copy()
    extra = sorted(set(out["interval_id"]) - keep)
    if extra:
        out = out.loc[out["interval_id"].isin(keep)].copy()
    return out.reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-tsv", type=Path)
    parser.add_argument("--from-firecloud", action="store_true")
    parser.add_argument("--namespace")
    parser.add_argument("--workspace")
    parser.add_argument("--entity-type")
    parser.add_argument("--merge-dir", type=Path, help="Existing shard stats directory; skip runs")
    parser.add_argument("--shard-dir", type=Path, default=Path("snv_stats/bcftools_shards"))
    parser.add_argument("--out-dir", type=Path, default=Path("summaries/manuscript"))
    parser.add_argument("--jobs", type=int, default=int(os.environ.get("SNV_JOBS", "8")))
    parser.add_argument("--threads", type=int, default=int(os.environ.get("SNV_BCFTOOLS_THREADS", "2")))
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--pull-stats",
        action="store_true",
        help="Download the GL_INTERVAL_set stats column into the shard directory.",
    )
    parser.add_argument("--stats-column", default=STATS_COLUMN)
    parser.add_argument("--covariates", type=Path)
    parser.add_argument(
        "--pass-only",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="If set, pass -f PASS,. (GLnexus FILTER is usually '.', not PASS).",
    )
    parser.add_argument("--autosomes-only", action="store_true")
    args = parser.parse_args()

    shard_dir = args.merge_dir or args.shard_dir
    shard_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.merge_dir and not args.manifest_tsv and not args.from_firecloud:
        intervals = sorted(
            p.stem.replace(".stats", "") if p.name.endswith(".stats.txt") else p.stem
            for p in shard_dir.glob("*.stats.txt")
        )
        intervals = [i[: -len(".stats")] if i.endswith(".stats") else i for i in intervals]
        chrom_manifest = pd.DataFrame({"interval_id": intervals})
    else:
        chrom_manifest = load_manifest(args)
        intervals = chrom_manifest["interval_id"].tolist()

    print(f"shards: {len(intervals)}  jobs={args.jobs}  threads/job={args.threads}", flush=True)
    if args.pull_stats:
        pull = pull_stats_from_table(
            chrom_manifest,
            shard_dir,
            stats_column=args.stats_column,
            jobs=args.jobs,
            force=args.force,
        )
        print(f"pulled stats: {pull}", flush=True)
    if args.run and "VCF" in chrom_manifest.columns:
        jobs = [
            {
                "interval_id": row["interval_id"],
                "vcf_uri": row["VCF"],
                "out_txt": str(shard_dir / f"{row['interval_id']}.stats.txt"),
                "threads": args.threads,
                "pass_only": args.pass_only,
                "force": args.force,
            }
            for _, row in chrom_manifest.iterrows()
        ]
        if args.jobs <= 1:
            results = [_run_one_shard_star(job) for job in jobs]
        else:
            results = []
            with ProcessPoolExecutor(max_workers=args.jobs) as pool:
                futures = {pool.submit(_run_one_shard_star, job): job["interval_id"] for job in jobs}
                for fut in as_completed(futures):
                    interval_id = futures[fut]
                    try:
                        result = fut.result()
                    except Exception as exc:  # noqa: BLE001
                        raise RuntimeError(f"{interval_id} failed: {exc}") from exc
                    print(f"  [{result['status']}] {result['interval_id']}", flush=True)
                    results.append(result)
        skipped = sum(1 for r in results if r["status"] == "skip")
        print(f"finished {len(results)} shards ({skipped} resumed)", flush=True)
    elif args.run:
        print("No VCF column; merge only.", flush=True)

    intervals = select_merge_intervals(intervals, autosomes_only=args.autosomes_only)
    catalog, sample_df = merge_shards(shard_dir, intervals)
    write_outputs(catalog, sample_df, args.out_dir, covariates=args.covariates)


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
