#!/usr/bin/env python3
"""Reproduce the final downstream atlas exploration.

This is intentionally the only executable script in the atlas package. It starts
from the 1,587-GWAS G-Atlas score-atlas parquet and package-local public annotation cache,
then regenerates all outputs.
"""

from __future__ import annotations

import concurrent.futures
import bisect
import builtins
from collections import Counter, defaultdict
import csv
from datetime import datetime
import gzip
import json
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import tarfile
import time
import urllib.request
import warnings
import zipfile

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import requests

warnings.filterwarnings("ignore", category=FutureWarning)

try:  # Optional; public MPRA element overlap still runs without liftover.
    from pyliftover import LiftOver
except ImportError:  # pragma: no cover
    LiftOver = None

try:  # Optional; only needed if the MPRAVarDB cache is absent.
    import websocket
except ImportError:  # pragma: no cover
    websocket = None


def print(*args: object, sep: str = " ", end: str = "\n", file: object | None = None, flush: bool = True) -> None:
    """Timestamp every log line written by this reproduction script."""
    target = file if file is not None else None
    text = sep.join(str(arg) for arg in args)
    lines = text.splitlines() or [""]
    for i, line in enumerate(lines):
        line_end = "\n" if i < len(lines) - 1 else end
        builtins.print(f"[{datetime.now():%H:%M:%S}] {line}", end=line_end, file=target, flush=flush)


def section(step: int, total: int, title: str, kind: str = "MAIN") -> None:
    print("=" * 70)
    print(f"ATLAS EXPLORATION {kind} STEP {step}/{total}: {title}")
    print("=" * 70)


def display_path(path: Path) -> str:
    return os.path.relpath(path, ROOT)


ROOT = Path(__file__).resolve().parent

RAW_ATLAS_INPUT = ROOT / "inputs" / "no_window_gatlas_1587_gwas.parquet"
OUT = ROOT / "outputs"
SUPP = OUT / "supplementary"
MECH = OUT / "mechanism_readiness"
INTERMEDIATE = OUT / "intermediate"
INPUT = INTERMEDIATE / "atlas_evidence_disjoint_screen_variants.tsv"

MAIN_RESULTS_OUT = OUT / "atlas_main_results.tsv"
TIERED_STACK_OUT = OUT / "atlas_tiered_mechanism_stack.tsv"
VARIANT_EVIDENCE_OUT = OUT / "atlas_variant_evidence.tsv"
PRIORITIZED_CANDIDATES_OUT = OUT / "atlas_prioritized_locus_gene_candidates.tsv"
REPORT_MD_OUT = OUT / "atlas_report.md"
MECH_REPORT_OUT = MECH / "atlas_mechanism_readiness_report.md"
MECH_DISCOVERY_TAXONOMY_OUT = MECH / "mechanism_discovery_taxonomy.tsv"
MECH_TARGET_GENE_RECURRENCE_OUT = MECH / "target_gene_recurrence.tsv"
MECH_LOCUS_AMBIGUITY_OUT = MECH / "locus_ambiguity.tsv"
MECH_READY_LOCI_OUT = MECH / "mechanism_ready_loci.tsv"
MECH_EXPERIMENT_CANDIDATES_OUT = MECH / "experiment_ready_candidates.tsv"
MECH_PUBLIC_FUNCTIONAL_SUPPORT_OUT = MECH / "public_functional_support.tsv"
MECH_PUBLIC_FUNCTIONAL_SUMMARY_OUT = MECH / "public_functional_support_summary.tsv"
MECH_PRIORITIZED_CANDIDATES_OUT = MECH / "prioritized_candidate_mechanisms.tsv"
MECH_CANDIDATE_DOSSIERS_OUT = MECH / "candidate_mechanism_dossiers.md"

SUMMARY_OUT = SUPP / "locus_convergence_summary.tsv"
NULL_DRAWS_OUT = SUPP / "locus_convergence_null_draws.tsv"
DENOMINATORS_OUT = SUPP / "locus_convergence_denominators.tsv"
LOCUS_REPORT_OUT = SUPP / "locus_convergence_report.txt"
OT_SUMMARY_OUT = SUPP / "opentargets26_l2g_e2g_summary.tsv"
OT_LABEL_NULL_OUT = SUPP / "opentargets26_credible_set_coordinate_null.tsv"
OT_GENE_NULL_OUT = SUPP / "opentargets26_l2g_e2g_gene_null.tsv"
OT_VARIANT_METRICS_OUT = SUPP / "opentargets26_l2g_e2g_variant_metrics.tsv"
OT_CREDIBLE_HITS_OUT = SUPP / "opentargets26_credible_set_hits.tsv"
OT_LD_PROXY_HITS_OUT = SUPP / "opentargets26_ld_proxy_hits.tsv"
OT_NOVELTY_SUMMARY_OUT = SUPP / "opentargets26_ld_proxy_novelty_summary.tsv"
OT_REPORT_OUT = SUPP / "opentargets26_l2g_e2g_report.txt"
OT_COLOC_SUMMARY_OUT = SUPP / "opentargets26_colocalisation_summary.tsv"
OT_COLOC_GENE_NULL_OUT = SUPP / "opentargets26_colocalisation_gene_null.tsv"
OT_COLOC_HITS_OUT = SUPP / "opentargets26_colocalisation_hits.tsv"
OT_COLOC_REPORT_OUT = SUPP / "opentargets26_colocalisation_report.txt"
OT_E2G_INTERVAL_HITS_OUT = SUPP / "opentargets26_enhancer_to_gene_interval_hits.tsv"
OT_EXTENDED_SUMMARY_OUT = SUPP / "opentargets26_extended_evidence_summary.tsv"
OT_EXTENDED_VARIANT_METRICS_OUT = SUPP / "opentargets26_extended_variant_metrics.tsv"
EQTL_CATALOGUE_HITS_OUT = SUPP / "eqtl_catalogue_r8_beta_candidate_hits.tsv"
EQTL_CATALOGUE_SUMMARY_OUT = SUPP / "eqtl_catalogue_r8_beta_summary.tsv"
EQTL_CATALOGUE_GENE_NULL_OUT = SUPP / "eqtl_catalogue_r8_beta_gene_null.tsv"
EQTL_CATALOGUE_REPORT_OUT = SUPP / "eqtl_catalogue_r8_beta_report.txt"
ABC_NASSER_HITS_OUT = SUPP / "abc_nasser2021_candidate_hits.tsv"
ABC_NASSER_SUMMARY_OUT = SUPP / "abc_nasser2021_summary.tsv"
ABC_NASSER_GENE_NULL_OUT = SUPP / "abc_nasser2021_gene_null.tsv"
ABC_NASSER_REPORT_OUT = SUPP / "abc_nasser2021_report.txt"
LEAVE_ONE_OUT_PREDICTIONS_OUT = SUPP / "leave_one_resource_out_target_gene_predictions.tsv"
LEAVE_ONE_OUT_SUMMARY_OUT = SUPP / "leave_one_resource_out_target_gene_summary.tsv"
LEAVE_ONE_OUT_NULL_DRAWS_OUT = SUPP / "leave_one_resource_out_target_gene_null_draws.tsv"
LEAVE_ONE_OUT_REPORT_OUT = SUPP / "leave_one_resource_out_target_gene_report.txt"
CONDITIONAL_LOCUS_TABLE_OUT = SUPP / "conditional_recurrence_locus_table.tsv"
CONDITIONAL_LAYER_YIELDS_OUT = SUPP / "conditional_recurrence_layer_yields.tsv"
CONDITIONAL_MODEL_RESULTS_OUT = SUPP / "conditional_recurrence_model_results.tsv"
CONDITIONAL_REPORT_OUT = SUPP / "conditional_recurrence_report.txt"
PUBLIC_FUNCTIONAL_MPRAVARDB_HITS_OUT = SUPP / "public_functional_mpravardb_hits.tsv"
PUBLIC_FUNCTIONAL_MPRABASE_HITS_OUT = SUPP / "public_functional_mprabase_element_hits.tsv.gz"
PUBLIC_FUNCTIONAL_CRISPRI_HITS_OUT = SUPP / "public_functional_crispri_flowfish_hits.tsv"
PUBLIC_FUNCTIONAL_OT_CRISPR_HITS_OUT = SUPP / "public_functional_opentargets_crispr_gene_context_hits.tsv"
RAW_ATLAS_QC_OUT = INTERMEDIATE / "raw_atlas_qc.tsv"

N_PERMUTATIONS = 2000
RANDOM_SEED = 20260425
LOCUS_SIZE = 500_000
OPENTARGETS_RELEASE = "26.03"
OPENTARGETS_BASE_URL = f"https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{OPENTARGETS_RELEASE}/output"
CACHE_ROOT = (ROOT / "cache_v3").resolve()
PUBLIC_CACHE_BOOTSTRAP = True
CACHE_RELEASE_STAMP = "20260426T092123Z"
PUBLIC_CACHE_RELEASE_BASE_URL = f"https://gatlas.szags.uk/g-atlas/atlas-exploration-final/{CACHE_RELEASE_STAMP}/resources"
CACHE_RELEASE_FILE_MANIFEST = CACHE_ROOT / "cache_release_manifest.tsv"
CACHE_RELEASE_ARCHIVE_DIR = CACHE_ROOT / ".cache_release_archives"
CACHE_RESOURCE_NAMES = [
    "abc_nasser2021",
    "causaldb",
    "eqtl_catalogue_r8_beta_susie",
    "gtex_v8",
    "gwas_catalog",
    "opentargets_platform_26_03",
    "public_functional_support",
    "screen_ccre",
    "trait_family_map",
]
OPENTARGETS_CACHE = CACHE_ROOT / f"opentargets_platform_{OPENTARGETS_RELEASE.replace('.', '_')}"
PUBLIC_FUNCTIONAL_CACHE = CACHE_ROOT / "public_functional_support"
CRISPRI_FLOWFISH_FILE = PUBLIC_FUNCTIONAL_CACHE / "crispri_flowfish_all_cell_types.tsv"
MPRABASE_DB = PUBLIC_FUNCTIONAL_CACHE / "mprabase_v4_9.3.db"
MPRAVARDB_CSV = PUBLIC_FUNCTIONAL_CACHE / "mpravardb_all_mpra_data.csv"
HG38_TO_HG19_CHAIN = PUBLIC_FUNCTIONAL_CACHE / "hg38ToHg19.over.chain.gz"
HG19_TO_HG38_CHAIN = PUBLIC_FUNCTIONAL_CACHE / "hg19ToHg38.over.chain.gz"
OT_CRISPR_CACHE = OPENTARGETS_CACHE / "evidence_crispr"
OT_CRISPR_SCREEN_CACHE = OPENTARGETS_CACHE / "evidence_crispr_screen"
CRISPRI_FLOWFISH_URL = (
    "https://raw.githubusercontent.com/EngreitzLab/ABC-GWAS-Paper/main/"
    "comparePredictorsToCRISPRData/comparisonRuns/AllCellTypes-ABC_comparison/"
    "experimentalData/experimentalData.AllCellTypes.txt"
)
MPRABASE_URL = "https://zenodo.org/api/records/10920747/files/mprabase_v4_9.3.db/content"
MPRAVARDB_URL = "https://mpravardb.rc.ufl.edu/"
HG38_TO_HG19_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
HG19_TO_HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
EQTL_CATALOGUE_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/r8_beta/susie/"
EQTL_CATALOGUE_METADATA_URL = (
    "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/"
    "data_tables/dataset_metadata_r8_beta.tsv"
)
EQTL_CATALOGUE_CACHE = CACHE_ROOT / "eqtl_catalogue_r8_beta_susie"
ABC_NASSER_URL = (
    "https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/"
    "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
)
ABC_NASSER_CACHE = CACHE_ROOT / "abc_nasser2021"
ABC_NASSER_FILE = ABC_NASSER_CACHE / "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
GTEX_CACHE = CACHE_ROOT / "gtex_v8"
GTEX_EQTL_TAR = GTEX_CACHE / "GTEx_Analysis_v8_eQTL.tar"
GTEX_SQTL_TAR = GTEX_CACHE / "GTEx_Analysis_v8_sQTL.tar"
SCREEN_CACHE = CACHE_ROOT / "screen_ccre"
SCREEN_GENE_LINK_ZIP = SCREEN_CACHE / "Human-Gene-Links.zip"
GWAS_CATALOG_ZIP = CACHE_ROOT / "gwas_catalog" / "gwas-catalog-associations-full.zip"
TRAIT_FAMILY_MAP = CACHE_ROOT / "trait_family_map" / "gcst_trait_families.json"
CAUSALDB_HM3 = CACHE_ROOT / "causaldb" / "hm3_credible_sets.parquet"
CAUSALDB_OVERLAP = CACHE_ROOT / "causaldb" / "overlap_credible_sets.parquet"

BLOOD_FAMILIES = {"hematological", "blood_protein", "blood_pressure"}
NEURAL_FAMILIES = {"brain_imaging", "psychiatric", "neurological"}
SYSTEMIC_FAMILIES = {
    "metabolic",
    "cardiovascular",
    "renal",
    "lipids",
    "immune_markers",
    "autoimmune",
    "infection",
    "liver_gi",
    "anthropometric",
    "longevity",
    "reproductive",
    "allergy_respiratory",
    "cancer",
    "metabolomics",
}
NON_BLOOD_AXES = {
    "systemic_non_neural",
    "neuro_systemic",
    "balanced_cross_domain",
    "neural_enriched",
    "systemic_other",
}
SCREEN_CCRE_FILES = {
    "PLS": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.PLS.bed",
    "pELS": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.pELS.bed",
    "dELS": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.dELS.bed",
    "CA-CTCF": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.CA-CTCF.bed",
    "CA-H3K4me3": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.CA-H3K4me3.bed",
    "CA-TF": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.CA-TF.bed",
    "CA": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.CA.bed",
    "TF": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.TF.bed",
    "CTCF-bound": "https://downloads.wenglab.org/GRCh38-cCREs.CTCF-bound.bed",
}
SCREEN_LINK_FILES = {
    "3D_chromatin": "V4-hg38.Gene-Links.3D-Chromatin.txt",
    "CRISPR": "V4-hg38.Gene-Links.CRISPR.txt",
    "SCREEN_eQTL": "V4-hg38.Gene-Links.eQTLs.txt",
}
LOW_INFO_GENE_PREFIXES = (
    "ENSG",
    "RP",
    "LINC",
    "LOC",
    "AC",
    "AL",
    "AP",
    "BX",
    "OR",
    "DEFB",
    "IG",
    "RNU",
    "MIR",
    "CR",
    "VN",
    "FAM90",
    "LA16",
    "WI2",
    "CT",
)
PSEUDOGENE_RE = re.compile(r".*P\d+$")


def reset_outputs() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    INTERMEDIATE.mkdir(parents=True, exist_ok=True)
    SUPP.mkdir(parents=True, exist_ok=True)
    MECH.mkdir(parents=True, exist_ok=True)


def require_file(path: Path, label: str) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing required {label}: {display_path(path)}")


def atlas_score_band(max_score: float) -> str:
    if max_score < 0.2:
        return "hidden_lt0.2"
    if max_score < 0.5:
        return "moderate_0.2_0.5"
    return "high_ge0.5"


def recurrence_bin(n_traits: int) -> str:
    if n_traits <= 4:
        return "r03_04"
    if n_traits <= 9:
        return "r05_09"
    if n_traits <= 19:
        return "r10_19"
    if n_traits <= 49:
        return "r20_49"
    return "r50_plus"


def score_bin(max_score: float) -> str:
    if max_score < 0.05:
        return "s00_0.05"
    if max_score < 0.10:
        return "s05_0.10"
    if max_score < 0.20:
        return "s10_0.20"
    if max_score < 0.50:
        return "s20_0.50"
    return "s50_plus"


def family_entropy(counter: Counter[str]) -> float:
    total = sum(counter.values())
    if total <= 0:
        return 0.0
    return float(-sum((value / total) * math.log2(value / total) for value in counter.values()))


def primary_axis_core(fams: Counter[str]) -> str:
    total = sum(fams.values())
    if total == 0:
        return "other"
    top_family, top_n = fams.most_common(1)[0]
    top_frac = top_n / total
    blood_traits = sum(fams[k] for k in BLOOD_FAMILIES)
    neural_traits = sum(fams[k] for k in NEURAL_FAMILIES)
    systemic_traits = sum(fams[k] for k in SYSTEMIC_FAMILIES)
    nonblood_families = sum(1 for k, value in fams.items() if k not in BLOOD_FAMILIES and value > 0)
    blood_frac = blood_traits / total
    if top_family in BLOOD_FAMILIES or blood_frac >= 0.50:
        return "blood_molecular"
    if neural_traits == 0 and systemic_traits >= 5 and nonblood_families >= 5 and top_family not in BLOOD_FAMILIES and top_frac <= 0.60:
        return "systemic_non_neural"
    if neural_traits >= 2 and systemic_traits >= 3 and nonblood_families >= 5 and top_frac <= 0.60:
        return "neuro_systemic"
    if len(fams) >= 7 and top_frac <= 0.40:
        return "balanced_cross_domain"
    if neural_traits > 0:
        return "neural_enriched"
    if systemic_traits > 0:
        return "systemic_other"
    return "other"


def primary_axis_landscape(fams: Counter[str]) -> str:
    total = sum(fams.values())
    if total == 0:
        return "unmapped"
    top_family, top_n = fams.most_common(1)[0]
    top_frac = top_n / total
    blood_traits = sum(fams[k] for k in BLOOD_FAMILIES)
    neural_traits = sum(fams[k] for k in NEURAL_FAMILIES)
    systemic_traits = sum(fams[k] for k in SYSTEMIC_FAMILIES)
    nonblood_families = sum(1 for k, value in fams.items() if k not in BLOOD_FAMILIES and value > 0)
    if top_family in BLOOD_FAMILIES or (blood_traits / total) >= 0.50:
        return "blood_molecular"
    if neural_traits == 0 and systemic_traits >= 3 and nonblood_families >= 4 and top_family not in BLOOD_FAMILIES and top_frac <= 0.65:
        return "systemic_non_neural"
    if neural_traits >= 1 and systemic_traits >= 2 and nonblood_families >= 4 and top_frac <= 0.65:
        return "neuro_systemic"
    if len(fams) >= 6 and top_frac <= 0.45:
        return "balanced_cross_domain"
    if neural_traits > 0:
        return "neural_enriched"
    if systemic_traits > 0:
        return "systemic_other"
    return "other"


def pos_key(chrom: object, bp: object) -> bytes:
    return f"{str(chrom).replace('chr', '')}:{int(bp)}".encode()


def variant_key_from_gtex(variant_id: bytes) -> bytes | None:
    parts = variant_id.split(b"_", 2)
    if len(parts) < 2 or not parts[0].startswith(b"chr"):
        return None
    return parts[0][3:] + b":" + parts[1]


def tissue_from_member(member_name: str, qtl_type: str) -> str:
    base = Path(member_name).name
    if qtl_type == "eQTL":
        return base.replace(".v8.signif_variant_gene_pairs.txt.gz", "")
    return base.replace(".v8.sqtl_signifpairs.txt.gz", "")


def tissue_group(tissue: str) -> str:
    if tissue.startswith("Brain_"):
        return "brain"
    if tissue.startswith("Artery_") or tissue.startswith("Heart_"):
        return "cardiovascular"
    if tissue.startswith("Adipose_") or tissue in {"Liver", "Pancreas"}:
        return "metabolic"
    if tissue == "Kidney_Cortex":
        return "renal"
    if tissue in {"Spleen", "Whole_Blood", "Cells_EBV-transformed_lymphocytes"}:
        return "immune_blood"
    if tissue in {"Lung", "Small_Intestine_Terminal_Ileum", "Colon_Sigmoid", "Colon_Transverse", "Stomach"}:
        return "epithelial_gi_lung"
    if tissue.startswith("Skin_") or tissue == "Cells_Cultured_fibroblasts":
        return "skin_fibroblast"
    if tissue in {"Muscle_Skeletal", "Nerve_Tibial"}:
        return "muscle_nerve"
    if tissue in {"Testis", "Ovary", "Uterus", "Vagina", "Prostate", "Breast_Mammary_Tissue"}:
        return "reproductive"
    return "other"


def is_clean_gene(gene: object) -> bool:
    if not isinstance(gene, str):
        return False
    text = gene.strip()
    if not text or text == "nan":
        return False
    return not any(text.startswith(prefix) for prefix in LOW_INFO_GENE_PREFIXES)


def is_clean_qtl_gene(gene: object) -> bool:
    if not is_clean_gene(gene):
        return False
    return not bool(PSEUDOGENE_RE.match(str(gene).strip()))


def split_gene_text(value: object) -> set[str]:
    if not isinstance(value, str):
        return set()
    return {item.strip() for item in value.split(";") if item.strip() and item.strip() != "nan"}


def join_gene_set(values: set[str]) -> str:
    return ";".join(sorted(values))


def list_values(value: object) -> list[object]:
    """Return a plain Python list for scalar, list, tuple, or numpy array cells."""
    if value is None:
        return []
    if isinstance(value, float) and pd.isna(value):
        return []
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (list, tuple, set)):
        return list(value)
    return [value]


def public_content_length(url: str) -> int:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        length = response.headers.get("Content-Length")
    if not length:
        raise RuntimeError(f"missing Content-Length for {url}")
    return int(length)


def download_public_archive(url: str, dest: Path, label: str) -> None:
    expected_size = public_content_length(url)
    if dest.exists():
        actual_size = dest.stat().st_size
        if actual_size == expected_size:
            print(f"{label}: cached {display_path(dest)} ({actual_size:,} bytes)")
            return
        print(f"{label}: deleting size-mismatched archive {display_path(dest)} ({actual_size:,} != {expected_size:,} bytes)")
        dest.unlink()

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".tmp")
    tmp.unlink(missing_ok=True)
    last_error: Exception | None = None
    for attempt in range(1, 6):
        try:
            print(f"{label}: downloading {url} ({expected_size:,} bytes)")
            request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(request, timeout=600) as response, tmp.open("wb") as handle:
                shutil.copyfileobj(response, handle, length=16 * 1024 * 1024)
            actual_size = tmp.stat().st_size
            if actual_size != expected_size:
                raise RuntimeError(f"downloaded size mismatch: {actual_size:,} != {expected_size:,} bytes")
            os.replace(tmp, dest)
            return
        except Exception as exc:
            last_error = exc
            tmp.unlink(missing_ok=True)
            if attempt == 5:
                break
            print(f"{label}: retrying after failed download attempt {attempt}: {exc}")
            time.sleep(2 * attempt)
    raise RuntimeError(f"{label}: failed public HTTPS download after 5 attempts: {url}") from last_error


def load_cache_release_manifest() -> dict[str, list[dict[str, object]]]:
    if not CACHE_RELEASE_FILE_MANIFEST.exists():
        return {}
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    with CACHE_RELEASE_FILE_MANIFEST.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            relative = str(row.get("relative_path", "")).strip()
            if not relative or "/" not in relative:
                continue
            resource = relative.split("/", 1)[0]
            if resource not in CACHE_RESOURCE_NAMES:
                continue
            grouped[resource].append(
                {
                    "relative_path": relative,
                    "size_bytes": int(row.get("size_bytes") or 0),
                }
            )
    return grouped


def cache_resource_missing_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    missing: list[dict[str, object]] = []
    for row in rows:
        path = CACHE_ROOT / str(row["relative_path"])
        expected_size = int(row["size_bytes"])
        if not path.exists() or path.stat().st_size != expected_size:
            missing.append(row)
    return missing


def restore_public_cache_resource(resource: str, expected_rows: list[dict[str, object]]) -> None:
    archive_name = f"cache_{resource}.tar.zst"
    archive = CACHE_RELEASE_ARCHIVE_DIR / archive_name
    url = f"{PUBLIC_CACHE_RELEASE_BASE_URL}/{archive_name}"
    download_public_archive(url, archive, f"public cache {resource}")
    print(f"public cache {resource}: extracting {display_path(archive)}")
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    subprocess.run(["tar", "-I", "zstd", "-xf", str(archive), "-C", str(CACHE_ROOT)], check=True)
    still_missing = cache_resource_missing_rows(expected_rows)
    if still_missing:
        examples = ", ".join(str(row["relative_path"]) for row in still_missing[:5])
        raise RuntimeError(f"public cache {resource}: {len(still_missing)} files still missing or size-mismatched after extraction; examples: {examples}")
    archive.unlink(missing_ok=True)
    print(f"public cache {resource}: ready ({len(expected_rows):,} files)")


def ensure_public_cache_release() -> None:
    if not PUBLIC_CACHE_BOOTSTRAP:
        print("Public cache bootstrap: disabled")
        return
    manifest = load_cache_release_manifest()
    if not manifest:
        print("Public cache bootstrap: no local frozen cache manifest found; skipping resource restore")
        return
    print(f"Public cache bootstrap: checking {display_path(CACHE_ROOT)}")
    for resource in CACHE_RESOURCE_NAMES:
        rows = manifest.get(resource, [])
        if not rows:
            continue
        missing = cache_resource_missing_rows(rows)
        if not missing:
            print(f"public cache {resource}: ready ({len(rows):,} files)")
            continue
        print(f"public cache {resource}: {len(missing):,}/{len(rows):,} files missing or size-mismatched")
        restore_public_cache_resource(resource, rows)


def ensure_upstream_annotation_cache() -> None:
    required = [
        (RAW_ATLAS_INPUT, "1,587-GWAS G-Atlas score-atlas parquet"),
        (TRAIT_FAMILY_MAP, "GCST trait-family annotation cache"),
        (GTEX_EQTL_TAR, "GTEx v8 eQTL tar"),
        (GTEX_SQTL_TAR, "GTEx v8 sQTL tar"),
        (SCREEN_GENE_LINK_ZIP, "SCREEN Human-Gene-Links zip"),
        (GWAS_CATALOG_ZIP, "GWAS Catalog association zip"),
        (CAUSALDB_HM3, "CAUSALdb HapMap3 credible-set parquet"),
        (CAUSALDB_OVERLAP, "CAUSALdb overlap credible-set parquet"),
    ]
    for path, label in required:
        require_file(path, label)
    for ccre_class in SCREEN_CCRE_FILES:
        require_file(SCREEN_CACHE / f"{ccre_class}.bed", f"SCREEN {ccre_class} cCRE BED")


def load_raw_atlas() -> pd.DataFrame:
    raw = pd.read_parquet(RAW_ATLAS_INPUT)
    raw = raw.rename(columns={"snp": "SNP"})
    raw["CHR"] = raw["CHR"].astype(str).str.replace("chr", "", regex=False)
    raw["BP"] = raw["BP"].astype(int)
    raw["position"] = raw["CHR"].astype(str) + ":" + raw["BP"].astype(str)
    return raw


def write_raw_atlas_qc(raw: pd.DataFrame) -> None:
    traits: set[str] = set()
    for values in raw["traits"]:
        traits.update(str(value) for value in list_values(values) if str(value))
    rows = [
        {"metric": "raw_atlas_variants", "value": len(raw)},
        {"metric": "unique_gwas_traits", "value": len(traits)},
        {"metric": "recurrent_variants_n_traits_ge_3", "value": int((raw["n_traits"] >= 3).sum())},
        {"metric": "core_recurrent_variants_n_traits_ge_10", "value": int((raw["n_traits"] >= 10).sum())},
        {"metric": "unique_snps", "value": raw["SNP"].nunique()},
    ]
    pd.DataFrame(rows).to_csv(RAW_ATLAS_QC_OUT, sep="\t", index=False)
    expected = {
        "raw_atlas_variants": 1_730_224,
        "unique_gwas_traits": 1_587,
        "recurrent_variants_n_traits_ge_3": 312_794,
        "core_recurrent_variants_n_traits_ge_10": 21_027,
    }
    observed = {row["metric"]: int(row["value"]) for row in rows}
    failures = [f"{key}: expected {value:,}, observed {observed.get(key):,}" for key, value in expected.items() if observed.get(key) != value]
    if failures:
        raise RuntimeError("Raw atlas QC failed: " + "; ".join(failures))


def build_family_features(raw: pd.DataFrame, core_only: bool) -> pd.DataFrame:
    fam_map = json.load(open(TRAIT_FAMILY_MAP, encoding="utf-8"))
    gcst_family = {key: value["family"] for key, value in fam_map.items()}
    gcst_trait = {key: value["mapped_trait"] for key, value in fam_map.items()}
    work = raw[raw["n_traits"].ge(10 if core_only else 3)].copy()
    rows: list[dict[str, object]] = []
    for row in work.itertuples(index=False):
        traits = [str(value) for value in list_values(getattr(row, "traits"))]
        scores = np.array(list_values(getattr(row, "scores")), dtype=float)
        fams = Counter(gcst_family.get(trait, "unmapped") for trait in traits)
        fams.pop("unmapped", None)
        total = sum(fams.values())
        top_family, top_n = ("", 0) if not fams else fams.most_common(1)[0]
        top_frac = float(top_n / total) if total else np.nan
        max_score = float(row.max_score)
        base = {
            "SNP": row.SNP,
            "CHR": str(row.CHR).replace("chr", ""),
            "BP": int(row.BP),
            "position": row.position,
            "n_traits": int(row.n_traits),
            "max_score": max_score,
            "mean_score": float(row.mean_score),
            "score_band": atlas_score_band(max_score),
            "n_families": len(fams),
            "family_entropy": family_entropy(fams),
            "top_family": top_family,
            "top_family_frac": top_frac,
            "family_distribution": "; ".join(f"{key}:{value}" for key, value in fams.most_common(20)),
            "example_traits": "; ".join(sorted({gcst_trait.get(trait, trait) for trait in traits if trait in gcst_trait})[:16 if core_only else 12]),
        }
        if core_only:
            neural_traits = sum(fams[key] for key in NEURAL_FAMILIES)
            blood_traits = sum(fams[key] for key in BLOOD_FAMILIES)
            systemic_traits = sum(fams[key] for key in SYSTEMIC_FAMILIES)
            nonblood_families = sum(1 for key, value in fams.items() if key not in BLOOD_FAMILIES and value > 0)
            base.update(
                {
                    "blood_traits": blood_traits,
                    "neural_traits": neural_traits,
                    "systemic_traits": systemic_traits,
                    "nonblood_families": nonblood_families,
                    "primary_axis": primary_axis_core(fams),
                }
            )
        else:
            base.update(
                {
                    "median_score": float(np.median(scores)) if len(scores) else 0.0,
                    "p90_score": float(np.quantile(scores, 0.9)) if len(scores) else 0.0,
                    "n_scores_ge_0.1": int((scores >= 0.1).sum()),
                    "n_scores_ge_0.2": int((scores >= 0.2).sum()),
                    "n_scores_ge_0.5": int((scores >= 0.5).sum()),
                    "recurrence_bin": recurrence_bin(int(row.n_traits)),
                    "score_bin": score_bin(max_score),
                    "primary_axis": primary_axis_landscape(fams),
                }
            )
        rows.append(base)
    return pd.DataFrame(rows)


def assign_lower_slice(row: pd.Series) -> str:
    n_traits = int(row["n_traits"])
    max_score = float(row["max_score"])
    n_families = int(row["n_families"])
    if 3 <= n_traits <= 4 and max_score >= 0.5:
        return "r03_04_high_score"
    if 5 <= n_traits <= 9 and max_score >= 0.5:
        return "r05_09_high_score"
    if 3 <= n_traits <= 4 and 0.2 <= max_score < 0.5 and n_families >= 3:
        return "r03_04_moderate_cross_family"
    if 5 <= n_traits <= 9 and 0.2 <= max_score < 0.5 and n_families >= 4:
        return "r05_09_moderate_cross_family"
    if 20 <= n_traits <= 49 and max_score < 0.2 and n_families >= 5:
        return "r20_49_hidden_cross_family"
    if n_traits >= 50 and max_score < 0.2 and n_families >= 5:
        return "r50_plus_hidden_cross_family"
    return ""


def load_gwas_catalog_rsids() -> set[str]:
    rs_re = re.compile(r"rs\d+")
    with zipfile.ZipFile(GWAS_CATALOG_ZIP) as zf:
        member = zf.namelist()[0]
        gwas = pd.read_csv(zf.open(member), sep="\t", usecols=["SNPS"], low_memory=False)
    out: set[str] = set()
    for value in gwas["SNPS"].dropna().astype(str):
        out.update(rs_re.findall(value))
    return out


def load_causal_rsids(path: Path) -> set[str]:
    return set(pd.read_parquet(path, columns=["rsid"])["rsid"].dropna().astype(str))


def add_novelty_flags(features: pd.DataFrame, gwas: set[str], hm3: set[str], overlap: set[str]) -> pd.DataFrame:
    out = features.copy()
    out["gwas_catalog_match"] = out["SNP"].astype(str).isin(gwas)
    out["causaldb_hm3_match"] = out["SNP"].astype(str).isin(hm3)
    out["causaldb_overlap_match"] = out["SNP"].astype(str).isin(overlap)
    out["absent_all_three"] = ~(out["gwas_catalog_match"] | out["causaldb_hm3_match"] | out["causaldb_overlap_match"])
    return out


def load_gtex_gene_names() -> dict[str, str]:
    names: dict[str, str] = {}
    specs = [(GTEX_EQTL_TAR, ".v8.egenes.txt.gz", "gene_id"), (GTEX_SQTL_TAR, ".v8.sgenes.txt.gz", "gene_id")]
    for tar_path, suffix, gene_col in specs:
        with tarfile.open(tar_path) as tar:
            members = [member for member in tar.getmembers() if member.name.endswith(suffix)]
            for member in members:
                fh = tar.extractfile(member)
                if fh is None:
                    continue
                with gzip.GzipFile(fileobj=fh) as gz:
                    header = gz.readline().decode().rstrip("\n").split("\t")
                    try:
                        gene_idx = header.index(gene_col)
                        name_idx = header.index("gene_name")
                    except ValueError:
                        continue
                    for line in gz:
                        parts = line.decode(errors="replace").rstrip("\n").split("\t")
                        if len(parts) <= max(gene_idx, name_idx):
                            continue
                        gene_id = parts[gene_idx].split(".")[0]
                        gene_name = parts[name_idx]
                        if gene_id and gene_name:
                            names.setdefault(gene_id, gene_name)
    return names


def stream_gtex(target_positions: set[bytes], gene_names: dict[str, str]) -> pd.DataFrame:
    variant_rows: dict[bytes, dict[str, object]] = {}
    gene_stats: dict[tuple[bytes, str], dict[str, object]] = {}
    specs = [
        ("eQTL", GTEX_EQTL_TAR, ".v8.signif_variant_gene_pairs.txt.gz", "gene_id"),
        ("sQTL", GTEX_SQTL_TAR, ".v8.sqtl_signifpairs.txt.gz", "phenotype_id"),
    ]
    for qtl_type, tar_path, suffix, gene_col in specs:
        with tarfile.open(tar_path) as tar:
            members = [member for member in tar.getmembers() if member.name.endswith(suffix)]
            for idx, member in enumerate(members, start=1):
                tissue = tissue_from_member(member.name, qtl_type)
                group = tissue_group(tissue)
                if idx % 10 == 1 or idx == len(members):
                    print(f"{qtl_type}: scanning tissue {idx}/{len(members)}")
                fh = tar.extractfile(member)
                if fh is None:
                    continue
                with gzip.GzipFile(fileobj=fh) as gz:
                    header = gz.readline().decode().rstrip("\n").split("\t")
                    variant_idx = header.index("variant_id")
                    gene_idx = header.index(gene_col)
                    for line in gz:
                        parts = line.rstrip(b"\n").split(b"\t")
                        key = variant_key_from_gtex(parts[variant_idx])
                        if key is None or key not in target_positions:
                            continue
                        if qtl_type == "eQTL":
                            gene_id = parts[gene_idx].decode(errors="replace").split(".")[0]
                        else:
                            gene_id = parts[gene_idx].decode(errors="replace").rsplit(":", 1)[-1].split(".")[0]
                        gene_name = gene_names.get(gene_id, gene_id)
                        row = variant_rows.setdefault(
                            key,
                            {
                                "position": key.decode(),
                                "n_eqtl_pairs": 0,
                                "n_sqtl_pairs": 0,
                                "eqtl_tissues": set(),
                                "sqtl_tissues": set(),
                                "eqtl_tissue_groups": set(),
                                "sqtl_tissue_groups": set(),
                                "eqtl_genes": set(),
                                "sqtl_genes": set(),
                            },
                        )
                        if qtl_type == "eQTL":
                            row["n_eqtl_pairs"] += 1
                            row["eqtl_tissues"].add(tissue)
                            row["eqtl_tissue_groups"].add(group)
                            row["eqtl_genes"].add(gene_name)
                        else:
                            row["n_sqtl_pairs"] += 1
                            row["sqtl_tissues"].add(tissue)
                            row["sqtl_tissue_groups"].add(group)
                            row["sqtl_genes"].add(gene_name)
                        if is_clean_qtl_gene(gene_name):
                            stats = gene_stats.setdefault(
                                (key, gene_name),
                                {
                                    "position": key.decode(),
                                    "gene_name": gene_name,
                                    "qtl_types": set(),
                                    "tissues": set(),
                                    "tissue_groups": set(),
                                    "n_links": 0,
                                },
                            )
                            stats["qtl_types"].add(qtl_type)
                            stats["tissues"].add(tissue)
                            stats["tissue_groups"].add(group)
                            stats["n_links"] += 1
    core_qtl_top5: dict[str, str] = {}
    if gene_stats:
        gene_rows = []
        for stats in gene_stats.values():
            gene_rows.append(
                {
                    "position": stats["position"],
                    "gene_name": stats["gene_name"],
                    "gene_n_tissue_groups": len(stats["tissue_groups"]),
                    "gene_n_tissues": len(stats["tissues"]),
                    "both_eqtl_sqtl": {"eQTL", "sQTL"}.issubset(stats["qtl_types"]),
                    "gene_n_links": int(stats["n_links"]),
                }
            )
        gene_df = pd.DataFrame(gene_rows).sort_values(
            ["position", "gene_n_tissue_groups", "gene_n_tissues", "both_eqtl_sqtl", "gene_n_links", "gene_name"],
            ascending=[True, False, False, False, False, True],
        )
        for position, group in gene_df.groupby("position", sort=False):
            core_qtl_top5[position] = ";".join(group.head(5)["gene_name"].tolist())
    rows: list[dict[str, object]] = []
    for row in variant_rows.values():
        eqtl_groups = set(row["eqtl_tissue_groups"])
        sqtl_groups = set(row["sqtl_tissue_groups"])
        all_groups = eqtl_groups | sqtl_groups
        nonbrain_groups = all_groups - {"brain"}
        rows.append(
            {
                "position": row["position"],
                "n_eqtl_pairs": row["n_eqtl_pairs"],
                "n_sqtl_pairs": row["n_sqtl_pairs"],
                "n_eqtl_tissues": len(row["eqtl_tissues"]),
                "n_sqtl_tissues": len(row["sqtl_tissues"]),
                "n_qtl_tissue_groups": len(all_groups),
                "n_nonbrain_qtl_tissue_groups": len(nonbrain_groups),
                "has_brain_qtl": "brain" in all_groups,
                "has_nonbrain_qtl": bool(nonbrain_groups),
                "has_multiorgan_qtl": len(nonbrain_groups) >= 3,
                "all_tissue_groups": ";".join(sorted(all_groups)),
                "eqtl_genes": ";".join(sorted(row["eqtl_genes"])[:60]),
                "sqtl_genes": ";".join(sorted(row["sqtl_genes"])[:60]),
                "core_qtl_gene_candidates": core_qtl_top5.get(row["position"], ""),
            }
        )
    return pd.DataFrame(rows)


def add_gtex_to_features(features: pd.DataFrame, qtl_variants: pd.DataFrame) -> pd.DataFrame:
    qcols = [
        "position",
        "n_eqtl_pairs",
        "n_sqtl_pairs",
        "n_eqtl_tissues",
        "n_sqtl_tissues",
        "n_qtl_tissue_groups",
        "n_nonbrain_qtl_tissue_groups",
        "has_brain_qtl",
        "has_nonbrain_qtl",
        "has_multiorgan_qtl",
        "all_tissue_groups",
        "eqtl_genes",
        "sqtl_genes",
        "core_qtl_gene_candidates",
    ]
    out = features.merge(qtl_variants[qcols], on="position", how="left") if not qtl_variants.empty else features.copy()
    for col in ["n_eqtl_pairs", "n_sqtl_pairs", "n_eqtl_tissues", "n_sqtl_tissues", "n_qtl_tissue_groups", "n_nonbrain_qtl_tissue_groups"]:
        out[col] = out.get(col, 0)
        out[col] = out[col].fillna(0).astype(int)
    for col in ["has_brain_qtl", "has_nonbrain_qtl", "has_multiorgan_qtl"]:
        out[col] = out.get(col, False)
        out[col] = out[col].astype(object).where(out[col].notna(), False).astype(bool)
    for col in ["all_tissue_groups", "eqtl_genes", "sqtl_genes", "core_qtl_gene_candidates"]:
        out[col] = out.get(col, "")
        out[col] = out[col].fillna("")
    out["any_gtex_qtl"] = (out["n_eqtl_pairs"] + out["n_sqtl_pairs"]) > 0
    return out


def build_position_index(features: pd.DataFrame) -> tuple[dict[str, list[int]], dict[tuple[str, int], list[int]]]:
    by_chr: dict[str, list[int]] = {}
    row_lookup: dict[tuple[str, int], list[int]] = {}
    for idx, row in features.iterrows():
        chrom = f"chr{str(row['CHR']).replace('chr', '')}"
        bp0 = int(row["BP"]) - 1
        by_chr.setdefault(chrom, []).append(bp0)
        row_lookup.setdefault((chrom, bp0), []).append(idx)
    for chrom in by_chr:
        by_chr[chrom] = sorted(set(by_chr[chrom]))
    return by_chr, row_lookup


def annotate_screen_class(
    features: pd.DataFrame,
    bed_path: Path,
    by_chr: dict[str, list[int]],
    row_lookup: dict[tuple[str, int], list[int]],
) -> tuple[set[int], dict[int, set[str]]]:
    hits: set[int] = set()
    hit_ids: dict[int, set[str]] = defaultdict(set)
    with bed_path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if chrom not in by_chr:
                continue
            start = int(parts[1])
            end = int(parts[2])
            ccre_id = parts[4] if len(parts) > 4 else parts[3]
            positions = by_chr[chrom]
            left = bisect.bisect_left(positions, start)
            right = bisect.bisect_left(positions, end)
            for bp0 in positions[left:right]:
                for row_idx in row_lookup.get((chrom, bp0), []):
                    hits.add(row_idx)
                    hit_ids[row_idx].add(ccre_id)
    return hits, hit_ids


def add_screen_ccres(features: pd.DataFrame) -> pd.DataFrame:
    out = features.reset_index(drop=True).copy()
    by_chr, row_lookup = build_position_index(out)
    hit_classes = [set() for _ in range(len(out))]
    hit_ids_by_row: dict[int, set[str]] = defaultdict(set)
    class_cols = []
    for ccre_class in SCREEN_CCRE_FILES:
        col = f"screen_{ccre_class.replace('-', '_').replace(' ', '_')}"
        class_cols.append(col)
        hits, hit_ids = annotate_screen_class(out, SCREEN_CACHE / f"{ccre_class}.bed", by_chr, row_lookup)
        out[col] = out.index.isin(hits)
        for idx in hits:
            hit_classes[idx].add(ccre_class)
        for idx, ids in hit_ids.items():
            hit_ids_by_row[idx].update(ids)
    out["any_screen_ccre"] = out[class_cols].any(axis=1)
    out["screen_ccre_classes"] = [";".join(sorted(values)) for values in hit_classes]
    out["screen_ccre_ids"] = [";".join(sorted(hit_ids_by_row.get(idx, set()))) for idx in range(len(out))]
    return out


def variant_ccre_table(features: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in features.itertuples(index=False):
        for ccre_id in str(row.screen_ccre_ids).split(";"):
            if ccre_id and ccre_id != "nan":
                rows.append({"SNP": row.SNP, "position": row.position, "ccre_id": ccre_id, "screen_ccre_classes": row.screen_ccre_classes})
    return pd.DataFrame(rows)


def stream_screen_gene_links(target_ccres: set[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with zipfile.ZipFile(SCREEN_GENE_LINK_ZIP) as zf:
        for link_type, member in SCREEN_LINK_FILES.items():
            print(f"SCREEN: scanning {member}")
            with zf.open(member) as handle:
                for raw in handle:
                    parts = raw.decode(errors="replace").rstrip("\n").split("\t")
                    if len(parts) < 8:
                        continue
                    ccre_id = parts[0]
                    if ccre_id not in target_ccres:
                        continue
                    gene_name = parts[2].strip()
                    gene_type = parts[3].strip()
                    if gene_type != "protein_coding" or not is_clean_gene(gene_name):
                        continue
                    rows.append(
                        {
                            "ccre_id": ccre_id,
                            "screen_link_type": link_type,
                            "screen_gene_id": parts[1],
                            "screen_gene_name": gene_name,
                            "screen_gene_type": gene_type,
                        }
                    )
    if not rows:
        return pd.DataFrame(columns=["ccre_id", "screen_link_type", "screen_gene_id", "screen_gene_name", "screen_gene_type"])
    return pd.DataFrame(rows).drop_duplicates()


def add_target_concordance(features: pd.DataFrame, screen_links: pd.DataFrame) -> pd.DataFrame:
    variant_ccres = variant_ccre_table(features)
    if variant_ccres.empty or screen_links.empty:
        out = features.copy()
        out["screen_linked_genes"] = ""
        out["screen_linked_gene_count"] = 0
        out["screen_link_types"] = ""
        out["qtl_gene_candidates"] = ""
        out["screen_gtex_gene_overlap"] = ""
        out["screen_gtex_gene_overlap_count"] = 0
        out["has_target_concordance"] = False
        return out
    variant_links = variant_ccres.merge(screen_links, on="ccre_id", how="inner")
    grouped_rows = []
    for position, group in variant_links.groupby("position"):
        genes = sorted(group["screen_gene_name"].dropna().astype(str).unique())
        grouped_rows.append(
            {
                "position": position,
                "screen_linked_genes": ";".join(genes[:80]),
                "screen_linked_gene_count": len(genes),
                "screen_link_types": ";".join(sorted(group["screen_link_type"].dropna().astype(str).unique())),
            }
        )
    grouped = pd.DataFrame(grouped_rows)
    out = features.merge(grouped, on="position", how="left")
    out["screen_linked_genes"] = out["screen_linked_genes"].fillna("")
    out["screen_link_types"] = out["screen_link_types"].fillna("")
    out["screen_linked_gene_count"] = out["screen_linked_gene_count"].fillna(0).astype(int)
    overlaps = []
    for row in out.itertuples(index=False):
        screen_genes = split_gene_text(row.screen_linked_genes)
        qtl_genes = split_gene_text(row.eqtl_genes) | split_gene_text(row.sqtl_genes)
        overlap = sorted(screen_genes & qtl_genes)
        overlaps.append(
            {
                "qtl_gene_candidates": ";".join(sorted(qtl_genes)[:80]),
                "screen_gtex_gene_overlap": ";".join(overlap),
                "screen_gtex_gene_overlap_count": len(overlap),
            }
        )
    out = pd.concat([out, pd.DataFrame(overlaps)], axis=1)
    out["has_target_concordance"] = out["screen_gtex_gene_overlap_count"].gt(0)
    return out


def screen_gene_sets_by_position(variant_links: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for position, group in variant_links.groupby("position"):
        all_genes = set(group["screen_gene_name"])
        non_eqtl = set(group.loc[group["screen_link_type"].ne("SCREEN_eQTL"), "screen_gene_name"])
        three_d = set(group.loc[group["screen_link_type"].eq("3D_chromatin"), "screen_gene_name"])
        crispr = set(group.loc[group["screen_link_type"].eq("CRISPR"), "screen_gene_name"])
        rows.append(
            {
                "position": position,
                "screen_all_genes_disjoint_script": join_gene_set(all_genes),
                "screen_non_eqtl_genes": join_gene_set(non_eqtl),
                "screen_3d_genes": join_gene_set(three_d),
                "screen_crispr_genes": join_gene_set(crispr),
                "screen_all_gene_count_disjoint_script": len(all_genes),
                "screen_non_eqtl_gene_count": len(non_eqtl),
                "screen_3d_gene_count": len(three_d),
                "screen_crispr_gene_count": len(crispr),
            }
        )
    return pd.DataFrame(rows)


def build_disjoint_rows(source_set: str, features: pd.DataFrame, screen_sets: pd.DataFrame) -> pd.DataFrame:
    out = features.copy()
    out["source_set"] = source_set
    out = out.merge(screen_sets, on="position", how="left")
    for col in ["screen_all_genes_disjoint_script", "screen_non_eqtl_genes", "screen_3d_genes", "screen_crispr_genes"]:
        out[col] = out[col].fillna("")
    for col in ["screen_all_gene_count_disjoint_script", "screen_non_eqtl_gene_count", "screen_3d_gene_count", "screen_crispr_gene_count"]:
        out[col] = out[col].fillna(0).astype(int)
    overlaps = []
    for row in out.itertuples(index=False):
        qtl = split_gene_text(row.qtl_gene_candidates)
        all_screen = split_gene_text(row.screen_all_genes_disjoint_script)
        non_eqtl = split_gene_text(row.screen_non_eqtl_genes)
        three_d = split_gene_text(row.screen_3d_genes)
        crispr = split_gene_text(row.screen_crispr_genes)
        overlaps.append(
            {
                "qtl_gene_count_for_disjoint": len(qtl),
                "overlap_all_screen_genes": join_gene_set(qtl & all_screen),
                "overlap_all_screen_count": len(qtl & all_screen),
                "overlap_non_eqtl_screen_genes": join_gene_set(qtl & non_eqtl),
                "overlap_non_eqtl_screen_count": len(qtl & non_eqtl),
                "overlap_3d_screen_genes": join_gene_set(qtl & three_d),
                "overlap_3d_screen_count": len(qtl & three_d),
                "overlap_crispr_screen_genes": join_gene_set(qtl & crispr),
                "overlap_crispr_screen_count": len(qtl & crispr),
            }
        )
    out = pd.concat([out, pd.DataFrame(overlaps)], axis=1)
    out["eligible_strict_nonblood_hidden_moderate"] = (
        out["absent_all_three"].astype(bool)
        & out["score_band"].isin(["hidden_lt0.2", "moderate_0.2_0.5"])
        & out["primary_axis"].isin(NON_BLOOD_AXES)
        & out["qtl_gene_count_for_disjoint"].gt(0)
        & out["screen_all_gene_count_disjoint_script"].gt(0)
    )
    out["exact_all_screen_recomputed"] = out["eligible_strict_nonblood_hidden_moderate"] & out["overlap_all_screen_count"].gt(0)
    out["exact_non_eqtl_screen"] = out["eligible_strict_nonblood_hidden_moderate"] & out["overlap_non_eqtl_screen_count"].gt(0)
    out["exact_3d_screen"] = out["eligible_strict_nonblood_hidden_moderate"] & out["overlap_3d_screen_count"].gt(0)
    out["exact_crispr_screen"] = out["eligible_strict_nonblood_hidden_moderate"] & out["overlap_crispr_screen_count"].gt(0)
    return out


def build_candidate_table_from_raw_atlas() -> pd.DataFrame:
    ensure_upstream_annotation_cache()
    print("Loading 1,587-GWAS G-Atlas score-atlas parquet")
    raw = load_raw_atlas()
    write_raw_atlas_qc(raw)
    print("Raw atlas QC passed")

    print("Building trait-family and score layers from raw atlas")
    core = build_family_features(raw, core_only=True)
    landscape = build_family_features(raw, core_only=False)
    lower = landscape.copy()
    lower["slice_label"] = lower.apply(assign_lower_slice, axis=1)
    lower = lower[lower["slice_label"].ne("")].copy()
    print(f"Core candidates: {len(core):,}; lower-recurrence candidates: {len(lower):,}")

    print("Adding exact-rsID catalog/CAUSALdb flags")
    gwas = load_gwas_catalog_rsids()
    hm3 = load_causal_rsids(CAUSALDB_HM3)
    overlap = load_causal_rsids(CAUSALDB_OVERLAP)
    core = add_novelty_flags(core, gwas, hm3, overlap)
    lower = add_novelty_flags(lower, gwas, hm3, overlap)

    print("Streaming GTEx v8 significant eQTL/sQTL target genes")
    target_positions = {pos_key(row.CHR, row.BP) for row in pd.concat([core[["CHR", "BP"]], lower[["CHR", "BP"]]], ignore_index=True).itertuples(index=False)}
    gene_names = load_gtex_gene_names()
    qtl_variants = stream_gtex(target_positions, gene_names)
    core = add_gtex_to_features(core, qtl_variants)
    lower = add_gtex_to_features(lower, qtl_variants)

    print("Annotating SCREEN cCRE overlaps")
    combined = pd.concat(
        [
            core.assign(_candidate_source="core_exact_n10"),
            lower.assign(_candidate_source="lower_recurrence_exact"),
        ],
        ignore_index=True,
        sort=False,
    )
    combined = add_screen_ccres(combined)
    variant_ccres = variant_ccre_table(combined)
    screen_links = stream_screen_gene_links(set(variant_ccres["ccre_id"])) if not variant_ccres.empty else pd.DataFrame()
    variant_links = variant_ccres.merge(screen_links, on="ccre_id", how="inner") if not variant_ccres.empty and not screen_links.empty else pd.DataFrame()
    screen_sets = screen_gene_sets_by_position(variant_links) if not variant_links.empty else pd.DataFrame(columns=["position"])

    print("Joining GTEx and SCREEN target-gene evidence")
    combined = add_target_concordance(combined, screen_links)
    core_annotated = combined[combined["_candidate_source"].eq("core_exact_n10")].drop(columns=["_candidate_source"])
    lower_annotated = combined[combined["_candidate_source"].eq("lower_recurrence_exact")].drop(columns=["_candidate_source"])

    core_disjoint = build_disjoint_rows("core_exact_n10", core_annotated, screen_sets)
    lower_disjoint = build_disjoint_rows("lower_recurrence_exact", lower_annotated, screen_sets)
    merged = pd.concat([core_disjoint, lower_disjoint], ignore_index=True, sort=False)
    columns = [
        "source_set",
        "SNP",
        "position",
        "n_traits",
        "max_score",
        "score_band",
        "primary_axis",
        "n_families",
        "absent_all_three",
        "qtl_gene_candidates",
        "qtl_gene_count_for_disjoint",
        "screen_all_gene_count_disjoint_script",
        "screen_non_eqtl_genes",
        "screen_non_eqtl_gene_count",
        "screen_3d_genes",
        "screen_crispr_genes",
        "overlap_all_screen_genes",
        "overlap_all_screen_count",
        "overlap_non_eqtl_screen_genes",
        "overlap_non_eqtl_screen_count",
        "overlap_3d_screen_genes",
        "overlap_3d_screen_count",
        "overlap_crispr_screen_genes",
        "overlap_crispr_screen_count",
        "eligible_strict_nonblood_hidden_moderate",
        "exact_all_screen_recomputed",
        "exact_non_eqtl_screen",
        "exact_3d_screen",
        "exact_crispr_screen",
    ]
    merged = merged[columns].copy()
    for col in ["eligible_strict_nonblood_hidden_moderate", "exact_all_screen_recomputed", "exact_non_eqtl_screen", "exact_3d_screen", "exact_crispr_screen"]:
        merged[col] = merged[col].astype(bool)
    merged.to_csv(INPUT, sep="\t", index=False)
    print(f"Wrote generated candidate table: {display_path(INPUT)}")
    return merged


def gene_set(value: object) -> frozenset[str]:
    if not isinstance(value, str):
        return frozenset()
    return frozenset(x.strip().upper() for x in value.split(";") if x.strip() and x.strip() != "nan")


def size_bin(n: int) -> str:
    if n <= 1:
        return "1"
    if n <= 3:
        return "2-3"
    if n <= 6:
        return "4-6"
    if n <= 12:
        return "7-12"
    return "13+"


def traits_bin(n: int) -> str:
    if n <= 4:
        return "3-4"
    if n <= 9:
        return "5-9"
    if n <= 19:
        return "10-19"
    if n <= 49:
        return "20-49"
    return "50+"


def family_bin(n: int) -> str:
    if n <= 2:
        return "1-2"
    if n <= 4:
        return "3-4"
    if n <= 7:
        return "5-7"
    return "8+"


def parse_variant_id_position(variant_id: object) -> tuple[str, int] | None:
    if not isinstance(variant_id, str):
        return None
    fields = variant_id.split("_")
    if len(fields) < 2:
        return None
    try:
        return fields[0], int(fields[1])
    except ValueError:
        return None


def normalize_coordinates(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    if ("CHR" not in work.columns or "BP" not in work.columns) and "position" in work.columns:
        parts = work["position"].astype(str).str.replace("chr", "", regex=False).str.split(":", expand=True)
        work["CHR"] = parts[0]
        work["BP"] = pd.to_numeric(parts[1], errors="coerce").astype("Int64")
    work = work.dropna(subset=["CHR", "BP"])
    work["CHR"] = work["CHR"].astype(str).str.replace("chr", "", regex=False)
    work["BP"] = work["BP"].astype(int)
    work["locus_500kb"] = work["CHR"] + ":" + ((work["BP"] // LOCUS_SIZE) * LOCUS_SIZE).astype(str)
    return work


def prepare_work(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = normalize_coordinates(df)
    if "primary_axis" in df.columns:
        df["primary_axis"] = df["primary_axis"].replace({"balanced_cross_domain": "mixed_trait_axis"})
    df["qtl_gene_set"] = df["qtl_gene_candidates"].map(gene_set)
    df["screen_non_eqtl_set"] = df["screen_non_eqtl_genes"].map(gene_set)
    df["qtl_gene_count_check"] = df["qtl_gene_set"].map(len)
    df["screen_non_eqtl_gene_count_check"] = df["screen_non_eqtl_set"].map(len)
    df["has_non_eqtl_overlap"] = [bool(a & b) for a, b in zip(df["qtl_gene_set"], df["screen_non_eqtl_set"])]

    strict = df["eligible_strict_nonblood_hidden_moderate"].astype(bool)
    perm_tested = strict & df["qtl_gene_count_check"].gt(0) & df["screen_non_eqtl_gene_count_check"].gt(0)

    denominators = pd.DataFrame(
        [
            {
                "denominator": "strict_eligible_all_screen_context",
                "definition": "eligible_strict_nonblood_hidden_moderate is true; includes variants without non-eQTL SCREEN genes",
                "variants": int(strict.sum()),
                "loci_500kb": int(df.loc[strict, "locus_500kb"].nunique()),
            },
            {
                "denominator": "permutation_tested_non_eqtl",
                "definition": "strict eligible plus at least one GTEx gene and at least one non-eQTL SCREEN gene",
                "variants": int(perm_tested.sum()),
                "loci_500kb": int(df.loc[perm_tested, "locus_500kb"].nunique()),
            },
            {
                "denominator": "observed_exact_non_eqtl",
                "definition": "permutation-tested variants with GTEx/non-eQTL SCREEN gene overlap",
                "variants": int((perm_tested & df["has_non_eqtl_overlap"]).sum()),
                "loci_500kb": int(df.loc[perm_tested & df["has_non_eqtl_overlap"], "locus_500kb"].nunique()),
            },
        ]
    )

    work = df.loc[perm_tested].copy().reset_index(drop=True)
    work["qtl_bin"] = work["qtl_gene_count_check"].map(size_bin)
    work["screen_bin"] = work["screen_non_eqtl_gene_count_check"].map(size_bin)
    work["traits_bin"] = work["n_traits"].astype(int).map(traits_bin)
    work["family_bin"] = work["n_families"].astype(int).map(family_bin)
    work["match_bin"] = (
        work["source_set"].astype(str)
        + "|"
        + work["primary_axis"].astype(str)
        + "|"
        + work["score_band"].astype(str)
        + "|t"
        + work["traits_bin"].astype(str)
        + "|f"
        + work["family_bin"].astype(str)
        + "|q"
        + work["qtl_bin"].astype(str)
        + "|s"
        + work["screen_bin"].astype(str)
    )
    return work, denominators


def list_opentargets_parts(dataset: str) -> list[str]:
    url = f"{OPENTARGETS_BASE_URL}/{dataset}/"
    html = requests.get(url, timeout=30).text
    parts = re.findall(r'href="([^"]+\.parquet)"', html)
    if not parts:
        raise RuntimeError(f"No parquet parts found for Open Targets dataset {dataset}: {url}")
    return sorted(parts)


def download_opentargets_dataset(dataset: str) -> None:
    parts = list_opentargets_parts(dataset)
    out_dir = OPENTARGETS_CACHE / dataset
    out_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    for part in parts:
        out = out_dir / part
        if not out.exists() or out.stat().st_size == 0:
            tasks.append(part)
    if not tasks:
        print(f"Open Targets {dataset}: {len(parts)} cached parts")
        return

    print(f"Open Targets {dataset}: downloading {len(tasks)} missing parts")

    def fetch(part: str) -> str:
        out = out_dir / part
        tmp = out.with_suffix(out.suffix + ".tmp")
        urllib.request.urlretrieve(f"{OPENTARGETS_BASE_URL}/{dataset}/{part}", tmp)
        os.replace(tmp, out)
        return part

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        for i, _ in enumerate(executor.map(fetch, tasks), start=1):
            if i % 25 == 0 or i == len(tasks):
                print(f"  downloaded {i}/{len(tasks)} missing {dataset} parts")


def ensure_opentargets_cache() -> None:
    for dataset in ["credible_set", "l2g_prediction", "target", "study", "colocalisation"]:
        download_opentargets_dataset(dataset)


def ensure_opentargets_extended_cache() -> None:
    for dataset in [
        "enhancer_to_gene",
    ]:
        download_opentargets_dataset(dataset)


def parse_size_label(value: str) -> int:
    text = value.strip()
    if not text or text == "-":
        return 0
    multiplier = 1
    if text[-1].upper() == "K":
        multiplier = 1024
        text = text[:-1]
    elif text[-1].upper() == "M":
        multiplier = 1024**2
        text = text[:-1]
    elif text[-1].upper() == "G":
        multiplier = 1024**3
        text = text[:-1]
    return int(float(text) * multiplier)


def list_eqtl_catalogue_credible_sets() -> list[tuple[str, str, str, str, int]]:
    root_html = requests.get(EQTL_CATALOGUE_BASE_URL, timeout=30).text
    studies = re.findall(r'href="(QTS\d+/)"', root_html)
    rows: list[tuple[str, str, str, str, int]] = []
    for study in studies:
        study_url = urljoin(EQTL_CATALOGUE_BASE_URL, study)
        study_html = requests.get(study_url, timeout=30).text
        datasets = re.findall(r'href="(QTD\d+/)"', study_html)
        for dataset in datasets:
            dataset_url = urljoin(study_url, dataset)
            html = requests.get(dataset_url, timeout=30).text
            match = re.search(
                r'href="([^"]*credible_sets\.parquet)"[^>]*>[^<]*</a></td>'
                r'<td align="right">[^<]*</td><td align="right">\s*([^<]+)</td>',
                html,
            )
            if not match:
                continue
            filename, size_label = match.groups()
            rows.append((study.strip("/"), dataset.strip("/"), urljoin(dataset_url, filename), filename, parse_size_label(size_label)))
    return rows


def ensure_eqtl_catalogue_cache() -> None:
    EQTL_CATALOGUE_CACHE.mkdir(parents=True, exist_ok=True)
    manifest = EQTL_CATALOGUE_CACHE / "manifest.tsv"
    if manifest.exists():
        rows = []
        for line in manifest.read_text().splitlines()[1:]:
            fields = line.split("\t")
            if len(fields) == 5:
                study, dataset, url, filename, size = fields
                rows.append((study, dataset, url, filename, int(size)))
            elif len(fields) == 4:
                study, dataset, url, filename = fields
                rows.append((study, dataset, url, filename, 0))
            else:
                raise ValueError(f"Unexpected eQTL Catalogue manifest row: {line}")
    else:
        print("eQTL Catalogue r8 beta: discovering SuSiE credible-set files...")
        rows = list_eqtl_catalogue_credible_sets()
        manifest.write_text(
            "study_id\tdataset_id\turl\tfilename\tsize_bytes\n"
            + "\n".join("\t".join(map(str, row)) for row in rows)
            + "\n"
        )

    missing = []
    for row in rows:
        _, dataset, _, _, expected_size = row
        out = EQTL_CATALOGUE_CACHE / f"{dataset}.credible_sets.parquet"
        if not out.exists() or out.stat().st_size == 0 or (expected_size > 0 and out.stat().st_size < min(expected_size, 1024)):
            missing.append(row)
    if not missing:
        print(f"eQTL Catalogue r8 beta: {len(rows)} credible-set files cached")
        return

    print(f"eQTL Catalogue r8 beta: downloading {len(missing)} missing credible-set files")

    def fetch(row: tuple[str, str, str, str, int]) -> str:
        _, dataset, url, _, _ = row
        out = EQTL_CATALOGUE_CACHE / f"{dataset}.credible_sets.parquet"
        for attempt in range(1, 6):
            try:
                tmp = out.with_suffix(out.suffix + ".tmp")
                if tmp.exists():
                    tmp.unlink()
                urllib.request.urlretrieve(url, tmp)
                if tmp.stat().st_size <= 0:
                    raise RuntimeError("empty download")
                os.replace(tmp, out)
                return dataset
            except Exception:
                if attempt == 5:
                    raise
                time.sleep(2 * attempt)
        return dataset

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        for i, _ in enumerate(executor.map(fetch, missing), start=1):
            if i % 25 == 0 or i == len(missing):
                print(f"  downloaded {i}/{len(missing)} missing eQTL Catalogue files")


def ensure_abc_nasser_cache() -> None:
    ABC_NASSER_CACHE.mkdir(parents=True, exist_ok=True)
    if ABC_NASSER_FILE.exists() and ABC_NASSER_FILE.stat().st_size > 0:
        print("ABC Nasser 2021: predictions cached")
        return
    print("ABC Nasser 2021: downloading 131-biosample enhancer-gene prediction file")
    tmp = ABC_NASSER_FILE.with_suffix(ABC_NASSER_FILE.suffix + ".tmp")
    urllib.request.urlretrieve(ABC_NASSER_URL, tmp)
    os.replace(tmp, ABC_NASSER_FILE)


def stage_rank(stage: object) -> int:
    ranks = {
        "APPROVAL": 8,
        "PREAPPROVAL": 7,
        "PHASE_3": 6,
        "PHASE_2_3": 5,
        "PHASE_2": 4,
        "PHASE_1_2": 3,
        "PHASE_1": 2,
        "EARLY_PHASE_1": 1,
        "IND": 1,
        "PRECLINICAL": 0,
        "UNKNOWN": -1,
    }
    return ranks.get(str(stage), -1)


def int_flag(value: object) -> int:
    if value is None:
        return 0
    try:
        if pd.isna(value):
            return 0
    except TypeError:
        pass
    return int(value)


def has_position_in_interval(pos_by_chr: dict[str, list[int]], chrom: str, start: int, end: int) -> bool:
    import bisect

    positions = pos_by_chr.get(str(chrom), [])
    if not positions:
        return False
    left = bisect.bisect_left(positions, int(start))
    return left < len(positions) and positions[left] <= int(end)


def overlapping_loci(loci: set[str], chrom: str, start: int, end: int) -> list[str]:
    if pd.isna(start) or pd.isna(end):
        return []
    start_i, end_i = int(start), int(end)
    if end_i < start_i:
        start_i, end_i = end_i, start_i
    out = []
    first = (start_i // LOCUS_SIZE) * LOCUS_SIZE
    last = (end_i // LOCUS_SIZE) * LOCUS_SIZE
    for locus_start in range(first, last + LOCUS_SIZE, LOCUS_SIZE):
        key = f"{chrom}:{locus_start}"
        if key in loci:
            out.append(key)
    return out


def feature_value(features: object, name: str) -> float:
    if features is None:
        return 0.0
    for feature in features:
        if feature.get("name") == name:
            value = feature.get("value")
            return float(value) if value is not None else 0.0
    return 0.0


def set_metrics(work: pd.DataFrame, metrics: pd.DataFrame, name: str, mask: np.ndarray) -> dict[str, object]:
    sub_work = work.reset_index(drop=True).loc[mask]
    sub_metrics = metrics.reset_index(drop=True).loc[mask]
    row: dict[str, object] = {
        "set": name,
        "variants": int(len(sub_work)),
        "loci_500kb": int(sub_work["locus_500kb"].nunique()),
    }
    indicators = [
        "ot_gwas_cs_locus_overlap_any",
        "ot_gwas_cs_position_any",
        "ot_gwas_cs_pip_ge_0_01",
        "ot_gwas_cs_pip_ge_0_05",
        "ot_gwas_cs_pip_ge_0_10",
        "ot_gwas_cs_top_pip_ge_0_50",
        "ot_l2g_gene_agree_score_ge_0_05",
        "ot_l2g_top_gene_agree",
        "ot_l2g_e2g_gene_agree",
        "ot_l2g_coloc_gene_agree",
    ]
    for indicator in indicators:
        indicator_mask = sub_metrics[indicator].to_numpy(dtype=bool)
        row[f"{indicator}_variants"] = int(indicator_mask.sum())
        row[f"{indicator}_loci"] = int(sub_work.loc[indicator_mask, "locus_500kb"].nunique())
    return row


def label_permutation_null(
    work: pd.DataFrame,
    metrics: pd.DataFrame,
    sets: dict[str, np.ndarray],
    indicators: list[str],
) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_SEED + 777)
    groups = [np.array(list(idx), dtype=np.int32) for idx in work.groupby("match_bin").indices.values()]
    loci = work["locus_500kb"].astype(str).to_numpy()
    rows: list[dict[str, object]] = []

    for set_name, mask in sets.items():
        mask = np.asarray(mask, dtype=bool)
        if not mask.any():
            continue
        for indicator in indicators:
            values = metrics[indicator].to_numpy(dtype=bool)
            observed = int(pd.Series(loci[mask & values]).nunique())
            null = np.zeros(N_PERMUTATIONS, dtype=np.int32)
            for p in range(N_PERMUTATIONS):
                shuffled = np.zeros(len(work), dtype=bool)
                for idx in groups:
                    shuffled[idx] = values[idx] if len(idx) <= 1 else values[np.random.default_rng(rng.integers(0, 2**32 - 1)).permutation(idx)]
                null[p] = int(pd.Series(loci[mask & shuffled]).nunique())
            mean = float(null.mean())
            sd = float(null.std(ddof=1))
            rows.append(
                {
                    "set": set_name,
                    "indicator": indicator,
                    "observed_loci": observed,
                    "null_locus_mean": mean,
                    "null_locus_sd": sd,
                    "fold_enrichment": observed / mean if mean > 0 else math.inf,
                    "empirical_p_upper": (float((null >= observed).sum()) + 1.0) / (len(null) + 1.0),
                    "n_permutations": N_PERMUTATIONS,
                }
            )
    return pd.DataFrame(rows)


def gene_label_permutation_null(
    work: pd.DataFrame,
    l2g_gene_sets: dict[str, list[set[str]]],
    sets: dict[str, np.ndarray],
) -> pd.DataFrame:
    exact = work["has_non_eqtl_overlap"].to_numpy(dtype=bool)
    exact_indices = np.flatnonzero(exact)
    if len(exact_indices) == 0:
        return pd.DataFrame()

    rng = np.random.default_rng(RANDOM_SEED + 1777)
    candidate_sets = work["overlap_gene_set"].tolist()
    loci = work["locus_500kb"].astype(str).to_numpy()
    overlap_bins = work["overlap_gene_set"].map(len).map(size_bin)
    group_key = (
        work["source_set"].astype(str)
        + "|"
        + work["primary_axis"].astype(str)
        + "|"
        + work["score_band"].astype(str)
        + "|g"
        + overlap_bins.astype(str)
    )
    group_series = pd.Series(group_key.to_numpy(), index=np.arange(len(work)))
    groups = [idx.to_numpy(dtype=np.int32) for _, idx in group_series.loc[exact_indices].groupby(group_series.loc[exact_indices]).groups.items()]

    rows: list[dict[str, object]] = []
    for set_name, mask in sets.items():
        mask = np.asarray(mask, dtype=bool) & exact
        if not mask.any():
            continue
        for indicator, target_sets in l2g_gene_sets.items():
            observed_hits = np.array([bool(candidate_sets[i] & target_sets[i]) for i in range(len(work))], dtype=bool)
            observed = int(pd.Series(loci[mask & observed_hits]).nunique())
            null = np.zeros(N_PERMUTATIONS, dtype=np.int32)
            for p in range(N_PERMUTATIONS):
                shuffled_sets = list(candidate_sets)
                for idx in groups:
                    if len(idx) <= 1:
                        continue
                    permuted = rng.permutation(idx)
                    for left, right in zip(idx, permuted):
                        shuffled_sets[int(left)] = candidate_sets[int(right)]
                shuffled_hits = np.array([bool(shuffled_sets[i] & target_sets[i]) for i in range(len(work))], dtype=bool)
                null[p] = int(pd.Series(loci[mask & shuffled_hits]).nunique())
            mean = float(null.mean())
            sd = float(null.std(ddof=1))
            rows.append(
                {
                    "set": set_name,
                    "indicator": indicator,
                    "observed_loci": observed,
                    "null_locus_mean": mean,
                    "null_locus_sd": sd,
                    "fold_enrichment": observed / mean if mean > 0 else math.inf,
                    "empirical_p_upper": (float((null >= observed).sum()) + 1.0) / (len(null) + 1.0),
                    "n_permutations": N_PERMUTATIONS,
                }
            )
    return pd.DataFrame(rows)


def run_locus_null(work: pd.DataFrame, seed: int) -> tuple[int, int, np.ndarray, np.ndarray]:
    qtl_sets = work["qtl_gene_set"].tolist()
    screen_sets = work["screen_non_eqtl_set"].tolist()
    loci = work["locus_500kb"].astype(str).to_numpy()
    observed_mask = np.array([bool(a & b) for a, b in zip(qtl_sets, screen_sets)], dtype=bool)
    observed_variants = int(observed_mask.sum())
    observed_loci = int(pd.Series(loci[observed_mask]).nunique())

    rng = np.random.default_rng(seed)
    null_variant_counts = np.zeros(N_PERMUTATIONS, dtype=np.int32)
    null_locus_counts = np.zeros(N_PERMUTATIONS, dtype=np.int32)
    grouped_indices = [np.array(list(idx), dtype=np.int32) for idx in work.groupby("match_bin").indices.values()]

    for p in range(N_PERMUTATIONS):
        concordant = np.zeros(len(work), dtype=bool)
        for idx in grouped_indices:
            if len(idx) <= 1:
                i = int(idx[0])
                concordant[i] = bool(qtl_sets[i] & screen_sets[i])
                continue
            permuted = rng.permutation(idx)
            for left_i, right_i in zip(idx, permuted):
                concordant[int(left_i)] = bool(qtl_sets[int(left_i)] & screen_sets[int(right_i)])
        null_variant_counts[p] = int(concordant.sum())
        null_locus_counts[p] = int(pd.Series(loci[concordant]).nunique())

    return observed_variants, observed_loci, null_variant_counts, null_locus_counts


def summarize(name: str, work: pd.DataFrame, seed: int) -> tuple[dict[str, object], pd.DataFrame]:
    observed_variants, observed_loci, null_variants, null_loci = run_locus_null(work, seed)
    mean_loci = float(null_loci.mean())
    sd_loci = float(null_loci.std(ddof=1))
    p_loci = (float((null_loci >= observed_loci).sum()) + 1.0) / (len(null_loci) + 1.0)
    fold_loci = observed_loci / mean_loci if mean_loci > 0 else math.inf
    z_loci = (observed_loci - mean_loci) / sd_loci if sd_loci > 0 else math.inf

    mean_variants = float(null_variants.mean())
    sd_variants = float(null_variants.std(ddof=1))
    p_variants = (float((null_variants >= observed_variants).sum()) + 1.0) / (len(null_variants) + 1.0)
    fold_variants = observed_variants / mean_variants if mean_variants > 0 else math.inf

    row = {
        "set": name,
        "tested_variants": int(len(work)),
        "tested_loci_500kb": int(work["locus_500kb"].nunique()),
        "observed_concordant_variants": observed_variants,
        "observed_concordant_loci_500kb": observed_loci,
        "null_variant_mean": mean_variants,
        "null_variant_sd": sd_variants,
        "variant_fold_enrichment": fold_variants,
        "variant_empirical_p_upper": p_variants,
        "null_locus_mean": mean_loci,
        "null_locus_sd": sd_loci,
        "locus_fold_enrichment": fold_loci,
        "locus_z_score": z_loci,
        "locus_empirical_p_upper": p_loci,
        "n_permutations": N_PERMUTATIONS,
    }
    draws = pd.DataFrame(
        {
            "set": name,
            "permutation": np.arange(N_PERMUTATIONS, dtype=np.int32),
            "null_concordant_variants": null_variants,
            "null_concordant_loci_500kb": null_loci,
        }
    )
    return row, draws


def run_analysis() -> pd.DataFrame:
    OUT.mkdir(parents=True, exist_ok=True)
    SUPP.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT, sep="\t", low_memory=False)
    work, denominators = prepare_work(df)
    denominators.to_csv(DENOMINATORS_OUT, sep="\t", index=False)

    sets: list[tuple[str, pd.DataFrame]] = [
        ("overall_non_eqtl_tested", work),
        ("source_core_exact_n10", work[work["source_set"].eq("core_exact_n10")]),
        ("source_lower_recurrence_exact", work[work["source_set"].eq("lower_recurrence_exact")]),
        ("score_hidden_lt0.2", work[work["score_band"].eq("hidden_lt0.2")]),
        ("score_moderate_0.2_0.5", work[work["score_band"].eq("moderate_0.2_0.5")]),
        ("axis_systemic_non_neural", work[work["primary_axis"].eq("systemic_non_neural")]),
        ("axis_neuro_systemic", work[work["primary_axis"].eq("neuro_systemic")]),
        ("axis_neural_enriched", work[work["primary_axis"].eq("neural_enriched")]),
        ("axis_systemic_other", work[work["primary_axis"].eq("systemic_other")]),
    ]

    rows = []
    draw_frames = []
    for i, (name, sub) in enumerate(sets):
        if sub.empty:
            continue
        row, draws = summarize(name, sub.reset_index(drop=True), RANDOM_SEED + 1000 + i)
        rows.append(row)
        draw_frames.append(draws)

    summary = pd.DataFrame(rows)
    summary.to_csv(SUMMARY_OUT, sep="\t", index=False)
    pd.concat(draw_frames, ignore_index=True).to_csv(NULL_DRAWS_OUT, sep="\t", index=False)

    def line(name: str) -> str:
        row = summary[summary["set"].eq(name)].iloc[0]
        return (
            f"- {name}: {int(row.observed_concordant_loci_500kb):,} observed loci vs "
            f"{row.null_locus_mean:.1f} +/- {row.null_locus_sd:.1f} null; "
            f"{row.locus_fold_enrichment:.2f}x, empirical p={row.locus_empirical_p_upper:.4g}."
        )

    report = "\n".join(
        [
            "Atlas exploration locus-collapsed evidence-disjoint permutation",
            "======================================================",
            "",
            "Denominators",
            "- strict_eligible_all_screen_context: broad strict eligible table; includes variants without non-eQTL SCREEN genes.",
            "- permutation_tested_non_eqtl: strict eligible plus at least one GTEx gene and at least one non-eQTL SCREEN gene.",
            "- observed_exact_non_eqtl: permutation-tested variants with GTEx/non-eQTL SCREEN target-gene overlap.",
            "",
            "Key locus-level results",
            line("overall_non_eqtl_tested"),
            line("score_moderate_0.2_0.5"),
            line("source_lower_recurrence_exact"),
            line("axis_systemic_non_neural"),
            "",
            "Reporting scope",
            "- 500 kb loci are broad clustering units used to audit locus dispersion.",
            "- This section reports evidence-disjoint target-gene convergence after broad-locus collapse.",
            "",
            f"Outputs: {display_path(SUMMARY_OUT)}; {display_path(NULL_DRAWS_OUT)}; {display_path(DENOMINATORS_OUT)}",
            "",
        ]
    )
    LOCUS_REPORT_OUT.write_text(report)
    print(f"Report: {display_path(LOCUS_REPORT_OUT)}")
    return summary


def run_opentargets_anchor() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Anchor final atlas candidates against Open Targets GWAS credible sets and L2G/e2G."""
    start = time.time()
    ensure_opentargets_cache()

    df = pd.read_csv(INPUT, sep="\t", low_memory=False)
    df = normalize_coordinates(df)
    df["variant_key"] = df["CHR"].astype(str) + ":" + df["BP"].astype(str)
    df["qtl_gene_set"] = df["qtl_gene_candidates"].map(gene_set)
    df["screen_non_eqtl_set"] = df["screen_non_eqtl_genes"].map(gene_set)
    df["overlap_gene_set"] = [a & b for a, b in zip(df["qtl_gene_set"], df["screen_non_eqtl_set"])]

    work, _ = prepare_work(df)
    work["variant_key"] = work["CHR"].astype(str) + ":" + work["BP"].astype(str)
    work["overlap_gene_set"] = [a & b for a, b in zip(work["qtl_gene_set"], work["screen_non_eqtl_set"])]

    pos_by_chr = {str(chrom): sorted(set(group["BP"].astype(int).tolist())) for chrom, group in df.groupby("CHR")}
    work_loci = set(work["locus_500kb"].astype(str))
    variant_to_work_indices: dict[tuple[str, int], list[int]] = {}
    for i, row in work[["CHR", "BP"]].iterrows():
        variant_to_work_indices.setdefault((str(row.CHR), int(row.BP)), []).append(i)

    metrics = pd.DataFrame({"variant_key": work["variant_key"], "locus_500kb": work["locus_500kb"]})
    bool_columns = [
        "ot_gwas_cs_position_any",
        "ot_gwas_cs_pip_ge_0_01",
        "ot_gwas_cs_pip_ge_0_05",
        "ot_gwas_cs_pip_ge_0_10",
        "ot_gwas_cs_top_pip_ge_0_50",
        "ot_gwas_cs_locus_overlap_any",
        "ot_gwas_ld_proxy_r2_ge_0_6",
        "ot_gwas_ld_proxy_r2_ge_0_8",
        "ot_l2g_gene_agree_score_ge_0_05",
        "ot_l2g_top_gene_agree",
        "ot_l2g_e2g_gene_agree",
        "ot_l2g_coloc_gene_agree",
    ]
    for column in bool_columns:
        metrics[column] = False
    metrics["ot_gwas_cs_max_pip"] = 0.0
    metrics["ot_gwas_cs_hit_count"] = 0
    metrics["ot_gwas_cs_study_loci"] = ""
    metrics["ot_gwas_ld_proxy_max_r2"] = 0.0
    metrics["ot_gwas_ld_proxy_study_loci"] = ""
    metrics["ot_l2g_max_score_for_candidate_gene"] = 0.0
    metrics["ot_l2g_agree_symbols"] = ""
    metrics["ot_l2g_score_gene_symbols"] = ""
    metrics["ot_l2g_top_gene_symbols"] = ""
    metrics["ot_l2g_e2g_gene_symbols"] = ""
    metrics["ot_l2g_coloc_feature_gene_symbols"] = ""

    study_locus_to_work_indices: dict[str, set[int]] = {}
    credible_hit_rows: list[dict[str, object]] = []
    ld_proxy_hit_rows: list[dict[str, object]] = []
    locus_overlap_counts = {locus: 0 for locus in work_loci}

    print(f"Scanning Open Targets {OPENTARGETS_RELEASE} GWAS credible sets...")
    credible_files = sorted((OPENTARGETS_CACHE / "credible_set").glob("*.parquet"))
    for file_i, file_path in enumerate(credible_files, start=1):
        table = pq.read_table(
            file_path,
            columns=[
                "studyLocusId",
                "studyId",
                "chromosome",
                "position",
                "locusStart",
                "locusEnd",
                "confidence",
                "studyType",
                "finemappingMethod",
                "variantId",
                "ldSet",
                "locus",
            ],
        )
        part = table.to_pandas()
        part = part[part["studyType"].eq("gwas")]
        if part.empty:
            continue
        for row in part.itertuples(index=False):
            chrom = str(row.chromosome)
            start_pos = row.locusStart if pd.notna(row.locusStart) else row.position
            end_pos = row.locusEnd if pd.notna(row.locusEnd) else row.position
            locus_ids = overlapping_loci(work_loci, chrom, start_pos, end_pos)
            for locus_id in locus_ids:
                locus_overlap_counts[locus_id] = locus_overlap_counts.get(locus_id, 0) + 1
            if not locus_ids and not has_position_in_interval(pos_by_chr, chrom, int(start_pos), int(end_pos)):
                continue
            if row.locus is None:
                continue
            for item in row.locus:
                if not bool(item.get("is95CredibleSet")):
                    continue
                variant_pos = parse_variant_id_position(item.get("variantId"))
                if variant_pos not in variant_to_work_indices:
                    continue
                pip = float(item.get("posteriorProbability") or 0.0)
                for work_i in variant_to_work_indices[variant_pos]:
                    metrics.at[work_i, "ot_gwas_cs_position_any"] = True
                    metrics.at[work_i, "ot_gwas_cs_max_pip"] = max(float(metrics.at[work_i, "ot_gwas_cs_max_pip"]), pip)
                    metrics.at[work_i, "ot_gwas_cs_hit_count"] = int(metrics.at[work_i, "ot_gwas_cs_hit_count"]) + 1
                    current = metrics.at[work_i, "ot_gwas_cs_study_loci"]
                    metrics.at[work_i, "ot_gwas_cs_study_loci"] = (current + ";" + row.studyLocusId).strip(";") if current else row.studyLocusId
                    study_locus_to_work_indices.setdefault(row.studyLocusId, set()).add(work_i)
                    credible_hit_rows.append(
                        {
                            "variant_key": work.at[work_i, "variant_key"],
                            "SNP": work.at[work_i, "SNP"],
                            "locus_500kb": work.at[work_i, "locus_500kb"],
                            "studyLocusId": row.studyLocusId,
                            "studyId": row.studyId,
                            "pip": pip,
                            "external_variantId": item.get("variantId"),
                            "confidence": row.confidence,
                            "finemappingMethod": row.finemappingMethod,
                            "score_band": work.at[work_i, "score_band"],
                            "source_set": work.at[work_i, "source_set"],
                            "primary_axis": work.at[work_i, "primary_axis"],
                            "overlap_non_eqtl_genes": ";".join(sorted(work.at[work_i, "overlap_gene_set"])),
                        }
                    )
            if row.ldSet is not None:
                for item in row.ldSet:
                    tag_pos = parse_variant_id_position(item.get("tagVariantId"))
                    if tag_pos not in variant_to_work_indices:
                        continue
                    r2 = float(item.get("r2Overall") or 0.0)
                    for work_i in variant_to_work_indices[tag_pos]:
                        metrics.at[work_i, "ot_gwas_ld_proxy_max_r2"] = max(float(metrics.at[work_i, "ot_gwas_ld_proxy_max_r2"]), r2)
                        current = metrics.at[work_i, "ot_gwas_ld_proxy_study_loci"]
                        metrics.at[work_i, "ot_gwas_ld_proxy_study_loci"] = (current + ";" + row.studyLocusId).strip(";") if current else row.studyLocusId
                        if r2 >= 0.6:
                            ld_proxy_hit_rows.append(
                                {
                                    "variant_key": work.at[work_i, "variant_key"],
                                    "SNP": work.at[work_i, "SNP"],
                                    "locus_500kb": work.at[work_i, "locus_500kb"],
                                    "studyLocusId": row.studyLocusId,
                                    "studyId": row.studyId,
                                    "tagVariantId": item.get("tagVariantId"),
                                    "r2Overall": r2,
                                    "lead_variantId": row.variantId,
                                    "confidence": row.confidence,
                                    "finemappingMethod": row.finemappingMethod,
                                    "score_band": work.at[work_i, "score_band"],
                                    "source_set": work.at[work_i, "source_set"],
                                    "primary_axis": work.at[work_i, "primary_axis"],
                                    "has_non_eqtl_overlap": bool(work.at[work_i, "has_non_eqtl_overlap"]),
                                }
                            )
        if file_i % 50 == 0:
            print(f"  credible_set parts {file_i}/{len(credible_files)}; coordinate hit rows {len(credible_hit_rows):,}")

    for i, locus_id in enumerate(metrics["locus_500kb"]):
        metrics.at[i, "ot_gwas_cs_locus_overlap_any"] = locus_overlap_counts.get(locus_id, 0) > 0
    metrics["ot_gwas_cs_pip_ge_0_01"] = metrics["ot_gwas_cs_max_pip"].ge(0.01)
    metrics["ot_gwas_cs_pip_ge_0_05"] = metrics["ot_gwas_cs_max_pip"].ge(0.05)
    metrics["ot_gwas_cs_pip_ge_0_10"] = metrics["ot_gwas_cs_max_pip"].ge(0.10)
    metrics["ot_gwas_cs_top_pip_ge_0_50"] = metrics["ot_gwas_cs_max_pip"].ge(0.50)
    metrics["ot_gwas_ld_proxy_r2_ge_0_6"] = metrics["ot_gwas_ld_proxy_max_r2"].ge(0.6)
    metrics["ot_gwas_ld_proxy_r2_ge_0_8"] = metrics["ot_gwas_ld_proxy_max_r2"].ge(0.8)

    print(
        "Open Targets coordinate hits: "
        f"{int(metrics['ot_gwas_cs_position_any'].sum()):,} variants, "
        f"{metrics.loc[metrics['ot_gwas_cs_position_any'], 'locus_500kb'].nunique():,} loci"
    )

    target_table = ds.dataset(str(OPENTARGETS_CACHE / "target"), format="parquet").to_table(columns=["id", "approvedSymbol"]).to_pandas()
    id_to_symbol = dict(zip(target_table["id"], target_table["approvedSymbol"].astype(str).str.upper()))
    study_loci = set(study_locus_to_work_indices)
    l2g_score_sets = [set() for _ in range(len(work))]
    l2g_top_sets = [set() for _ in range(len(work))]
    l2g_e2g_sets = [set() for _ in range(len(work))]
    l2g_coloc_sets = [set() for _ in range(len(work))]

    print(f"Scanning Open Targets L2G for {len(study_loci):,} matched study loci...")
    by_study: dict[str, list[tuple[str, float, float, float]]] = {}
    for file_path in sorted((OPENTARGETS_CACHE / "l2g_prediction").glob("*.parquet")):
        table = pq.read_table(file_path, columns=["studyLocusId", "geneId", "score", "features"])
        part = table.to_pandas()
        part = part[part["studyLocusId"].isin(study_loci)]
        if part.empty:
            continue
        for row in part.itertuples(index=False):
            symbol = id_to_symbol.get(row.geneId, "")
            if not symbol:
                continue
            e2g_score = max(feature_value(row.features, "e2gMean"), feature_value(row.features, "e2gMeanNeighbourhood"))
            coloc_score = max(
                feature_value(row.features, "eQtlColocH4Maximum"),
                feature_value(row.features, "sQtlColocH4Maximum"),
                feature_value(row.features, "pQtlColocH4Maximum"),
                feature_value(row.features, "eQtlColocClppMaximum"),
                feature_value(row.features, "sQtlColocClppMaximum"),
                feature_value(row.features, "pQtlColocClppMaximum"),
            )
            by_study.setdefault(row.studyLocusId, []).append((symbol, float(row.score), e2g_score, coloc_score))

    for study_locus_id, rows in by_study.items():
        top_score = max(row[1] for row in rows)
        top_symbols = {symbol for symbol, score, _, _ in rows if score == top_score}
        score_symbols = {symbol for symbol, score, _, _ in rows if score >= 0.05}
        e2g_symbols = {symbol for symbol, _, e2g_score, _ in rows if e2g_score > 0}
        coloc_symbols = {symbol for symbol, _, _, coloc_score in rows if coloc_score > 0.01}
        for work_i in study_locus_to_work_indices.get(study_locus_id, set()):
            l2g_score_sets[work_i].update(score_symbols)
            l2g_top_sets[work_i].update(top_symbols)
            l2g_e2g_sets[work_i].update(e2g_symbols)
            l2g_coloc_sets[work_i].update(coloc_symbols)

    for i, candidate_genes in enumerate(work["overlap_gene_set"]):
        score_agree = candidate_genes & l2g_score_sets[i]
        top_agree = candidate_genes & l2g_top_sets[i]
        e2g_agree = candidate_genes & l2g_e2g_sets[i]
        coloc_agree = candidate_genes & l2g_coloc_sets[i]
        metrics.at[i, "ot_l2g_gene_agree_score_ge_0_05"] = bool(score_agree)
        metrics.at[i, "ot_l2g_top_gene_agree"] = bool(top_agree)
        metrics.at[i, "ot_l2g_e2g_gene_agree"] = bool(e2g_agree)
        metrics.at[i, "ot_l2g_coloc_gene_agree"] = bool(coloc_agree)
        metrics.at[i, "ot_l2g_agree_symbols"] = ";".join(sorted(score_agree | top_agree | e2g_agree | coloc_agree))
        metrics.at[i, "ot_l2g_score_gene_symbols"] = ";".join(sorted(score_agree))
        metrics.at[i, "ot_l2g_top_gene_symbols"] = ";".join(sorted(top_agree))
        metrics.at[i, "ot_l2g_e2g_gene_symbols"] = ";".join(sorted(e2g_agree))
        metrics.at[i, "ot_l2g_coloc_feature_gene_symbols"] = ";".join(sorted(coloc_agree))
        if score_agree:
            matching_scores = [
                score
                for study_locus_id in str(metrics.at[i, "ot_gwas_cs_study_loci"]).split(";")
                for symbol, score, _, _ in by_study.get(study_locus_id, [])
                if symbol in candidate_genes
            ]
            metrics.at[i, "ot_l2g_max_score_for_candidate_gene"] = max(matching_scores) if matching_scores else 0.0

    sets = {
        "permutation_tested_non_eqtl": np.ones(len(work), dtype=bool),
        "observed_exact_non_eqtl": work["has_non_eqtl_overlap"].to_numpy(dtype=bool),
        "non_exact_tested": (~work["has_non_eqtl_overlap"]).to_numpy(dtype=bool),
        "moderate_exact_non_eqtl": (work["has_non_eqtl_overlap"] & work["score_band"].eq("moderate_0.2_0.5")).to_numpy(dtype=bool),
        "hidden_exact_non_eqtl": (work["has_non_eqtl_overlap"] & work["score_band"].eq("hidden_lt0.2")).to_numpy(dtype=bool),
        "lower_recurrence_exact": (work["has_non_eqtl_overlap"] & work["source_set"].eq("lower_recurrence_exact")).to_numpy(dtype=bool),
        "core_exact_n10": (work["has_non_eqtl_overlap"] & work["source_set"].eq("core_exact_n10")).to_numpy(dtype=bool),
        "systemic_non_neural_exact": (work["has_non_eqtl_overlap"] & work["primary_axis"].eq("systemic_non_neural")).to_numpy(dtype=bool),
    }
    summary = pd.DataFrame([set_metrics(work, metrics, name, mask) for name, mask in sets.items()])
    summary.to_csv(OT_SUMMARY_OUT, sep="\t", index=False)
    pd.DataFrame(ld_proxy_hit_rows).drop_duplicates().to_csv(OT_LD_PROXY_HITS_OUT, sep="\t", index=False)

    novelty_rows: list[dict[str, object]] = []
    for set_name, mask in sets.items():
        mask = np.asarray(mask, dtype=bool)
        sub_work = work.loc[mask].copy()
        sub_metrics = metrics.loc[mask].copy()
        exact_coordinate = sub_metrics["ot_gwas_cs_position_any"].to_numpy(dtype=bool)
        ld_proxy_08 = sub_metrics["ot_gwas_ld_proxy_r2_ge_0_8"].to_numpy(dtype=bool)
        ld_proxy_06 = sub_metrics["ot_gwas_ld_proxy_r2_ge_0_6"].to_numpy(dtype=bool)
        broad_locus = sub_metrics["ot_gwas_cs_locus_overlap_any"].to_numpy(dtype=bool)
        novelty_rows.append(
            {
                "set": set_name,
                "variants": int(len(sub_work)),
                "loci_500kb": int(sub_work["locus_500kb"].nunique()),
                "exact_coordinate_hit_loci": int(sub_work.loc[exact_coordinate, "locus_500kb"].nunique()),
                "ld_proxy_r2_ge_0_8_loci": int(sub_work.loc[ld_proxy_08, "locus_500kb"].nunique()),
                "ld_proxy_r2_ge_0_6_loci": int(sub_work.loc[ld_proxy_06, "locus_500kb"].nunique()),
                "not_exact_or_ld_proxy_r2_ge_0_8_loci": int(sub_work.loc[~(exact_coordinate | ld_proxy_08), "locus_500kb"].nunique()),
                "not_exact_or_ld_proxy_r2_ge_0_6_loci": int(sub_work.loc[~(exact_coordinate | ld_proxy_06), "locus_500kb"].nunique()),
                "no_external_gwas_credible_set_broad_locus_overlap_loci": int(sub_work.loc[~broad_locus, "locus_500kb"].nunique()),
            }
        )
    novelty_summary = pd.DataFrame(novelty_rows)
    novelty_summary.to_csv(OT_NOVELTY_SUMMARY_OUT, sep="\t", index=False)

    label_null = label_permutation_null(
        work,
        metrics,
        sets,
        [
            "ot_gwas_cs_position_any",
            "ot_gwas_ld_proxy_r2_ge_0_8",
            "ot_gwas_ld_proxy_r2_ge_0_6",
            "ot_gwas_cs_pip_ge_0_01",
            "ot_gwas_cs_pip_ge_0_05",
            "ot_gwas_cs_pip_ge_0_10",
        ],
    )
    label_null.to_csv(OT_LABEL_NULL_OUT, sep="\t", index=False)

    gene_null = gene_label_permutation_null(
        work,
        {
            "ot_l2g_gene_agree_score_ge_0_05": l2g_score_sets,
            "ot_l2g_top_gene_agree": l2g_top_sets,
            "ot_l2g_e2g_gene_agree": l2g_e2g_sets,
            "ot_l2g_coloc_gene_agree": l2g_coloc_sets,
        },
        sets,
    )
    gene_null.to_csv(OT_GENE_NULL_OUT, sep="\t", index=False)

    variant_metrics = pd.concat([work.reset_index(drop=True), metrics.drop(columns=["variant_key", "locus_500kb"])], axis=1)
    keep_columns = [
        "SNP",
        "variant_key",
        "locus_500kb",
        "source_set",
        "score_band",
        "primary_axis",
        "n_traits",
        "n_families",
        "has_non_eqtl_overlap",
        "ot_gwas_cs_max_pip",
        "ot_gwas_cs_hit_count",
        "ot_gwas_cs_position_any",
        "ot_gwas_cs_pip_ge_0_01",
        "ot_gwas_cs_pip_ge_0_05",
        "ot_gwas_cs_pip_ge_0_10",
        "ot_gwas_ld_proxy_max_r2",
        "ot_gwas_ld_proxy_r2_ge_0_6",
        "ot_gwas_ld_proxy_r2_ge_0_8",
        "ot_l2g_gene_agree_score_ge_0_05",
        "ot_l2g_top_gene_agree",
        "ot_l2g_e2g_gene_agree",
        "ot_l2g_coloc_gene_agree",
        "ot_l2g_max_score_for_candidate_gene",
        "ot_l2g_agree_symbols",
        "ot_l2g_score_gene_symbols",
        "ot_l2g_top_gene_symbols",
        "ot_l2g_e2g_gene_symbols",
        "ot_l2g_coloc_feature_gene_symbols",
    ]
    variant_metrics[keep_columns].to_csv(OT_VARIANT_METRICS_OUT, sep="\t", index=False)
    pd.DataFrame(credible_hit_rows).to_csv(OT_CREDIBLE_HITS_OUT, sep="\t", index=False)

    def summary_line(set_name: str) -> str:
        row = summary[summary["set"].eq(set_name)].iloc[0]
        return (
            f"- {set_name}: {int(row.ot_gwas_cs_position_any_loci):,} GWAS credible-set coordinate loci; "
            f"{int(row.ot_l2g_gene_agree_score_ge_0_05_loci):,} L2G-agreeing loci; "
            f"{int(row.ot_l2g_e2g_gene_agree_loci):,} e2G-agreeing loci."
        )

    def null_line(null_df: pd.DataFrame, set_name: str, indicator: str) -> str:
        row = null_df[(null_df["set"].eq(set_name)) & (null_df["indicator"].eq(indicator))].iloc[0]
        return (
            f"- {set_name}/{indicator}: {int(row.observed_loci):,} observed loci vs "
            f"{row.null_locus_mean:.1f} null; {row.fold_enrichment:.2f}x; "
            f"p={row.empirical_p_upper:.4g}."
        )

    def novelty_line(set_name: str) -> str:
        row = novelty_summary[novelty_summary["set"].eq(set_name)].iloc[0]
        return (
            f"- {set_name}: {int(row.not_exact_or_ld_proxy_r2_ge_0_8_loci):,}/{int(row.loci_500kb):,} loci "
            f"are not exact external GWAS credible-set coordinate hits and not r2>=0.8 LD proxies; "
            f"{int(row.no_external_gwas_credible_set_broad_locus_overlap_loci):,} loci have no broad 500 kb external GWAS credible-set overlap."
        )

    report = "\n".join(
        [
            f"Open Targets {OPENTARGETS_RELEASE} external genetics anchor",
            "================================================",
            "",
            f"External source: {OPENTARGETS_BASE_URL}",
            f"External cache: {display_path(OPENTARGETS_CACHE)}",
            "",
            "Scope",
            "- GWAS credible-set coordinate overlap uses Open Targets 95% credible-set variants by chromosome/position.",
            "- Target-gene agreement tests whether GTEx/non-eQTL SCREEN candidate genes match Open Targets L2G genes at matched GWAS credible-set coordinates.",
            "- L2G/e2G agreement is reported as an orthogonal target-assignment anchor alongside separate colocalisation and functional-support layers.",
            "",
            "Descriptive anchor counts",
            summary_line("observed_exact_non_eqtl"),
            summary_line("moderate_exact_non_eqtl"),
            summary_line("lower_recurrence_exact"),
            summary_line("systemic_non_neural_exact"),
            "",
            "Matched coordinate-overlap null",
            null_line(label_null, "observed_exact_non_eqtl", "ot_gwas_cs_position_any"),
            null_line(label_null, "moderate_exact_non_eqtl", "ot_gwas_cs_position_any"),
            null_line(label_null, "systemic_non_neural_exact", "ot_gwas_cs_position_any"),
            "",
            "LD-proxy knownness audit",
            novelty_line("observed_exact_non_eqtl"),
            novelty_line("moderate_exact_non_eqtl"),
            novelty_line("systemic_non_neural_exact"),
            "",
            "Gene-label null among exact candidate genes",
            null_line(gene_null, "observed_exact_non_eqtl", "ot_l2g_gene_agree_score_ge_0_05"),
            null_line(gene_null, "moderate_exact_non_eqtl", "ot_l2g_gene_agree_score_ge_0_05"),
            null_line(gene_null, "systemic_non_neural_exact", "ot_l2g_e2g_gene_agree"),
            "",
            "Reporting interpretation",
            "- External GWAS credible-set coordinate overlap is common; target-gene agreement is the strongest matched signal.",
            "- Atlas GTEx/non-eQTL candidate genes recur at external GWAS credible-set coordinates and agree with Open Targets L2G/e2G target-gene predictions above gene-label null.",
            "- These records prioritize regulatory target-gene hypotheses for full-stack evidence integration.",
            "",
            f"Elapsed seconds: {time.time() - start:.1f}",
            "",
        ]
    )
    OT_REPORT_OUT.write_text(report)
    print(f"Report: {display_path(OT_REPORT_OUT)}")
    return summary, label_null, gene_null


def run_opentargets_colocalisation() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Add same-gene Open Targets GWAS-molQTL colocalisation support."""
    start = time.time()

    df = pd.read_csv(INPUT, sep="\t", low_memory=False)
    df = normalize_coordinates(df)
    df["variant_key"] = df["CHR"].astype(str) + ":" + df["BP"].astype(str)
    df["qtl_gene_set"] = df["qtl_gene_candidates"].map(gene_set)
    df["screen_non_eqtl_set"] = df["screen_non_eqtl_genes"].map(gene_set)
    df["overlap_gene_set"] = [a & b for a, b in zip(df["qtl_gene_set"], df["screen_non_eqtl_set"])]
    work, _ = prepare_work(df)
    work["variant_key"] = work["CHR"].astype(str) + ":" + work["BP"].astype(str)
    work["overlap_gene_set"] = [a & b for a, b in zip(work["qtl_gene_set"], work["screen_non_eqtl_set"])]

    credible_hits = pd.read_csv(OT_CREDIBLE_HITS_OUT, sep="\t")
    left_study_loci = set(credible_hits["studyLocusId"].astype(str))
    variant_to_lefts: dict[str, set[str]] = {}
    left_to_candidate_genes: dict[str, set[str]] = {}
    left_to_loci: dict[str, set[str]] = {}
    for row in credible_hits.itertuples(index=False):
        left_id = str(row.studyLocusId)
        variant_to_lefts.setdefault(str(row.variant_key), set()).add(left_id)
        genes = gene_set(row.overlap_non_eqtl_genes)
        left_to_candidate_genes.setdefault(left_id, set()).update(genes)
        left_to_loci.setdefault(left_id, set()).add(str(row.locus_500kb))

    molecular_types = {"eqtl", "sqtl", "pqtl", "sceqtl", "tuqtl"}
    right_study_loci: set[str] = set()
    coloc_rows: list[tuple[str, str, str, float, float, float, int]] = []

    print(f"Scanning Open Targets {OPENTARGETS_RELEASE} colocalisation for {len(left_study_loci):,} GWAS study loci...")
    coloc_files = sorted((OPENTARGETS_CACHE / "colocalisation").glob("*.parquet"))
    for file_i, file_path in enumerate(coloc_files, start=1):
        table = pq.read_table(
            file_path,
            columns=[
                "leftStudyLocusId",
                "rightStudyLocusId",
                "rightStudyType",
                "h4",
                "clpp",
                "betaRatioSignAverage",
                "numberColocalisingVariants",
            ],
        )
        part = table.to_pandas()
        part = part[part["leftStudyLocusId"].isin(left_study_loci) & part["rightStudyType"].isin(molecular_types)]
        if not part.empty:
            for row in part.itertuples(index=False):
                right_id = str(row.rightStudyLocusId)
                right_study_loci.add(right_id)
                coloc_rows.append(
                    (
                        str(row.leftStudyLocusId),
                        right_id,
                        str(row.rightStudyType),
                        float(row.h4 or 0.0),
                        float(row.clpp or 0.0),
                        float(row.betaRatioSignAverage) if pd.notna(row.betaRatioSignAverage) else float("nan"),
                        int(row.numberColocalisingVariants),
                    )
                )
        if file_i % 50 == 0:
            print(f"  colocalisation parts {file_i}/{len(coloc_files)}; rows {len(coloc_rows):,}")

    right_to_study: dict[str, tuple[str, str]] = {}
    for file_path in sorted((OPENTARGETS_CACHE / "credible_set").glob("*.parquet")):
        table = pq.read_table(file_path, columns=["studyLocusId", "studyId", "studyType"])
        part = table.to_pandas()
        part = part[part["studyLocusId"].isin(right_study_loci)]
        for row in part.itertuples(index=False):
            right_to_study[str(row.studyLocusId)] = (str(row.studyId), str(row.studyType))

    target_table = ds.dataset(str(OPENTARGETS_CACHE / "target"), format="parquet").to_table(columns=["id", "approvedSymbol"]).to_pandas()
    id_to_symbol = dict(zip(target_table["id"], target_table["approvedSymbol"].astype(str).str.upper()))
    study_table = ds.dataset(str(OPENTARGETS_CACHE / "study"), format="parquet").to_table(columns=["studyId", "geneId"]).to_pandas()
    study_to_symbol = {
        str(row.studyId): id_to_symbol.get(row.geneId, "")
        for row in study_table.itertuples(index=False)
        if pd.notna(row.geneId)
    }

    left_maps: dict[str, dict[str, set[str]]] = {
        "ot_coloc_any_gene_agree": {},
        "ot_coloc_h4_ge_0_5_gene_agree": {},
        "ot_coloc_h4_ge_0_8_gene_agree": {},
        "ot_coloc_clpp_ge_0_01_gene_agree": {},
    }
    hit_rows: list[dict[str, object]] = []
    for left_id, right_id, right_type, h4, clpp, beta_ratio_sign_average, n_coloc_variants in coloc_rows:
        study_id, _ = right_to_study.get(right_id, ("", ""))
        symbol = study_to_symbol.get(study_id, "")
        if not symbol:
            continue
        left_maps["ot_coloc_any_gene_agree"].setdefault(left_id, set()).add(symbol)
        if h4 >= 0.5:
            left_maps["ot_coloc_h4_ge_0_5_gene_agree"].setdefault(left_id, set()).add(symbol)
        if h4 >= 0.8:
            left_maps["ot_coloc_h4_ge_0_8_gene_agree"].setdefault(left_id, set()).add(symbol)
        if clpp >= 0.01:
            left_maps["ot_coloc_clpp_ge_0_01_gene_agree"].setdefault(left_id, set()).add(symbol)
        if symbol in left_to_candidate_genes.get(left_id, set()):
            for locus_id in left_to_loci.get(left_id, set()):
                hit_rows.append(
                    {
                        "leftStudyLocusId": left_id,
                        "rightStudyLocusId": right_id,
                        "rightStudyType": right_type,
                        "gene_symbol": symbol,
                        "locus_500kb": locus_id,
                        "h4": h4,
                        "clpp": clpp,
                        "beta_ratio_sign_average": beta_ratio_sign_average,
                        "number_colocalising_variants": n_coloc_variants,
                        "h4_ge_0_5": h4 >= 0.5,
                        "h4_ge_0_8": h4 >= 0.8,
                        "clpp_ge_0_01": clpp >= 0.01,
                    }
                )

    per_work_sets: dict[str, list[set[str]]] = {name: [] for name in left_maps}
    for row in work.itertuples(index=False):
        lefts = variant_to_lefts.get(str(row.variant_key), set())
        for name, left_map in left_maps.items():
            genes: set[str] = set()
            for left_id in lefts:
                genes.update(left_map.get(left_id, set()))
            per_work_sets[name].append(genes)

    metrics = pd.DataFrame({"variant_key": work["variant_key"], "locus_500kb": work["locus_500kb"]})
    for name, gene_sets in per_work_sets.items():
        metrics[name] = [bool(candidate & support) for candidate, support in zip(work["overlap_gene_set"], gene_sets)]

    sets = {
        "observed_exact_non_eqtl": work["has_non_eqtl_overlap"].to_numpy(dtype=bool),
        "moderate_exact_non_eqtl": (work["has_non_eqtl_overlap"] & work["score_band"].eq("moderate_0.2_0.5")).to_numpy(dtype=bool),
        "hidden_exact_non_eqtl": (work["has_non_eqtl_overlap"] & work["score_band"].eq("hidden_lt0.2")).to_numpy(dtype=bool),
        "lower_recurrence_exact": (work["has_non_eqtl_overlap"] & work["source_set"].eq("lower_recurrence_exact")).to_numpy(dtype=bool),
        "core_exact_n10": (work["has_non_eqtl_overlap"] & work["source_set"].eq("core_exact_n10")).to_numpy(dtype=bool),
        "systemic_non_neural_exact": (work["has_non_eqtl_overlap"] & work["primary_axis"].eq("systemic_non_neural")).to_numpy(dtype=bool),
    }

    summary_rows: list[dict[str, object]] = []
    for set_name, mask in sets.items():
        sub_work = work.loc[mask].copy()
        sub_metrics = metrics.loc[mask].copy()
        row: dict[str, object] = {
            "set": set_name,
            "variants": int(len(sub_work)),
            "loci_500kb": int(sub_work["locus_500kb"].nunique()),
        }
        for indicator in per_work_sets:
            indicator_mask = sub_metrics[indicator].to_numpy(dtype=bool)
            row[f"{indicator}_variants"] = int(indicator_mask.sum())
            row[f"{indicator}_loci"] = int(sub_work.loc[indicator_mask, "locus_500kb"].nunique())
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OT_COLOC_SUMMARY_OUT, sep="\t", index=False)

    gene_null = gene_label_permutation_null(work, per_work_sets, sets)
    gene_null.to_csv(OT_COLOC_GENE_NULL_OUT, sep="\t", index=False)
    pd.DataFrame(hit_rows).to_csv(OT_COLOC_HITS_OUT, sep="\t", index=False)

    def null_line(set_name: str, indicator: str) -> str:
        row = gene_null[(gene_null["set"].eq(set_name)) & (gene_null["indicator"].eq(indicator))].iloc[0]
        return (
            f"- {set_name}/{indicator}: {int(row.observed_loci):,} observed loci vs "
            f"{row.null_locus_mean:.1f} null; {row.fold_enrichment:.2f}x; "
            f"p={row.empirical_p_upper:.4g}."
        )

    report = "\n".join(
        [
            f"Open Targets {OPENTARGETS_RELEASE} same-gene colocalisation anchor",
            "============================================================",
            "",
            "Scope",
            "- Uses Open Targets GWAS-molecular QTL colocalisation rows linked to atlas GWAS credible-set coordinate hits.",
            "- Counts candidate genes only when the colocalized molecular-QTL gene matches the GTEx/non-eQTL SCREEN candidate gene.",
            "- This is stronger than L2G/e2G coherence, but still precomputed public genetics support rather than experimental validation.",
            "",
            "Gene-label null",
            null_line("observed_exact_non_eqtl", "ot_coloc_h4_ge_0_8_gene_agree"),
            null_line("moderate_exact_non_eqtl", "ot_coloc_h4_ge_0_8_gene_agree"),
            null_line("lower_recurrence_exact", "ot_coloc_h4_ge_0_8_gene_agree"),
            null_line("systemic_non_neural_exact", "ot_coloc_h4_ge_0_8_gene_agree"),
            null_line("moderate_exact_non_eqtl", "ot_coloc_clpp_ge_0_01_gene_agree"),
            "",
            f"Elapsed seconds: {time.time() - start:.1f}",
            "",
        ]
    )
    OT_COLOC_REPORT_OUT.write_text(report)
    print(f"Report: {display_path(OT_COLOC_REPORT_OUT)}")
    return summary, gene_null, metrics


def run_opentargets_extended_layers() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Add supplementary enhancer-gene, pQTL, and target-context support."""
    start = time.time()
    ensure_opentargets_extended_cache()

    df = pd.read_csv(INPUT, sep="\t", low_memory=False)
    df = normalize_coordinates(df)
    df["qtl_gene_set"] = df["qtl_gene_candidates"].map(gene_set)
    df["screen_non_eqtl_set"] = df["screen_non_eqtl_genes"].map(gene_set)
    df["overlap_gene_set"] = [a & b for a, b in zip(df["qtl_gene_set"], df["screen_non_eqtl_set"])]
    work, _ = prepare_work(df)
    work["variant_key"] = work["CHR"].astype(str) + ":" + work["BP"].astype(str)
    work["overlap_gene_set"] = [a & b for a, b in zip(work["qtl_gene_set"], work["screen_non_eqtl_set"])]
    exact = work[work["has_non_eqtl_overlap"]].copy()

    target_table = ds.dataset(str(OPENTARGETS_CACHE / "target"), format="parquet").to_table(
        columns=["id", "approvedSymbol"]
    ).to_pandas()
    id_to_symbol = dict(zip(target_table["id"], target_table["approvedSymbol"].astype(str).str.upper()))
    symbol_to_id = {symbol: target_id for target_id, symbol in id_to_symbol.items() if symbol and symbol != "NAN"}
    candidate_genes = set().union(*exact["overlap_gene_set"].tolist()) if not exact.empty else set()
    candidate_ids = {symbol_to_id[gene] for gene in candidate_genes if gene in symbol_to_id}

    # Direct interval-level E2G support: candidate variant position falls inside an
    # Open Targets E2G interval linked to the same candidate target gene.
    pos_by_chr: dict[str, list[tuple[int, int, set[str], str, str]]] = {}
    for i, row in exact.iterrows():
        pos_by_chr.setdefault(str(row.CHR), []).append(
            (int(row.BP), int(i), set(row.overlap_gene_set), str(row.locus_500kb), str(row.SNP))
        )
    for chrom in pos_by_chr:
        pos_by_chr[chrom].sort(key=lambda x: x[0])
    pos_lists = {chrom: [entry[0] for entry in entries] for chrom, entries in pos_by_chr.items()}

    e2g_hits: list[dict[str, object]] = []
    e2g_by_work: dict[int, dict[str, object]] = {}
    print("Scanning standalone Open Targets enhancer_to_gene intervals...")
    for file_path in sorted((OPENTARGETS_CACHE / "enhancer_to_gene").glob("*.parquet")):
        pf = pq.ParquetFile(file_path)
        for batch in pf.iter_batches(
            batch_size=100_000,
            columns=[
                "geneId",
                "chromosome",
                "start",
                "end",
                "score",
                "datasourceId",
                "intervalType",
                "biosampleName",
            ],
        ):
            part = batch.to_pandas()
            for row in part.itertuples(index=False):
                chrom = str(row.chromosome)
                positions = pos_lists.get(chrom)
                if not positions:
                    continue
                left = bisect.bisect_left(positions, int(row.start))
                right = bisect.bisect_right(positions, int(row.end))
                if left == right:
                    continue
                symbol = id_to_symbol.get(row.geneId, "")
                if not symbol:
                    continue
                for bp, work_i, genes, locus_id, snp in pos_by_chr[chrom][left:right]:
                    if symbol not in genes:
                        continue
                    score = float(row.score or 0.0)
                    e2g_hits.append(
                        {
                            "variant_key": f"{chrom}:{bp}",
                            "SNP": snp,
                            "locus_500kb": locus_id,
                            "gene_symbol": symbol,
                            "e2g_score": score,
                            "datasourceId": row.datasourceId,
                            "intervalType": row.intervalType,
                            "biosampleName": row.biosampleName,
                        }
                    )
                    entry = e2g_by_work.setdefault(
                        work_i,
                        {
                            "genes": set(),
                            "max_score": 0.0,
                            "biosamples": set(),
                            "interval_types": set(),
                        },
                    )
                    entry["genes"].add(symbol)
                    entry["max_score"] = max(float(entry["max_score"]), score)
                    if isinstance(row.biosampleName, str):
                        entry["biosamples"].add(row.biosampleName)
                    if isinstance(row.intervalType, str):
                        entry["interval_types"].add(row.intervalType)
    pd.DataFrame(e2g_hits).to_csv(OT_E2G_INTERVAL_HITS_OUT, sep="\t", index=False)

    target_support = pd.DataFrame(columns=["gene_symbol", "targetId"])

    coloc_hits = pd.read_csv(OT_COLOC_HITS_OUT, sep="\t")
    pqtl_hits = coloc_hits[coloc_hits["rightStudyType"].astype(str).eq("pqtl")].copy()
    pqtl_h4_pairs = set(
        zip(
            pqtl_hits.loc[pqtl_hits["h4_ge_0_8"].astype(bool), "locus_500kb"].astype(str),
            pqtl_hits.loc[pqtl_hits["h4_ge_0_8"].astype(bool), "gene_symbol"].astype(str),
        )
    )
    pqtl_clpp_pairs = set(
        zip(
            pqtl_hits.loc[pqtl_hits["clpp_ge_0_01"].astype(bool), "locus_500kb"].astype(str),
            pqtl_hits.loc[pqtl_hits["clpp_ge_0_01"].astype(bool), "gene_symbol"].astype(str),
        )
    )

    variant_rows: list[dict[str, object]] = []
    for i, row in work.iterrows():
        genes = set(row.overlap_gene_set)
        locus_id = str(row.locus_500kb)
        e2g = e2g_by_work.get(i, {"genes": set(), "max_score": 0.0, "biosamples": set(), "interval_types": set()})
        pqtl_h4_genes = sorted(gene for gene in genes if (locus_id, gene) in pqtl_h4_pairs)
        pqtl_clpp_genes = sorted(gene for gene in genes if (locus_id, gene) in pqtl_clpp_pairs)
        variant_rows.append(
            {
                "variant_key": row.variant_key,
                "locus_500kb": locus_id,
                "SNP": row.SNP,
                "standalone_e2g_gene_agree": bool(e2g["genes"]),
                "standalone_e2g_gene_symbols": ";".join(sorted(e2g["genes"])),
                "standalone_e2g_max_score": float(e2g["max_score"]),
                "standalone_e2g_biosample_count": len(e2g["biosamples"]),
                "standalone_e2g_interval_types": ";".join(sorted(e2g["interval_types"])),
                "same_gene_pqtl_h4_ge_0_8": bool(pqtl_h4_genes),
                "same_gene_pqtl_h4_genes": ";".join(pqtl_h4_genes),
                "same_gene_pqtl_clpp_ge_0_01": bool(pqtl_clpp_genes),
                "same_gene_pqtl_clpp_genes": ";".join(pqtl_clpp_genes),
            }
        )
    variant_metrics = pd.DataFrame(variant_rows)
    variant_metrics.to_csv(OT_EXTENDED_VARIANT_METRICS_OUT, sep="\t", index=False)

    sets = {
        "observed_exact_non_eqtl": work["has_non_eqtl_overlap"].to_numpy(dtype=bool),
        "moderate_exact_non_eqtl": (work["has_non_eqtl_overlap"] & work["score_band"].eq("moderate_0.2_0.5")).to_numpy(dtype=bool),
        "lower_recurrence_exact": (work["has_non_eqtl_overlap"] & work["source_set"].eq("lower_recurrence_exact")).to_numpy(dtype=bool),
        "core_exact_n10": (work["has_non_eqtl_overlap"] & work["source_set"].eq("core_exact_n10")).to_numpy(dtype=bool),
        "systemic_non_neural_exact": (work["has_non_eqtl_overlap"] & work["primary_axis"].eq("systemic_non_neural")).to_numpy(dtype=bool),
    }
    indicators = [
        "standalone_e2g_gene_agree",
        "same_gene_pqtl_h4_ge_0_8",
        "same_gene_pqtl_clpp_ge_0_01",
    ]
    summary_rows: list[dict[str, object]] = []
    for set_name, mask in sets.items():
        sub_work = work.loc[mask].copy()
        sub_metrics = variant_metrics.loc[mask].copy()
        row: dict[str, object] = {
            "set": set_name,
            "variants": int(len(sub_work)),
            "loci_500kb": int(sub_work["locus_500kb"].nunique()),
        }
        for indicator in indicators:
            indicator_mask = sub_metrics[indicator].to_numpy(dtype=bool)
            row[f"{indicator}_variants"] = int(indicator_mask.sum())
            row[f"{indicator}_loci"] = int(sub_work.loc[indicator_mask, "locus_500kb"].nunique())
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OT_EXTENDED_SUMMARY_OUT, sep="\t", index=False)

    print(
        "Supplementary Open Targets context: "
        f"{int(variant_metrics['standalone_e2g_gene_agree'].sum()):,} variants with direct E2G interval support; "
        f"{int(variant_metrics['same_gene_pqtl_h4_ge_0_8'].sum()):,} variants with same-gene pQTL h4>=0.8."
    )
    print(f"Elapsed seconds: {time.time() - start:.1f}")
    return summary, variant_metrics, target_support


def run_eqtl_catalogue_anchor() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Add eQTL Catalogue r8 beta SuSiE fine-mapped molecular-QTL support."""
    start = time.time()
    ensure_eqtl_catalogue_cache()

    df = pd.read_csv(INPUT, sep="\t", low_memory=False)
    df = normalize_coordinates(df)
    df["qtl_gene_set"] = df["qtl_gene_candidates"].map(gene_set)
    df["screen_non_eqtl_set"] = df["screen_non_eqtl_genes"].map(gene_set)
    df["overlap_gene_set"] = [a & b for a, b in zip(df["qtl_gene_set"], df["screen_non_eqtl_set"])]
    work, _ = prepare_work(df)
    work["variant_key"] = work["CHR"].astype(str) + ":" + work["BP"].astype(str)
    work["overlap_gene_set"] = [a & b for a, b in zip(work["qtl_gene_set"], work["screen_non_eqtl_set"])]

    target_table = ds.dataset(str(OPENTARGETS_CACHE / "target"), format="parquet").to_table(
        columns=["id", "approvedSymbol"]
    ).to_pandas()
    id_to_symbol = dict(zip(target_table["id"].astype(str), target_table["approvedSymbol"].astype(str).str.upper()))

    metadata = pd.read_csv(EQTL_CATALOGUE_METADATA_URL, sep="\t")
    dataset_metadata = metadata.set_index("dataset_id").to_dict("index")

    exact = work[work["has_non_eqtl_overlap"]].copy()
    position_to_work_indices: dict[tuple[str, int], list[int]] = {}
    for i, row in exact.iterrows():
        position_to_work_indices.setdefault((str(row.CHR), int(row.BP)), []).append(i)

    support_sets = {
        "eqtl_catalogue_any_same_gene": [set() for _ in range(len(work))],
        "eqtl_catalogue_pip_ge_0_1_same_gene": [set() for _ in range(len(work))],
        "eqtl_catalogue_pip_ge_0_5_same_gene": [set() for _ in range(len(work))],
    }
    hit_rows: list[dict[str, object]] = []
    files = sorted(EQTL_CATALOGUE_CACHE.glob("*.credible_sets.parquet"))
    print(f"Scanning {len(files)} eQTL Catalogue r8 beta credible-set files...")
    for file_i, file_path in enumerate(files, start=1):
        dataset_id = file_path.name.split(".")[0]
        columns = [
            "chromosome",
            "position",
            "gene_id",
            "molecular_trait_id",
            "rsid",
            "cs_id",
            "cs_size",
            "pip",
            "pvalue",
            "beta",
            "z",
            "cs_min_r2",
            "type",
        ]
        part = pq.read_table(file_path, columns=columns).to_pandas()
        for row in part.itertuples(index=False):
            key = (str(row.chromosome).replace("chr", ""), int(row.position))
            work_indices = position_to_work_indices.get(key)
            if not work_indices:
                continue
            symbol = id_to_symbol.get(str(row.gene_id), str(row.gene_id).upper())
            pip = float(row.pip or 0.0)
            for work_i in work_indices:
                if symbol not in work.at[work_i, "overlap_gene_set"]:
                    continue
                support_sets["eqtl_catalogue_any_same_gene"][work_i].add(symbol)
                if pip >= 0.1:
                    support_sets["eqtl_catalogue_pip_ge_0_1_same_gene"][work_i].add(symbol)
                if pip >= 0.5:
                    support_sets["eqtl_catalogue_pip_ge_0_5_same_gene"][work_i].add(symbol)
                meta = dataset_metadata.get(dataset_id, {})
                hit_rows.append(
                    {
                        "variant_key": work.at[work_i, "variant_key"],
                        "SNP": work.at[work_i, "SNP"],
                        "locus_500kb": work.at[work_i, "locus_500kb"],
                        "gene_symbol": symbol,
                        "dataset_id": dataset_id,
                        "study_id": meta.get("study_id", ""),
                        "study_label": meta.get("study_label", ""),
                        "tissue_label": meta.get("tissue_label", ""),
                        "condition_label": meta.get("condition_label", ""),
                        "quant_method": meta.get("quant_method", ""),
                        "study_type": meta.get("study_type", ""),
                        "pip": pip,
                        "pvalue": float(row.pvalue) if pd.notna(row.pvalue) else "",
                        "beta": float(row.beta) if pd.notna(row.beta) else "",
                        "z": float(row.z) if pd.notna(row.z) else "",
                        "cs_id": row.cs_id,
                        "cs_size": int(row.cs_size) if pd.notna(row.cs_size) else "",
                        "cs_min_r2": float(row.cs_min_r2) if pd.notna(row.cs_min_r2) else "",
                        "rsid": row.rsid,
                        "qtl_type": row.type,
                        "score_band": work.at[work_i, "score_band"],
                        "primary_axis": work.at[work_i, "primary_axis"],
                        "source_set": work.at[work_i, "source_set"],
                    }
                )
        if file_i % 25 == 0:
            print(f"  eQTL Catalogue files {file_i}/{len(files)}; hit rows {len(hit_rows):,}")

    hits = pd.DataFrame(hit_rows).drop_duplicates()
    hits.to_csv(EQTL_CATALOGUE_HITS_OUT, sep="\t", index=False)

    metrics_rows: list[dict[str, object]] = []
    for i, row in work.iterrows():
        any_genes = support_sets["eqtl_catalogue_any_same_gene"][i]
        pip01_genes = support_sets["eqtl_catalogue_pip_ge_0_1_same_gene"][i]
        pip05_genes = support_sets["eqtl_catalogue_pip_ge_0_5_same_gene"][i]
        metrics_rows.append(
            {
                "variant_key": row.variant_key,
                "locus_500kb": row.locus_500kb,
                "eqtl_catalogue_any_same_gene": bool(any_genes),
                "eqtl_catalogue_any_genes": ";".join(sorted(any_genes)),
                "eqtl_catalogue_pip_ge_0_1_same_gene": bool(pip01_genes),
                "eqtl_catalogue_pip_ge_0_1_genes": ";".join(sorted(pip01_genes)),
                "eqtl_catalogue_pip_ge_0_5_same_gene": bool(pip05_genes),
                "eqtl_catalogue_pip_ge_0_5_genes": ";".join(sorted(pip05_genes)),
            }
        )
    metrics = pd.DataFrame(metrics_rows)

    sets = {
        "observed_exact_non_eqtl": work["has_non_eqtl_overlap"].to_numpy(dtype=bool),
        "moderate_exact_non_eqtl": (work["has_non_eqtl_overlap"] & work["score_band"].eq("moderate_0.2_0.5")).to_numpy(dtype=bool),
        "hidden_exact_non_eqtl": (work["has_non_eqtl_overlap"] & work["score_band"].eq("hidden_lt0.2")).to_numpy(dtype=bool),
        "lower_recurrence_exact": (work["has_non_eqtl_overlap"] & work["source_set"].eq("lower_recurrence_exact")).to_numpy(dtype=bool),
        "core_exact_n10": (work["has_non_eqtl_overlap"] & work["source_set"].eq("core_exact_n10")).to_numpy(dtype=bool),
        "systemic_non_neural_exact": (work["has_non_eqtl_overlap"] & work["primary_axis"].eq("systemic_non_neural")).to_numpy(dtype=bool),
    }
    summary_rows: list[dict[str, object]] = []
    for set_name, mask in sets.items():
        sub_work = work.loc[mask].copy()
        sub_metrics = metrics.loc[mask].copy()
        row: dict[str, object] = {
            "set": set_name,
            "variants": int(len(sub_work)),
            "loci_500kb": int(sub_work["locus_500kb"].nunique()),
        }
        for indicator in support_sets:
            indicator_mask = sub_metrics[indicator].to_numpy(dtype=bool)
            row[f"{indicator}_variants"] = int(indicator_mask.sum())
            row[f"{indicator}_loci"] = int(sub_work.loc[indicator_mask, "locus_500kb"].nunique())
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(EQTL_CATALOGUE_SUMMARY_OUT, sep="\t", index=False)

    gene_null = gene_label_permutation_null(work, support_sets, sets)
    gene_null.to_csv(EQTL_CATALOGUE_GENE_NULL_OUT, sep="\t", index=False)

    def null_line(set_name: str, indicator: str) -> str:
        row = gene_null[(gene_null["set"].eq(set_name)) & (gene_null["indicator"].eq(indicator))].iloc[0]
        return (
            f"- {set_name}/{indicator}: {int(row.observed_loci):,} observed loci vs "
            f"{row.null_locus_mean:.1f} null; {row.fold_enrichment:.2f}x; p={row.empirical_p_upper:.4g}."
        )

    report = "\n".join(
        [
            "eQTL Catalogue r8 beta fine-mapped QTL anchor",
            "================================================",
            "",
            "Scope",
            "- Uses eQTL Catalogue r8 beta SuSiE credible-set variants.",
            "- Counts support only when the atlas candidate variant is present in a molecular-QTL credible set for the same candidate gene.",
            "- This is separate from the GTEx/non-eQTL SCREEN convergence layer.",
            "",
            null_line("observed_exact_non_eqtl", "eqtl_catalogue_any_same_gene"),
            null_line("moderate_exact_non_eqtl", "eqtl_catalogue_any_same_gene"),
            null_line("moderate_exact_non_eqtl", "eqtl_catalogue_pip_ge_0_1_same_gene"),
            null_line("moderate_exact_non_eqtl", "eqtl_catalogue_pip_ge_0_5_same_gene"),
            null_line("systemic_non_neural_exact", "eqtl_catalogue_any_same_gene"),
            "",
            f"Elapsed seconds: {time.time() - start:.1f}",
            "",
        ]
    )
    EQTL_CATALOGUE_REPORT_OUT.write_text(report)
    print(f"Report: {display_path(EQTL_CATALOGUE_REPORT_OUT)}")
    return summary, gene_null, metrics


def run_abc_nasser_anchor() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Add ABC enhancer-gene map support from Nasser et al. 2021."""
    start = time.time()
    ensure_abc_nasser_cache()

    df = pd.read_csv(INPUT, sep="\t", low_memory=False)
    df = normalize_coordinates(df)
    df["qtl_gene_set"] = df["qtl_gene_candidates"].map(gene_set)
    df["screen_non_eqtl_set"] = df["screen_non_eqtl_genes"].map(gene_set)
    df["overlap_gene_set"] = [a & b for a, b in zip(df["qtl_gene_set"], df["screen_non_eqtl_set"])]
    work, _ = prepare_work(df)
    work["variant_key"] = work["CHR"].astype(str) + ":" + work["BP"].astype(str)
    work["overlap_gene_set"] = [a & b for a, b in zip(work["qtl_gene_set"], work["screen_non_eqtl_set"])]
    exact = work[work["has_non_eqtl_overlap"]].copy()

    pos_by_chr: dict[str, list[tuple[int, int, set[str], str, str, str, str, str]]] = {}
    for i, row in exact.iterrows():
        pos_by_chr.setdefault(str(row.CHR), []).append(
            (
                int(row.BP),
                int(i),
                set(row.overlap_gene_set),
                str(row.locus_500kb),
                str(row.SNP),
                str(row.score_band),
                str(row.primary_axis),
                str(row.source_set),
            )
        )
    for chrom in pos_by_chr:
        pos_by_chr[chrom].sort(key=lambda x: x[0])
    pos_lists = {chrom: [entry[0] for entry in entries] for chrom, entries in pos_by_chr.items()}

    support_sets = {
        "abc_nasser_any_same_gene": [set() for _ in range(len(work))],
        "abc_nasser_score_ge_0_05_same_gene": [set() for _ in range(len(work))],
        "abc_nasser_score_ge_0_10_same_gene": [set() for _ in range(len(work))],
    }
    max_score_by_work: dict[int, float] = {}
    biosamples_by_work: dict[int, set[str]] = {}
    hit_rows: list[dict[str, object]] = []

    print("Scanning ABC Nasser 2021 enhancer-gene maps...")
    usecols = ["chr", "start", "end", "class", "TargetGene", "ABC.Score", "CellType", "distance", "isSelfPromoter"]
    scanned = 0
    for chunk in pd.read_csv(ABC_NASSER_FILE, sep="\t", compression="gzip", usecols=usecols, chunksize=250_000):
        scanned += len(chunk)
        for chrom, start_i, end_i, klass, gene, score, cell_type, distance, is_self_promoter in zip(
            chunk["chr"],
            chunk["start"],
            chunk["end"],
            chunk["class"],
            chunk["TargetGene"],
            chunk["ABC.Score"],
            chunk["CellType"],
            chunk["distance"],
            chunk["isSelfPromoter"],
        ):
            chrom = str(chrom).replace("chr", "")
            positions = pos_lists.get(chrom)
            if not positions:
                continue
            left = bisect.bisect_left(positions, int(start_i))
            right = bisect.bisect_right(positions, int(end_i))
            if left == right:
                continue
            symbol = str(gene).upper()
            score_f = float(score or 0.0)
            for bp, work_i, genes, locus_id, snp, score_band, axis, source_set in pos_by_chr[chrom][left:right]:
                if symbol not in genes:
                    continue
                support_sets["abc_nasser_any_same_gene"][work_i].add(symbol)
                if score_f >= 0.05:
                    support_sets["abc_nasser_score_ge_0_05_same_gene"][work_i].add(symbol)
                if score_f >= 0.10:
                    support_sets["abc_nasser_score_ge_0_10_same_gene"][work_i].add(symbol)
                max_score_by_work[work_i] = max(max_score_by_work.get(work_i, 0.0), score_f)
                biosamples_by_work.setdefault(work_i, set()).add(str(cell_type))
                hit_rows.append(
                    {
                        "variant_key": f"{chrom}:{bp}",
                        "SNP": snp,
                        "locus_500kb": locus_id,
                        "gene_symbol": symbol,
                        "abc_score": score_f,
                        "cell_type": cell_type,
                        "abc_class": klass,
                        "distance": distance,
                        "is_self_promoter": is_self_promoter,
                        "score_band": score_band,
                        "primary_axis": axis,
                        "source_set": source_set,
                    }
                )
        if scanned % 1_000_000 == 0:
            print(f"  ABC rows scanned {scanned:,}; hit rows {len(hit_rows):,}")

    hits = pd.DataFrame(hit_rows).drop_duplicates()
    hits.to_csv(ABC_NASSER_HITS_OUT, sep="\t", index=False)

    metrics_rows: list[dict[str, object]] = []
    for i, row in work.iterrows():
        any_genes = support_sets["abc_nasser_any_same_gene"][i]
        score05_genes = support_sets["abc_nasser_score_ge_0_05_same_gene"][i]
        score10_genes = support_sets["abc_nasser_score_ge_0_10_same_gene"][i]
        metrics_rows.append(
            {
                "variant_key": row.variant_key,
                "locus_500kb": row.locus_500kb,
                "abc_nasser_any_same_gene": bool(any_genes),
                "abc_nasser_any_genes": ";".join(sorted(any_genes)),
                "abc_nasser_score_ge_0_05_same_gene": bool(score05_genes),
                "abc_nasser_score_ge_0_05_genes": ";".join(sorted(score05_genes)),
                "abc_nasser_score_ge_0_10_same_gene": bool(score10_genes),
                "abc_nasser_score_ge_0_10_genes": ";".join(sorted(score10_genes)),
                "abc_nasser_max_score": max_score_by_work.get(i, 0.0),
                "abc_nasser_biosample_count": len(biosamples_by_work.get(i, set())),
            }
        )
    metrics = pd.DataFrame(metrics_rows)

    sets = {
        "observed_exact_non_eqtl": work["has_non_eqtl_overlap"].to_numpy(dtype=bool),
        "moderate_exact_non_eqtl": (work["has_non_eqtl_overlap"] & work["score_band"].eq("moderate_0.2_0.5")).to_numpy(dtype=bool),
        "hidden_exact_non_eqtl": (work["has_non_eqtl_overlap"] & work["score_band"].eq("hidden_lt0.2")).to_numpy(dtype=bool),
        "lower_recurrence_exact": (work["has_non_eqtl_overlap"] & work["source_set"].eq("lower_recurrence_exact")).to_numpy(dtype=bool),
        "core_exact_n10": (work["has_non_eqtl_overlap"] & work["source_set"].eq("core_exact_n10")).to_numpy(dtype=bool),
        "systemic_non_neural_exact": (work["has_non_eqtl_overlap"] & work["primary_axis"].eq("systemic_non_neural")).to_numpy(dtype=bool),
    }
    summary_rows: list[dict[str, object]] = []
    for set_name, mask in sets.items():
        sub_work = work.loc[mask].copy()
        sub_metrics = metrics.loc[mask].copy()
        row: dict[str, object] = {
            "set": set_name,
            "variants": int(len(sub_work)),
            "loci_500kb": int(sub_work["locus_500kb"].nunique()),
        }
        for indicator in support_sets:
            indicator_mask = sub_metrics[indicator].to_numpy(dtype=bool)
            row[f"{indicator}_variants"] = int(indicator_mask.sum())
            row[f"{indicator}_loci"] = int(sub_work.loc[indicator_mask, "locus_500kb"].nunique())
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(ABC_NASSER_SUMMARY_OUT, sep="\t", index=False)

    gene_null = gene_label_permutation_null(work, support_sets, sets)
    gene_null.to_csv(ABC_NASSER_GENE_NULL_OUT, sep="\t", index=False)

    def null_line(set_name: str, indicator: str) -> str:
        row = gene_null[(gene_null["set"].eq(set_name)) & (gene_null["indicator"].eq(indicator))].iloc[0]
        return (
            f"- {set_name}/{indicator}: {int(row.observed_loci):,} observed loci vs "
            f"{row.null_locus_mean:.1f} null; {row.fold_enrichment:.2f}x; p={row.empirical_p_upper:.4g}."
        )

    report = "\n".join(
        [
            "ABC Nasser 2021 enhancer-gene map anchor",
            "=========================================",
            "",
            "Scope",
            "- Uses the published Nasser et al. 2021 ABC predictions across 131 biosamples.",
            "- Counts support only when the atlas candidate variant falls inside an ABC enhancer interval linked to the same candidate gene.",
            "- This is an independent enhancer-to-gene map, not another SCREEN or GTEx overlap.",
            "",
            null_line("observed_exact_non_eqtl", "abc_nasser_any_same_gene"),
            null_line("moderate_exact_non_eqtl", "abc_nasser_any_same_gene"),
            null_line("moderate_exact_non_eqtl", "abc_nasser_score_ge_0_05_same_gene"),
            null_line("moderate_exact_non_eqtl", "abc_nasser_score_ge_0_10_same_gene"),
            null_line("systemic_non_neural_exact", "abc_nasser_any_same_gene"),
            "",
            f"Elapsed seconds: {time.time() - start:.1f}",
            "",
        ]
    )
    ABC_NASSER_REPORT_OUT.write_text(report)
    print(f"Report: {display_path(ABC_NASSER_REPORT_OUT)}")
    return summary, gene_null, metrics


def run_leave_one_resource_out_target_gene_analysis() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Predict held-out target-gene support using the rest of the evidence stack."""
    start = time.time()

    df = pd.read_csv(INPUT, sep="\t", low_memory=False)
    df = normalize_coordinates(df)
    work, _ = prepare_work(df)
    work["variant_key"] = work["CHR"].astype(str) + ":" + work["BP"].astype(str)
    work["overlap_gene_set"] = [a & b for a, b in zip(work["qtl_gene_set"], work["screen_non_eqtl_set"])]
    work["overlap_3d_gene_set"] = work["overlap_3d_screen_genes"].map(gene_set)
    work["overlap_crispr_gene_set"] = work["overlap_crispr_screen_genes"].map(gene_set)

    target_table = ds.dataset(str(OPENTARGETS_CACHE / "target"), format="parquet").to_table(
        columns=["approvedSymbol", "genomicLocation"]
    ).to_pandas()
    gene_locations: dict[str, list[tuple[str, int]]] = {}
    for row in target_table.itertuples(index=False):
        symbol = str(row.approvedSymbol).upper()
        loc = row.genomicLocation
        if not symbol or symbol == "NAN" or not isinstance(loc, dict):
            continue
        chrom = str(loc.get("chromosome", "")).replace("chr", "")
        if not chrom:
            continue
        start_i = int(loc.get("start") or 0)
        end_i = int(loc.get("end") or start_i)
        strand = int(loc.get("strand") or 1)
        tss = start_i if strand >= 0 else end_i
        gene_locations.setdefault(symbol, []).append((chrom, tss))

    set_masks = {
        "observed_exact_non_eqtl": work["has_non_eqtl_overlap"].to_numpy(dtype=bool),
        "moderate_exact_non_eqtl": (work["has_non_eqtl_overlap"] & work["score_band"].eq("moderate_0.2_0.5")).to_numpy(dtype=bool),
        "lower_recurrence_exact": (work["has_non_eqtl_overlap"] & work["source_set"].eq("lower_recurrence_exact")).to_numpy(dtype=bool),
        "systemic_non_neural_exact": (work["has_non_eqtl_overlap"] & work["primary_axis"].eq("systemic_non_neural")).to_numpy(dtype=bool),
    }

    gene_rows: list[dict[str, object]] = []
    for set_name, mask in set_masks.items():
        for row in work.loc[mask].itertuples(index=False):
            genes = set(row.overlap_gene_set)
            if not genes:
                continue
            genes_3d = set(row.overlap_3d_gene_set)
            genes_crispr = set(row.overlap_crispr_gene_set)
            candidate_gene_count = len(genes)
            for gene in sorted(genes):
                distances = [
                    abs(int(row.BP) - tss)
                    for chrom, tss in gene_locations.get(gene, [])
                    if chrom == str(row.CHR)
                ]
                distance = min(distances) if distances else math.inf
                gene_rows.append(
                    {
                        "set": set_name,
                        "locus_500kb": str(row.locus_500kb),
                        "gene_symbol": gene,
                        "variant_key": str(row.variant_key),
                        "SNP": str(row.SNP),
                        "position": f"{row.CHR}:{int(row.BP)}",
                        "source_set": str(row.source_set),
                        "score_band": str(row.score_band),
                        "primary_axis": str(row.primary_axis),
                        "max_score": float(row.max_score),
                        "n_traits": int(row.n_traits),
                        "n_families": int(row.n_families),
                        "candidate_gene_count": candidate_gene_count,
                        "distance_to_tss": distance,
                        "has_3d_screen": gene in genes_3d,
                        "has_crispr_screen": gene in genes_crispr,
                    }
                )
    if not gene_rows:
        empty = pd.DataFrame()
        empty.to_csv(LEAVE_ONE_OUT_PREDICTIONS_OUT, sep="\t", index=False)
        empty.to_csv(LEAVE_ONE_OUT_SUMMARY_OUT, sep="\t", index=False)
        empty.to_csv(LEAVE_ONE_OUT_NULL_DRAWS_OUT, sep="\t", index=False)
        return empty, empty

    gene_long = pd.DataFrame(gene_rows)
    best_variant = (
        gene_long.sort_values(
            [
                "set",
                "locus_500kb",
                "gene_symbol",
                "max_score",
                "n_families",
                "n_traits",
                "distance_to_tss",
                "candidate_gene_count",
                "variant_key",
            ],
            ascending=[True, True, True, False, False, False, True, True, True],
        )
        .drop_duplicates(["set", "locus_500kb", "gene_symbol"])
        .rename(
            columns={
                "variant_key": "best_variant_key",
                "SNP": "best_snp",
                "position": "best_position",
                "max_score": "best_variant_max_score",
                "n_traits": "best_variant_n_traits",
                "n_families": "best_variant_n_families",
                "distance_to_tss": "best_variant_distance_to_tss",
            }
        )
    )
    best_variant = best_variant[
        [
            "set",
            "locus_500kb",
            "gene_symbol",
            "best_variant_key",
            "best_snp",
            "best_position",
            "best_variant_max_score",
            "best_variant_n_traits",
            "best_variant_n_families",
            "best_variant_distance_to_tss",
        ]
    ]

    candidate_genes = (
        gene_long.groupby(["set", "locus_500kb", "gene_symbol"], as_index=False)
        .agg(
            supporting_variant_count=("variant_key", "nunique"),
            screen_3d_supporting_variant_count=("has_3d_screen", "sum"),
            screen_crispr_supporting_variant_count=("has_crispr_screen", "sum"),
            max_score=("max_score", "max"),
            mean_score=("max_score", "mean"),
            max_n_traits=("n_traits", "max"),
            max_n_families=("n_families", "max"),
            min_candidate_gene_count=("candidate_gene_count", "min"),
            min_distance_to_tss=("distance_to_tss", "min"),
            source_sets=("source_set", lambda s: ";".join(sorted(set(map(str, s))))),
            score_bands=("score_band", lambda s: ";".join(sorted(set(map(str, s))))),
            primary_axes=("primary_axis", lambda s: ";".join(sorted(set(map(str, s))))),
        )
        .merge(best_variant, on=["set", "locus_500kb", "gene_symbol"], how="left")
    )
    candidate_genes["internal_score"] = (
        candidate_genes["supporting_variant_count"].astype(float) * 80.0
        + candidate_genes["screen_crispr_supporting_variant_count"].astype(float) * 30.0
        + candidate_genes["screen_3d_supporting_variant_count"].astype(float) * 20.0
        + candidate_genes["max_n_families"].astype(float) * 5.0
        + np.minimum(candidate_genes["max_n_traits"].astype(float), 200.0) / 10.0
        + candidate_genes["max_score"].astype(float)
        - candidate_genes["min_candidate_gene_count"].astype(float) * 0.25
    )

    variant_evidence = pd.read_csv(VARIANT_EVIDENCE_OUT, sep="\t", low_memory=False)

    def pairs_from_gene_column(column: str) -> set[tuple[str, str]]:
        if column not in variant_evidence.columns:
            return set()
        pairs: set[tuple[str, str]] = set()
        for locus_id, value in zip(variant_evidence["locus_500kb"].astype(str), variant_evidence[column]):
            for gene in gene_set(value):
                pairs.add((locus_id, gene))
        return pairs

    support_pairs = {
        "ot_l2g_score_ge_0_05": pairs_from_gene_column("ot_l2g_score_gene_symbols"),
        "ot_l2g_top_gene": pairs_from_gene_column("ot_l2g_top_gene_symbols"),
        "ot_e2g_feature": pairs_from_gene_column("ot_l2g_e2g_gene_symbols"),
        "standalone_e2g_interval": pairs_from_gene_column("standalone_e2g_gene_symbols"),
        "eqtl_catalogue_any": pairs_from_gene_column("eqtl_catalogue_any_genes"),
        "eqtl_catalogue_pip_ge_0_1": pairs_from_gene_column("eqtl_catalogue_pip_ge_0_1_genes"),
        "eqtl_catalogue_pip_ge_0_5": pairs_from_gene_column("eqtl_catalogue_pip_ge_0_5_genes"),
        "abc_nasser_any": pairs_from_gene_column("abc_nasser_any_genes"),
        "abc_nasser_score_ge_0_05": pairs_from_gene_column("abc_nasser_score_ge_0_05_genes"),
        "pqtl_h4_ge_0_8": pairs_from_gene_column("same_gene_pqtl_h4_genes"),
    }
    coloc_hits = pd.read_csv(OT_COLOC_HITS_OUT, sep="\t")
    support_pairs["ot_coloc_h4_ge_0_8"] = set(
        zip(
            coloc_hits.loc[coloc_hits["h4_ge_0_8"].astype(bool), "locus_500kb"].astype(str),
            coloc_hits.loc[coloc_hits["h4_ge_0_8"].astype(bool), "gene_symbol"].astype(str).str.upper(),
        )
    )
    support_pairs["ot_coloc_clpp_ge_0_01"] = set(
        zip(
            coloc_hits.loc[coloc_hits["clpp_ge_0_01"].astype(bool), "locus_500kb"].astype(str),
            coloc_hits.loc[coloc_hits["clpp_ge_0_01"].astype(bool), "gene_symbol"].astype(str).str.upper(),
        )
    )

    pair_keys = list(zip(candidate_genes["locus_500kb"].astype(str), candidate_genes["gene_symbol"].astype(str)))
    for name, pairs in support_pairs.items():
        candidate_genes[name] = [pair in pairs for pair in pair_keys]

    holdouts = [
        {
            "holdout": "ot_coloc_h4_ge_0_8",
            "outcome": "ot_coloc_h4_ge_0_8",
            "exclude_groups": {"opentargets_colocalisation", "pqtl"},
        },
        {
            "holdout": "eqtl_catalogue_pip_ge_0_1",
            "outcome": "eqtl_catalogue_pip_ge_0_1",
            "exclude_groups": {"eqtl_catalogue"},
        },
        {
            "holdout": "abc_nasser_any",
            "outcome": "abc_nasser_any",
            "exclude_groups": {"abc"},
        },
        {
            "holdout": "ot_l2g_score_ge_0_05",
            "outcome": "ot_l2g_score_ge_0_05",
            "exclude_groups": {"opentargets_l2g_e2g", "standalone_e2g"},
        },
        {
            "holdout": "standalone_e2g_interval",
            "outcome": "standalone_e2g_interval",
            "exclude_groups": {"opentargets_l2g_e2g", "standalone_e2g"},
        },
    ]

    def evidence_score(row: pd.Series, exclude_groups: set[str]) -> float:
        score = float(row["internal_score"])
        if "opentargets_l2g_e2g" not in exclude_groups:
            score += 160.0 * bool(row["ot_l2g_score_ge_0_05"])
            score += 100.0 * bool(row["ot_l2g_top_gene"])
            score += 120.0 * bool(row["ot_e2g_feature"])
        if "opentargets_colocalisation" not in exclude_groups:
            score += 180.0 * bool(row["ot_coloc_h4_ge_0_8"])
            score += 90.0 * bool(row["ot_coloc_clpp_ge_0_01"])
        if "eqtl_catalogue" not in exclude_groups:
            score += 80.0 * bool(row["eqtl_catalogue_any"])
            score += 140.0 * bool(row["eqtl_catalogue_pip_ge_0_1"])
            score += 200.0 * bool(row["eqtl_catalogue_pip_ge_0_5"])
        if "abc" not in exclude_groups:
            score += 110.0 * bool(row["abc_nasser_any"])
            score += 170.0 * bool(row["abc_nasser_score_ge_0_05"])
        if "standalone_e2g" not in exclude_groups:
            score += 80.0 * bool(row["standalone_e2g_interval"])
        if "pqtl" not in exclude_groups:
            score += 80.0 * bool(row["pqtl_h4_ge_0_8"])
        finite_distance = float(row["min_distance_to_tss"])
        if math.isfinite(finite_distance):
            score -= min(finite_distance / 1_000_000.0, 5.0)
        return score

    rng = np.random.default_rng(RANDOM_SEED + 12024)
    summary_rows: list[dict[str, object]] = []
    draw_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []

    for holdout in holdouts:
        holdout_name = holdout["holdout"]
        outcome = holdout["outcome"]
        score_column = f"score_without_{holdout_name}"
        candidate_genes[score_column] = candidate_genes.apply(
            lambda row: evidence_score(row, holdout["exclude_groups"]), axis=1
        )
        for set_name, set_table in candidate_genes.groupby("set", sort=False):
            grouped = []
            for locus_id, sub in set_table.groupby("locus_500kb", sort=False):
                sub = sub.sort_values(
                    [
                        score_column,
                        "internal_score",
                        "supporting_variant_count",
                        "max_n_families",
                        "max_n_traits",
                        "max_score",
                        "min_distance_to_tss",
                        "gene_symbol",
                    ],
                    ascending=[False, False, False, False, False, False, True, True],
                ).copy()
                sub["leave_one_out_rank"] = np.arange(1, len(sub) + 1)
                grouped.append((locus_id, sub))

            n_loci = len(grouped)
            multi_gene_loci = sum(1 for _, sub in grouped if len(sub) > 1)
            top1_observed = 0
            top3_observed = 0
            nearest_observed = 0
            outcome_any = 0

            gene_degree = set_table.groupby("gene_symbol")[outcome].sum().to_dict()
            null_uniform_top1 = np.zeros(N_PERMUTATIONS, dtype=np.int32)
            null_uniform_top3 = np.zeros(N_PERMUTATIONS, dtype=np.int32)
            null_degree_top1 = np.zeros(N_PERMUTATIONS, dtype=np.int32)
            null_distance_top1 = np.zeros(N_PERMUTATIONS, dtype=np.int32)

            values_by_locus: list[np.ndarray] = []
            degree_weights_by_locus: list[np.ndarray] = []
            distance_weights_by_locus: list[np.ndarray] = []
            for locus_id, sub in grouped:
                values = sub[outcome].to_numpy(dtype=bool)
                values_by_locus.append(values)
                degree_weights = np.array([float(gene_degree.get(gene, 0.0)) + 1.0 for gene in sub["gene_symbol"]])
                degree_weights_by_locus.append(degree_weights / degree_weights.sum())
                distances = sub["min_distance_to_tss"].replace([np.inf, -np.inf], np.nan)
                if distances.notna().any():
                    filled = distances.fillna(float(distances.max()) + 1_000_000.0).astype(float).to_numpy()
                    distance_weights = 1.0 / (filled + 10_000.0)
                    distance_weights_by_locus.append(distance_weights / distance_weights.sum())
                else:
                    distance_weights_by_locus.append(np.ones(len(sub), dtype=float) / len(sub))

                top1 = sub.iloc[0]
                top3 = sub.head(3)
                nearest = sub.sort_values(["min_distance_to_tss", "gene_symbol"], ascending=[True, True]).iloc[0]
                top1_observed += int(bool(top1[outcome]))
                top3_observed += int(top3[outcome].any())
                nearest_observed += int(bool(nearest[outcome]))
                outcome_any += int(values.any())
                prediction_rows.append(
                    {
                        "set": set_name,
                        "holdout": holdout_name,
                        "locus_500kb": locus_id,
                        "locus_candidate_gene_count": int(len(sub)),
                        "top1_gene": top1.gene_symbol,
                        "top1_score": float(top1[score_column]),
                        "top1_has_heldout_support": bool(top1[outcome]),
                        "top3_genes": ";".join(top3["gene_symbol"].astype(str).tolist()),
                        "top3_has_heldout_support": bool(top3[outcome].any()),
                        "nearest_gene": nearest.gene_symbol,
                        "nearest_has_heldout_support": bool(nearest[outcome]),
                        "any_candidate_has_heldout_support": bool(values.any()),
                    }
                )

            for p in range(N_PERMUTATIONS):
                uniform_top1 = 0
                uniform_top3 = 0
                degree_top1 = 0
                distance_top1 = 0
                for values, degree_weights, distance_weights in zip(
                    values_by_locus, degree_weights_by_locus, distance_weights_by_locus
                ):
                    if len(values) == 0:
                        continue
                    uniform_top1 += int(values[int(rng.integers(0, len(values)))])
                    k = min(3, len(values))
                    if k == len(values):
                        uniform_top3 += int(values.any())
                    else:
                        uniform_top3 += int(values[rng.choice(len(values), size=k, replace=False)].any())
                    degree_top1 += int(values[int(rng.choice(len(values), p=degree_weights))])
                    distance_top1 += int(values[int(rng.choice(len(values), p=distance_weights))])
                null_uniform_top1[p] = uniform_top1
                null_uniform_top3[p] = uniform_top3
                null_degree_top1[p] = degree_top1
                null_distance_top1[p] = distance_top1

            for mode, observed, null in [
                ("top1_vs_uniform_local_gene_null", top1_observed, null_uniform_top1),
                ("top1_vs_gene_degree_local_gene_null", top1_observed, null_degree_top1),
                ("top1_vs_distance_weighted_local_gene_null", top1_observed, null_distance_top1),
                ("top3_vs_uniform_local_gene_null", top3_observed, null_uniform_top3),
            ]:
                mean = float(null.mean())
                sd = float(null.std(ddof=1))
                summary_rows.append(
                    {
                        "set": set_name,
                        "holdout": holdout_name,
                        "mode": mode,
                        "loci": n_loci,
                        "multi_gene_loci": multi_gene_loci,
                        "heldout_supported_loci_any_candidate": outcome_any,
                        "observed_loci": int(observed),
                        "null_locus_mean": mean,
                        "null_locus_sd": sd,
                        "fold_enrichment": observed / mean if mean > 0 else math.inf,
                        "empirical_p_upper": (float((null >= observed).sum()) + 1.0) / (len(null) + 1.0),
                        "n_permutations": N_PERMUTATIONS,
                    }
                )
                for p, null_value in enumerate(null):
                    draw_rows.append(
                        {
                            "set": set_name,
                            "holdout": holdout_name,
                            "mode": mode,
                            "permutation": p,
                            "null_loci": int(null_value),
                        }
                    )
            summary_rows.append(
                {
                    "set": set_name,
                    "holdout": holdout_name,
                    "mode": "nearest_gene_baseline_descriptive",
                    "loci": n_loci,
                    "multi_gene_loci": multi_gene_loci,
                    "heldout_supported_loci_any_candidate": outcome_any,
                    "observed_loci": int(nearest_observed),
                    "null_locus_mean": "",
                    "null_locus_sd": "",
                    "fold_enrichment": "",
                    "empirical_p_upper": "",
                    "n_permutations": "",
                }
            )

    predictions = pd.DataFrame(prediction_rows)
    summary = pd.DataFrame(summary_rows)
    null_draws = pd.DataFrame(draw_rows)
    predictions.to_csv(LEAVE_ONE_OUT_PREDICTIONS_OUT, sep="\t", index=False)
    summary.to_csv(LEAVE_ONE_OUT_SUMMARY_OUT, sep="\t", index=False)
    null_draws.to_csv(LEAVE_ONE_OUT_NULL_DRAWS_OUT, sep="\t", index=False)

    def line(holdout_name: str, mode: str) -> str:
        row = summary[
            summary["set"].eq("moderate_exact_non_eqtl")
            & summary["holdout"].eq(holdout_name)
            & summary["mode"].eq(mode)
        ].iloc[0]
        return (
            f"| {holdout_name} | {mode} | {int(row.observed_loci):,} observed loci vs "
            f"{float(row.null_locus_mean):.1f} null; {float(row.fold_enrichment):.2f}x; "
            f"p={float(row.empirical_p_upper):.4g} |"
        )

    report = "\n".join(
        [
            "Leave-one-resource-out target-gene prediction",
            "=============================================",
            "",
            "Scope",
            "- Ranks candidate genes per 500 kb locus using the atlas plus all external evidence families except the held-out source.",
            "- Tests whether the held-out evidence source is recovered among top-ranked genes.",
            "- Nulls sample same-locus atlas-nominated genes, with uniform, gene-degree-weighted, and distance-weighted variants.",
            "- This is the integrated downstream-pipeline test for recovering held-out target-gene evidence.",
            "",
            "| held_out_resource | mode | result |",
            "|---|---|---|",
            line("ot_coloc_h4_ge_0_8", "top1_vs_uniform_local_gene_null"),
            line("ot_coloc_h4_ge_0_8", "top1_vs_gene_degree_local_gene_null"),
            line("eqtl_catalogue_pip_ge_0_1", "top1_vs_uniform_local_gene_null"),
            line("eqtl_catalogue_pip_ge_0_1", "top3_vs_uniform_local_gene_null"),
            line("abc_nasser_any", "top1_vs_uniform_local_gene_null"),
            line("ot_l2g_score_ge_0_05", "top1_vs_uniform_local_gene_null"),
            "",
            "Reporting interpretation",
            "- Passing rows mean the integrated atlas pipeline predicts held-out target-gene evidence among same-locus alternatives.",
            "- These rows quantify target-gene triangulation strength under local gene nulls.",
            "",
            f"Elapsed seconds: {time.time() - start:.1f}",
            "",
        ]
    )
    LEAVE_ONE_OUT_REPORT_OUT.write_text(report)
    print(f"Report: {display_path(LEAVE_ONE_OUT_REPORT_OUT)}")

    return summary, predictions


def stable_bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(0).astype(float).ne(0)
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def fit_logistic_terms(
    data: pd.DataFrame,
    outcome: str,
    predictors: list[str],
    model_name: str,
    unit: str,
) -> pd.DataFrame:
    from scipy.optimize import minimize
    from scipy.special import expit
    from scipy.stats import norm

    work = data[[outcome] + predictors].replace([np.inf, -np.inf], np.nan).dropna().copy()
    if work.empty or work[outcome].nunique() < 2:
        return pd.DataFrame()

    y = work[outcome].astype(float).to_numpy()
    x_columns: list[np.ndarray] = []
    term_names: list[str] = []
    term_scales: list[str] = []

    for predictor in predictors:
        values = work[predictor].astype(float).to_numpy()
        if len(np.unique(values)) < 2:
            continue
        unique_values = set(np.unique(values))
        if unique_values.issubset({0.0, 1.0}):
            x_columns.append(values)
            term_names.append(predictor)
            term_scales.append("binary")
        else:
            sd = float(values.std())
            if sd == 0:
                continue
            x_columns.append((values - float(values.mean())) / sd)
            term_names.append(predictor)
            term_scales.append("per_1_sd")

    if not x_columns:
        return pd.DataFrame()

    X = np.column_stack([np.ones(len(work))] + x_columns)

    def loss(beta: np.ndarray) -> float:
        eta = X @ beta
        return float(np.sum(np.logaddexp(0, eta) - y * eta) + 1e-6 * np.sum(beta[1:] ** 2))

    result = minimize(loss, np.zeros(X.shape[1], dtype=float), method="BFGS")
    beta = result.x
    fitted = expit(X @ beta)
    weights = fitted * (1.0 - fitted)
    hessian = X.T @ (weights[:, None] * X)
    hessian += np.eye(hessian.shape[0]) * 1e-6
    hessian[0, 0] -= 1e-6
    covariance = np.linalg.pinv(hessian)
    se = np.sqrt(np.clip(np.diag(covariance), 0, np.inf))

    rows: list[dict[str, object]] = []
    for term, scale, coef, coef_se in zip(term_names, term_scales, beta[1:], se[1:]):
        z_score = coef / coef_se if coef_se > 0 else np.nan
        p_value = float(2.0 * norm.sf(abs(z_score))) if np.isfinite(z_score) else np.nan
        rows.append(
            {
                "model": model_name,
                "unit": unit,
                "outcome": outcome,
                "term": term,
                "scale": scale,
                "n": int(len(work)),
                "events": int(y.sum()),
                "odds_ratio": float(np.exp(coef)),
                "ci95_low": float(np.exp(coef - 1.96 * coef_se)),
                "ci95_high": float(np.exp(coef + 1.96 * coef_se)),
                "wald_p": p_value,
                "converged": bool(result.success),
            }
        )
    return pd.DataFrame(rows)


def run_conditional_recurrence_contribution_model() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(
        INPUT,
        sep="\t",
        usecols=[
            "SNP",
            "position",
            "max_score",
            "qtl_gene_count_for_disjoint",
            "screen_non_eqtl_gene_count",
            "overlap_non_eqtl_screen_count",
            "exact_non_eqtl_screen",
            "eligible_strict_nonblood_hidden_moderate",
        ],
        low_memory=False,
    )
    raw = normalize_coordinates(raw)
    raw["variant_key"] = raw["CHR"].astype(str) + ":" + raw["BP"].astype(str)
    raw["exact_non_eqtl_screen"] = stable_bool_series(raw["exact_non_eqtl_screen"])
    raw["eligible_strict_nonblood_hidden_moderate"] = stable_bool_series(raw["eligible_strict_nonblood_hidden_moderate"])

    variant_evidence = pd.read_csv(VARIANT_EVIDENCE_OUT, sep="\t", low_memory=False)
    variant_evidence = variant_evidence.merge(
        raw[
            [
                "variant_key",
                "max_score",
                "qtl_gene_count_for_disjoint",
                "screen_non_eqtl_gene_count",
                "overlap_non_eqtl_screen_count",
                "exact_non_eqtl_screen",
            ]
        ],
        on="variant_key",
        how="left",
    )

    for column in [
        "has_non_eqtl_overlap",
        "ot_gwas_cs_position_any",
        "ot_gwas_ld_proxy_r2_ge_0_8",
        "ot_l2g_gene_agree_score_ge_0_05",
        "ot_l2g_e2g_gene_agree",
        "ot_coloc_h4_ge_0_8_gene_agree",
        "eqtl_catalogue_pip_ge_0_1_same_gene",
        "abc_nasser_any_same_gene",
    ]:
        variant_evidence[column] = stable_bool_series(variant_evidence[column])

    strict = raw.loc[raw["eligible_strict_nonblood_hidden_moderate"]].copy()
    strict = strict.merge(
        variant_evidence[
            [
                "variant_key",
                "source_set",
                "score_band",
                "primary_axis",
                "n_traits",
                "n_families",
                "has_non_eqtl_overlap",
                "ot_gwas_cs_position_any",
                "ot_gwas_ld_proxy_r2_ge_0_8",
                "ot_l2g_gene_agree_score_ge_0_05",
                "ot_l2g_e2g_gene_agree",
                "ot_coloc_h4_ge_0_8_gene_agree",
                "eqtl_catalogue_pip_ge_0_1_same_gene",
                "abc_nasser_any_same_gene",
            ]
        ],
        on="variant_key",
        how="left",
    )
    strict["has_non_eqtl_overlap"] = strict["has_non_eqtl_overlap"].where(
        strict["has_non_eqtl_overlap"].notna(),
        strict["exact_non_eqtl_screen"],
    )
    strict["has_non_eqtl_overlap"] = stable_bool_series(strict["has_non_eqtl_overlap"])
    for column in [
        "ot_gwas_cs_position_any",
        "ot_gwas_ld_proxy_r2_ge_0_8",
        "ot_l2g_gene_agree_score_ge_0_05",
        "ot_l2g_e2g_gene_agree",
        "ot_coloc_h4_ge_0_8_gene_agree",
        "eqtl_catalogue_pip_ge_0_1_same_gene",
        "abc_nasser_any_same_gene",
    ]:
        strict[column] = stable_bool_series(strict[column])
    strict["source_set"] = strict["source_set"].fillna("unannotated_strict")
    strict["score_band"] = strict["score_band"].fillna("not_permutation_tested")
    strict["primary_axis"] = strict["primary_axis"].fillna("not_permutation_tested")
    strict["n_traits"] = pd.to_numeric(strict["n_traits"], errors="coerce").fillna(0)
    strict["n_families"] = pd.to_numeric(strict["n_families"], errors="coerce").fillna(0)

    rows: list[dict[str, object]] = []
    for locus_id, group in strict.groupby("locus_500kb"):
        row = {
            "locus_500kb": locus_id,
            "strict_variant_count": int(len(group)),
            "exact_variant_count": int(group["has_non_eqtl_overlap"].sum()),
            "has_exact_convergence": int(group["has_non_eqtl_overlap"].any()),
            "moderate_locus": int(group["score_band"].eq("moderate_0.2_0.5").any()),
            "hidden_locus": int(group["score_band"].eq("hidden_lt0.2").any()),
            "high_locus": int(group["score_band"].eq("high_ge0.5").any()),
            "lower_recurrence_locus": int(group["source_set"].eq("lower_recurrence_exact").any()),
            "core_recurrence_locus": int(group["source_set"].eq("core_exact_n10").any()),
            "systemic_non_neural_locus": int(group["primary_axis"].eq("systemic_non_neural").any()),
            "neuro_systemic_locus": int(group["primary_axis"].eq("neuro_systemic").any()),
            "neural_enriched_locus": int(group["primary_axis"].eq("neural_enriched").any()),
            "systemic_other_locus": int(group["primary_axis"].eq("systemic_other").any()),
            "mixed_trait_axis_locus": int(group["primary_axis"].eq("mixed_trait_axis").any()),
            "max_n_families": float(group["n_families"].max()),
            "max_n_traits": float(group["n_traits"].max()),
            "max_score": float(group["max_score"].max()),
            "max_qtl_gene_count": float(group["qtl_gene_count_for_disjoint"].max()),
            "max_screen_gene_count": float(group["screen_non_eqtl_gene_count"].max()),
            "max_candidate_gene_count": float(group["overlap_non_eqtl_screen_count"].max()),
            "ot_coordinate_locus": int(group["ot_gwas_cs_position_any"].any()),
            "ot_ld_proxy_locus": int(group["ot_gwas_ld_proxy_r2_ge_0_8"].any()),
            "ot_l2g_locus": int(group["ot_l2g_gene_agree_score_ge_0_05"].any()),
            "ot_e2g_locus": int(group["ot_l2g_e2g_gene_agree"].any()),
            "ot_coloc_h4_locus": int(group["ot_coloc_h4_ge_0_8_gene_agree"].any()),
            "eqtl_catalogue_pip_ge_0_1_locus": int(group["eqtl_catalogue_pip_ge_0_1_same_gene"].any()),
            "abc_nasser_locus": int(group["abc_nasser_any_same_gene"].any()),
        }
        row["external_any_locus"] = int(
            row["ot_coloc_h4_locus"]
            or row["eqtl_catalogue_pip_ge_0_1_locus"]
            or row["abc_nasser_locus"]
        )
        row["external_multi_source_locus"] = int(
            row["ot_coloc_h4_locus"]
            and (row["eqtl_catalogue_pip_ge_0_1_locus"] or row["abc_nasser_locus"])
        )
        row["full_stack_locus"] = int(
            row["ot_l2g_locus"]
            and row["ot_coloc_h4_locus"]
            and (row["eqtl_catalogue_pip_ge_0_1_locus"] or row["abc_nasser_locus"])
        )
        rows.append(row)

    locus_table = pd.DataFrame(rows)
    locus_table["family_breadth_ge5"] = locus_table["max_n_families"].ge(5).astype(int)
    locus_table["family_breadth_ge6"] = locus_table["max_n_families"].ge(6).astype(int)
    locus_table.to_csv(CONDITIONAL_LOCUS_TABLE_OUT, sep="\t", index=False)

    layer_masks = {
        "all_strict_loci": locus_table.index == locus_table.index,
        "evidence_disjoint_exact_loci": locus_table["has_exact_convergence"].eq(1),
        "moderate_exact_loci": locus_table["has_exact_convergence"].eq(1) & locus_table["moderate_locus"].eq(1),
        "hidden_exact_loci": locus_table["has_exact_convergence"].eq(1) & locus_table["hidden_locus"].eq(1),
        "lower_recurrence_exact_loci": locus_table["has_exact_convergence"].eq(1) & locus_table["lower_recurrence_locus"].eq(1),
        "core_recurrence_exact_loci": locus_table["has_exact_convergence"].eq(1) & locus_table["core_recurrence_locus"].eq(1),
        "systemic_non_neural_exact_loci": locus_table["has_exact_convergence"].eq(1) & locus_table["systemic_non_neural_locus"].eq(1),
    }
    yield_rows: list[dict[str, object]] = []
    for name, mask in layer_masks.items():
        sub = locus_table.loc[mask].copy()
        if sub.empty:
            continue
        yield_rows.append(
            {
                "set": name,
                "loci": int(len(sub)),
                "exact_convergence_loci": int(sub["has_exact_convergence"].sum()),
                "ot_coloc_h4_loci": int(sub["ot_coloc_h4_locus"].sum()),
                "ot_coloc_h4_rate": float(sub["ot_coloc_h4_locus"].mean()),
                "eqtl_catalogue_pip_ge_0_1_loci": int(sub["eqtl_catalogue_pip_ge_0_1_locus"].sum()),
                "eqtl_catalogue_pip_ge_0_1_rate": float(sub["eqtl_catalogue_pip_ge_0_1_locus"].mean()),
                "abc_nasser_loci": int(sub["abc_nasser_locus"].sum()),
                "abc_nasser_rate": float(sub["abc_nasser_locus"].mean()),
                "external_multi_source_loci": int(sub["external_multi_source_locus"].sum()),
                "external_multi_source_rate": float(sub["external_multi_source_locus"].mean()),
                "full_stack_loci": int(sub["full_stack_locus"].sum()),
                "full_stack_rate": float(sub["full_stack_locus"].mean()),
                "median_max_n_families": float(sub["max_n_families"].median()),
                "median_max_score": float(sub["max_score"].median()),
                "median_candidate_gene_count": float(sub["max_candidate_gene_count"].median()),
            }
        )
    layer_yields = pd.DataFrame(yield_rows)
    layer_yields.to_csv(CONDITIONAL_LAYER_YIELDS_OUT, sep="\t", index=False)

    exact = locus_table[locus_table["has_exact_convergence"].eq(1)].copy()
    exact_predictors = [
        "moderate_locus",
        "family_breadth_ge5",
        "max_n_families",
        "max_n_traits",
        "lower_recurrence_locus",
        "systemic_non_neural_locus",
        "strict_variant_count",
        "max_qtl_gene_count",
        "max_screen_gene_count",
        "max_score",
    ]
    external_predictors = exact_predictors + ["max_candidate_gene_count", "ot_ld_proxy_locus"]
    model_frames = [
        fit_logistic_terms(
            locus_table,
            "has_exact_convergence",
            exact_predictors,
            "primary_evidence_disjoint_convergence",
            "strict_eligible_500kb_locus",
        )
    ]
    for outcome, model_name in [
        ("ot_coloc_h4_locus", "external_same_gene_colocalisation_h4"),
        ("eqtl_catalogue_pip_ge_0_1_locus", "external_eqtl_catalogue_pip_ge_0_1"),
        ("abc_nasser_locus", "external_abc_enhancer_gene"),
        ("external_multi_source_locus", "external_multi_source_support"),
        ("full_stack_locus", "external_full_stack_support"),
    ]:
        model_frames.append(
            fit_logistic_terms(
                exact,
                outcome,
                external_predictors,
                model_name,
                "evidence_disjoint_exact_500kb_locus",
            )
        )
    model_results = pd.concat([frame for frame in model_frames if not frame.empty], ignore_index=True)
    model_results.to_csv(CONDITIONAL_MODEL_RESULTS_OUT, sep="\t", index=False)

    def result_line(model: str, term: str) -> str:
        row = model_results[model_results["model"].eq(model) & model_results["term"].eq(term)]
        if row.empty:
            return f"- {model}/{term}: not estimable."
        r = row.iloc[0]
        return (
            f"- {model}/{term}: OR={r.odds_ratio:.2f} "
            f"(95% CI {r.ci95_low:.2f}-{r.ci95_high:.2f}), p={r.wald_p:.4g}."
        )

    report = "\n".join(
        [
            "Conditional atlas-recurrence contribution model",
            "===============================================",
            "",
            "Question",
            "- Does the family-broad moderate recurrence layer add predictive signal after broad annotation and coverage covariates?",
            "",
            "Units",
            "- Primary convergence model: all strict eligible 500 kb loci.",
            "- External support models: 500 kb loci with evidence-disjoint GTEx/non-eQTL SCREEN target-gene convergence.",
            "",
            "Main audit lines",
            result_line("primary_evidence_disjoint_convergence", "moderate_locus"),
            result_line("primary_evidence_disjoint_convergence", "family_breadth_ge5"),
            result_line("external_same_gene_colocalisation_h4", "moderate_locus"),
            result_line("external_same_gene_colocalisation_h4", "family_breadth_ge5"),
            result_line("external_eqtl_catalogue_pip_ge_0_1", "moderate_locus"),
            result_line("external_abc_enhancer_gene", "moderate_locus"),
            "",
            "Interpretation boundary",
            "- This audit tests recurrence as an independent covariate after annotation-rich-locus proxies.",
            "- A weak layer coefficient does not invalidate the atlas; it limits claims that recurrence alone resolves target genes.",
            "- The main atlas claim remains target-gene triangulation, not standalone one-gene causal resolution.",
            "",
            f"Outputs: {display_path(CONDITIONAL_LOCUS_TABLE_OUT)}; {display_path(CONDITIONAL_LAYER_YIELDS_OUT)}; {display_path(CONDITIONAL_MODEL_RESULTS_OUT)}",
            "",
        ]
    )
    CONDITIONAL_REPORT_OUT.write_text(report)
    print(f"Report: {display_path(CONDITIONAL_REPORT_OUT)}")
    return locus_table, layer_yields, model_results


def prepare_mechanism_input() -> pd.DataFrame:
    df = pd.read_csv(INPUT, sep="\t", low_memory=False)
    df = normalize_coordinates(df)
    df["primary_axis"] = df["primary_axis"].replace({"balanced_cross_domain": "mixed_trait_axis"})
    for column in [
        "absent_all_three",
        "eligible_strict_nonblood_hidden_moderate",
        "exact_non_eqtl_screen",
        "exact_3d_screen",
        "exact_crispr_screen",
    ]:
        df[column] = stable_bool_series(df[column])
    df["overlap_gene_set"] = df["overlap_non_eqtl_screen_genes"].map(gene_set)
    df["overlap_gene_count_clean"] = df["overlap_gene_set"].map(len)
    df["single_gene_candidate"] = df["overlap_gene_count_clean"].eq(1)
    df["multi_gene_candidate"] = df["overlap_gene_count_clean"].gt(1)
    df["has_non_eqtl_overlap"] = df["exact_non_eqtl_screen"]
    df["has_3d_overlap"] = df["exact_3d_screen"]
    df["has_crispr_overlap"] = df["exact_crispr_screen"]
    df["catalog_distinct_exact_rsid"] = df["absent_all_three"]
    df["family_breadth_ge5"] = pd.to_numeric(df["n_families"], errors="coerce").fillna(0).ge(5)
    df["family_breadth_ge6"] = pd.to_numeric(df["n_families"], errors="coerce").fillna(0).ge(6)
    return df


def run_mechanism_readiness_layer() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = prepare_mechanism_input()
    exact = df[df["has_non_eqtl_overlap"]].copy()

    taxonomy_rows: list[dict[str, object]] = []
    taxonomy_masks = {
        "all_evidence_disjoint_candidates": pd.Series(True, index=exact.index),
        "single_gene_candidates": exact["single_gene_candidate"],
        "multi_gene_candidates": exact["multi_gene_candidate"],
        "three_dimensional_link_candidates": exact["has_3d_overlap"],
        "exact_rsid_catalog_distinct_candidates": exact["catalog_distinct_exact_rsid"],
        "family_breadth_ge5_candidates": exact["family_breadth_ge5"],
        "family_breadth_ge6_candidates": exact["family_breadth_ge6"],
        "lower_recurrence_candidates": exact["source_set"].eq("lower_recurrence_exact"),
        "core_recurrence_candidates": exact["source_set"].eq("core_exact_n10"),
        "score_hidden_candidates": exact["score_band"].eq("hidden_lt0.2"),
        "score_moderate_candidates": exact["score_band"].eq("moderate_0.2_0.5"),
        "systemic_non_neural_candidates": exact["primary_axis"].eq("systemic_non_neural"),
        "neural_enriched_candidates": exact["primary_axis"].eq("neural_enriched"),
        "neuro_systemic_candidates": exact["primary_axis"].eq("neuro_systemic"),
        "systemic_other_candidates": exact["primary_axis"].eq("systemic_other"),
    }
    for name, mask in taxonomy_masks.items():
        sub = exact[mask]
        if sub.empty:
            continue
        genes = set().union(*sub["overlap_gene_set"].tolist())
        taxonomy_rows.append(
            {
                "set": name,
                "variants": int(len(sub)),
                "broad_loci_500kb": int(sub["locus_500kb"].nunique()),
                "genes": int(len(genes)),
                "single_gene_variants": int(sub["single_gene_candidate"].sum()),
                "multi_gene_variants": int(sub["multi_gene_candidate"].sum()),
                "three_dimensional_link_variants": int(sub["has_3d_overlap"].sum()),
                "crispr_link_variants": int(sub["has_crispr_overlap"].sum()),
                "catalog_distinct_exact_rsid_variants": int(sub["catalog_distinct_exact_rsid"].sum()),
                "median_n_traits": float(sub["n_traits"].median()),
                "median_n_families": float(sub["n_families"].median()),
                "median_max_score": float(sub["max_score"].median()),
                "median_candidate_genes": float(sub["overlap_gene_count_clean"].median()),
            }
        )
    taxonomy = pd.DataFrame(taxonomy_rows)
    taxonomy.to_csv(MECH_DISCOVERY_TAXONOMY_OUT, sep="\t", index=False)

    gene_rows: list[dict[str, object]] = []
    exploded = exact.explode("overlap_gene_set")
    for gene, idx in exploded.groupby("overlap_gene_set").groups.items():
        if not gene:
            continue
        sub = exact.loc[idx]
        axes = sorted(sub["primary_axis"].dropna().astype(str).unique())
        gene_rows.append(
            {
                "gene_symbol": gene,
                "variants": int(len(sub)),
                "broad_loci_500kb": int(sub["locus_500kb"].nunique()),
                "axes": ";".join(axes),
                "axis_count": int(len(axes)),
                "score_bands": ";".join(sorted(sub["score_band"].dropna().astype(str).unique())),
                "source_sets": ";".join(sorted(sub["source_set"].dropna().astype(str).unique())),
                "single_gene_variant_count": int(sub["single_gene_candidate"].sum()),
                "three_dimensional_link_count": int(sub["has_3d_overlap"].sum()),
                "crispr_link_count": int(sub["has_crispr_overlap"].sum()),
                "catalog_distinct_exact_rsid_count": int(sub["catalog_distinct_exact_rsid"].sum()),
                "max_n_traits": int(sub["n_traits"].max()),
                "max_n_families": int(sub["n_families"].max()),
                "median_max_score": float(sub["max_score"].median()),
                "example_variants": ";".join(sub["SNP"].astype(str).head(8)),
            }
        )
    target_genes = pd.DataFrame(gene_rows).sort_values(
        ["broad_loci_500kb", "axis_count", "variants", "single_gene_variant_count", "max_n_families"],
        ascending=[False, False, False, False, False],
    )
    target_genes.to_csv(MECH_TARGET_GENE_RECURRENCE_OUT, sep="\t", index=False)

    locus_rows: list[dict[str, object]] = []
    for locus, sub in exact.groupby("locus_500kb"):
        genes = sorted(set().union(*sub["overlap_gene_set"].tolist()))
        locus_rows.append(
            {
                "locus_500kb": locus,
                "variants": int(len(sub)),
                "genes": ";".join(genes),
                "gene_count": int(len(genes)),
                "single_gene_variant_count": int(sub["single_gene_candidate"].sum()),
                "multi_gene_variant_count": int(sub["multi_gene_candidate"].sum()),
                "three_dimensional_link_count": int(sub["has_3d_overlap"].sum()),
                "crispr_link_count": int(sub["has_crispr_overlap"].sum()),
                "catalog_distinct_exact_rsid_count": int(sub["catalog_distinct_exact_rsid"].sum()),
                "max_n_traits": int(sub["n_traits"].max()),
                "max_n_families": int(sub["n_families"].max()),
                "median_max_score": float(sub["max_score"].median()),
                "primary_axes": ";".join(sorted(sub["primary_axis"].dropna().astype(str).unique())),
                "score_bands": ";".join(sorted(sub["score_band"].dropna().astype(str).unique())),
                "source_sets": ";".join(sorted(sub["source_set"].dropna().astype(str).unique())),
                "clean_single_gene_locus": bool(len(genes) == 1 and sub["single_gene_candidate"].any()),
            }
        )
    locus_ambiguity = pd.DataFrame(locus_rows).sort_values(
        ["clean_single_gene_locus", "max_n_families", "catalog_distinct_exact_rsid_count", "variants"],
        ascending=[False, False, False, False],
    )
    locus_ambiguity.to_csv(MECH_LOCUS_AMBIGUITY_OUT, sep="\t", index=False)

    ready_loci = locus_ambiguity[
        locus_ambiguity["clean_single_gene_locus"]
        & locus_ambiguity["max_n_families"].ge(5)
        & locus_ambiguity["three_dimensional_link_count"].gt(0)
    ].copy()
    ready_rows: list[dict[str, object]] = []
    for row in ready_loci.itertuples(index=False):
        sub = exact[exact["locus_500kb"].eq(row.locus_500kb)].copy()
        best = sub.sort_values(
            ["catalog_distinct_exact_rsid", "has_3d_overlap", "n_families", "n_traits", "max_score"],
            ascending=[False, False, False, False, False],
        ).iloc[0]
        ready_rows.append(
            {
                "locus_500kb": row.locus_500kb,
                "gene_symbol": row.genes,
                "representative_variant": best.SNP,
                "position": best.position,
                "source_set": best.source_set,
                "score_band": best.score_band,
                "primary_axis": best.primary_axis,
                "n_traits": int(best.n_traits),
                "n_families": int(best.n_families),
                "max_score": float(best.max_score),
                "variants_in_locus": int(row.variants),
                "catalog_distinct_exact_rsid_count": int(row.catalog_distinct_exact_rsid_count),
                "three_dimensional_link_count": int(row.three_dimensional_link_count),
                "mechanism_readiness": "clean_single_gene_3d_supported",
            }
        )
    mechanism_ready = pd.DataFrame(ready_rows)
    if not mechanism_ready.empty:
        mechanism_ready = mechanism_ready.sort_values(
            ["n_families", "catalog_distinct_exact_rsid_count", "n_traits", "max_score"],
            ascending=[False, False, False, False],
        )
    mechanism_ready.to_csv(MECH_READY_LOCI_OUT, sep="\t", index=False)

    candidate_rows: list[dict[str, object]] = []
    experiment_candidates = exact[exact["single_gene_candidate"] & exact["family_breadth_ge5"] & exact["has_3d_overlap"]].copy()
    for row in experiment_candidates.itertuples(index=False):
        gene = next(iter(row.overlap_gene_set))
        priority_score = (
            100
            + 12 * int(row.catalog_distinct_exact_rsid)
            + min(int(row.n_families), 20)
            + min(int(row.n_traits) / 5.0, 20.0)
            + min(float(row.max_score) * 20.0, 10.0)
        )
        candidate_rows.append(
            {
                "gene_symbol": gene,
                "SNP": row.SNP,
                "position": row.position,
                "locus_500kb": row.locus_500kb,
                "source_set": row.source_set,
                "score_band": row.score_band,
                "primary_axis": row.primary_axis,
                "n_traits": int(row.n_traits),
                "n_families": int(row.n_families),
                "max_score": float(row.max_score),
                "catalog_distinct_exact_rsid": bool(row.catalog_distinct_exact_rsid),
                "has_3d_overlap": bool(row.has_3d_overlap),
                "has_crispr_overlap": bool(row.has_crispr_overlap),
                "qtl_gene_count": int(row.qtl_gene_count_for_disjoint),
                "screen_non_eqtl_gene_count": int(row.screen_non_eqtl_gene_count),
                "priority_score": float(priority_score),
                "validation_rationale": "single_gene_QTL_non_eqtl_chromatin_convergence",
            }
        )
    experiment_ready = pd.DataFrame(candidate_rows)
    if not experiment_ready.empty:
        experiment_ready = experiment_ready.sort_values(
            ["priority_score", "n_families", "n_traits", "max_score"],
            ascending=[False, False, False, False],
        )
    experiment_ready.to_csv(MECH_EXPERIMENT_CANDIDATES_OUT, sep="\t", index=False)

    all_row = taxonomy[taxonomy["set"].eq("all_evidence_disjoint_candidates")].iloc[0]
    report = f"""# Atlas Mechanism-Readiness Report

## Purpose

This section sits downstream of the public-database evidence stack. It separates clean candidate mechanisms from ambiguous loci and identifies candidates suitable for experimental follow-up.

## Main Results

| result | count |
|---|---:|
| evidence-disjoint candidate variants | {int(all_row.variants):,} |
| evidence-disjoint broad 500 kb loci | {int(all_row.broad_loci_500kb):,} |
| single-gene candidate variants | {int(all_row.single_gene_variants):,} |
| clean single-gene loci | {int(locus_ambiguity["clean_single_gene_locus"].sum()):,} |
| exact-rsID catalog-distinct candidate variants | {int(all_row.catalog_distinct_exact_rsid_variants):,} |
| candidates with CRISPR-linked SCREEN support | {int(all_row.crispr_link_variants):,} |
| mechanism-ready clean single-gene loci | {len(mechanism_ready):,} |
| experiment-ready candidate variants | {len(experiment_ready):,} |

## Interpretation

The atlas yields a practical follow-up layer: hundreds of candidate regulatory variant-gene hypotheses where the GTEx QTL gene and non-eQTL chromatin-linked gene converge on one gene. These loci are suitable for targeted colocalisation, MPRA, CRISPRi, base editing, or other functional validation.

The tables rank loci by regulatory convergence, locus ambiguity, catalog distinctness, and external functional support so follow-up experiments can start from the strongest candidates.
"""
    MECH_REPORT_OUT.write_text(report)
    if REPORT_MD_OUT.exists():
        base_report = REPORT_MD_OUT.read_text().split("\n## Mechanism-Ready Candidate Layer", 1)[0]
        main_report_addition = f"""

## Mechanism-Ready Candidate Layer
| result | count |
|---|---:|
| single-gene candidate variants | {int(all_row.single_gene_variants):,} |
| clean single-gene loci | {int(locus_ambiguity["clean_single_gene_locus"].sum()):,} |
| mechanism-ready clean single-gene loci | {len(mechanism_ready):,} |
| experiment-ready candidate variants | {len(experiment_ready):,} |

These counts define the practical follow-up layer: low-ambiguity candidate regulatory variant-gene hypotheses suitable for targeted colocalisation, MPRA, CRISPRi, base editing, and disease-context follow-up.
"""
        REPORT_MD_OUT.write_text(base_report + main_report_addition)
    print(f"Report: {display_path(MECH_REPORT_OUT)}")
    return taxonomy, target_genes, locus_ambiguity, mechanism_ready, experiment_ready


def download_public_file(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        return
    print(f"Downloading public functional-support resource: {display_path(dest)}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    urllib.request.urlretrieve(url, tmp)
    os.replace(tmp, dest)


def download_mpravardb(dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        return
    if websocket is None:
        raise RuntimeError("MPRAVarDB cache is absent and websocket-client is not installed")
    print(f"Downloading MPRAVarDB MPRA table: {display_path(dest)}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    ws = websocket.create_connection(
        "wss://mpravardb.rc.ufl.edu/websocket/",
        timeout=20,
        header=["Origin: https://mpravardb.rc.ufl.edu"],
    )
    try:
        _ = ws.recv()
        ws.send(
            json.dumps(
                {
                    "method": "init",
                    "data": {
                        "Database_Section": "Single Query",
                        "Study": "Select a study",
                        "Disease": "Select a disease",
                        "CellLine": "Select a cell line",
                        "chr": "",
                        "start": "",
                        "end": "",
                        ".clientdata_url_protocol": "https:",
                        ".clientdata_url_hostname": "mpravardb.rc.ufl.edu",
                        ".clientdata_url_pathname": "/",
                        ".clientdata_url_search": "",
                        ".clientdata_url_hash_initial": "",
                        ".clientdata_url_hash": "",
                        ".clientdata_pixelratio": 1,
                        ".clientdata_singletons": "",
                        ".clientdata_allowDataUriScheme": True,
                        ".clientdata_output_download_all_hidden": False,
                    },
                }
            )
        )
        ws.settimeout(10)
        download_url = None
        for _ in range(30):
            payload = json.loads(ws.recv())
            values = payload.get("values", {})
            if values.get("download_all"):
                download_url = values["download_all"]
                break
        if not download_url:
            raise RuntimeError("MPRAVarDB Shiny session did not return download_all")
        if not str(download_url).startswith("http"):
            download_url = MPRAVARDB_URL + str(download_url).lstrip("/")
        request = urllib.request.Request(
            download_url,
            headers={"User-Agent": "Mozilla/5.0", "Referer": MPRAVARDB_URL},
        )
        tmp = dest.with_suffix(dest.suffix + ".tmp")
        with urllib.request.urlopen(request, timeout=180) as response, tmp.open("wb") as handle:
            handle.write(response.read())
        os.replace(tmp, dest)
    finally:
        ws.close()


def ensure_public_functional_cache() -> None:
    download_public_file(CRISPRI_FLOWFISH_URL, CRISPRI_FLOWFISH_FILE)
    download_public_file(MPRABASE_URL, MPRABASE_DB)
    download_mpravardb(MPRAVARDB_CSV)
    download_public_file(HG38_TO_HG19_URL, HG38_TO_HG19_CHAIN)
    download_public_file(HG19_TO_HG38_URL, HG19_TO_HG38_CHAIN)


def normalize_chr_public(value: object) -> str:
    text = str(value).strip().replace("chr", "")
    if text == "23":
        return "X"
    if text == "24":
        return "Y"
    if text in {"MT", "M"}:
        return "M"
    return text


def parse_public_position(position: object) -> tuple[str, int]:
    chrom, pos = str(position).split(":")
    return normalize_chr_public(chrom), int(float(pos))


def make_liftover(chain: Path):
    if LiftOver is None or not chain.exists():
        return None
    return LiftOver(str(chain))


def lift_one(lo, chrom: str, pos: int) -> tuple[str | None, int | None]:
    if lo is None:
        return None, None
    query = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
    hits = lo.convert_coordinate(query, pos - 1)
    if not hits:
        return None, None
    out_chrom, out_zero, *_ = hits[0]
    return normalize_chr_public(out_chrom), int(out_zero) + 1


def parse_public_interval(text: object, fallback_build: object = "") -> tuple[str | None, int | None, int | None, str | None]:
    if text is None:
        return None, None, None, None
    value = str(text)
    if not value or value.lower() == "nan":
        return None, None, None, None
    pattern = re.compile(r"(?:(hg(?:18|19|38))[:_])?chr?([0-9XYM]+):(\d+)-(\d+)")
    match = pattern.search(value)
    if not match:
        return None, None, None, None
    build = match.group(1) or str(fallback_build or "").strip() or "unknown"
    chrom = normalize_chr_public(match.group(2))
    start = int(match.group(3))
    end = int(match.group(4))
    if start > end:
        start, end = end, start
    return chrom, start, end, build


def first_unique(values: pd.Series | list[object], limit: int = 8) -> str:
    out: list[str] = []
    seen: set[str] = set()
    for value in values:
        if isinstance(value, (list, tuple, set)):
            iterable = value
        else:
            iterable = [value]
        for item in iterable:
            text = "" if item is None else str(item)
            if not text or text.lower() == "nan" or text in seen:
                continue
            seen.add(text)
            out.append(text)
            if len(out) >= limit:
                return ";".join(out)
    return ";".join(out)


def load_public_functional_candidates() -> pd.DataFrame:
    prioritized = pd.read_csv(PRIORITIZED_CANDIDATES_OUT, sep="\t", low_memory=False)
    prioritized["candidate_source"] = "full_same_gene_stack"
    ready = pd.read_csv(MECH_EXPERIMENT_CANDIDATES_OUT, sep="\t", low_memory=False) if MECH_EXPERIMENT_CANDIDATES_OUT.exists() else pd.DataFrame()
    if not ready.empty:
        ready = ready.assign(
            candidate_source="clean_single_gene_followup",
            candidate_gene_count=1,
            same_gene_l2g=False,
            same_gene_top_l2g=False,
            same_gene_e2g=False,
            same_gene_l2g_coloc_feature=False,
            same_gene_ot_coloc_h4_ge_0_8=False,
            same_gene_ot_coloc_clpp_ge_0_01=False,
            same_gene_standalone_e2g_interval=False,
            same_gene_eqtl_catalogue_pip_ge_0_1=False,
            same_gene_eqtl_catalogue_pip_ge_0_5=False,
            same_gene_abc_nasser=False,
            same_gene_abc_nasser_score_ge_0_05=False,
            same_gene_pqtl_h4_ge_0_8=False,
            max_coloc_h4=0.0,
            max_coloc_clpp=0.0,
            n_h4_coloc_study_loci=0,
            n_h4_coloc_qtl_types=0,
            candidate_class="clean_single_gene_followup",
            ot_l2g_max_score_for_candidate_gene=0.0,
            ot_gwas_cs_max_pip=0.0,
            standalone_e2g_max_score=0.0,
            abc_nasser_max_score=0.0,
            abc_nasser_biosample_count=0,
        )
        ready = ready[[col for col in prioritized.columns if col in ready.columns]]
    candidates = pd.concat([prioritized, ready], ignore_index=True, sort=False)
    candidates = candidates.sort_values(["SNP", "gene_symbol", "candidate_source"]).drop_duplicates(["SNP", "gene_symbol"], keep="first")
    coords = candidates["position"].map(parse_public_position)
    candidates["raw_chr"] = [value[0] for value in coords]
    candidates["raw_pos"] = [value[1] for value in coords]
    return candidates.reset_index(drop=True)


def add_public_coordinate_modes(candidates: pd.DataFrame) -> pd.DataFrame:
    hg38_to_hg19 = make_liftover(HG38_TO_HG19_CHAIN)
    hg19_to_hg38 = make_liftover(HG19_TO_HG38_CHAIN)
    if LiftOver is None:
        print("pyliftover is not installed; public functional coordinate tests use as-provided coordinates only")
    rows: list[dict[str, object]] = []
    for row in candidates.itertuples(index=False):
        base = row._asdict()
        rows.append({**base, "coord_mode": "as_provided", "match_chr": row.raw_chr, "match_pos": int(row.raw_pos)})
        chrom, pos = lift_one(hg38_to_hg19, row.raw_chr, int(row.raw_pos))
        if chrom is not None:
            rows.append({**base, "coord_mode": "lift_hg38_to_hg19", "match_chr": chrom, "match_pos": pos})
        chrom, pos = lift_one(hg19_to_hg38, row.raw_chr, int(row.raw_pos))
        if chrom is not None:
            rows.append({**base, "coord_mode": "lift_hg19_to_hg38", "match_chr": chrom, "match_pos": pos})
    return pd.DataFrame(rows)


def public_interval_hits(candidates: pd.DataFrame, intervals: pd.DataFrame, same_gene: bool) -> pd.DataFrame:
    if intervals.empty:
        return pd.DataFrame()
    hits: list[pd.DataFrame] = []
    by_chr = {chrom: group.sort_values("start").reset_index(drop=True) for chrom, group in intervals.groupby("chrom", sort=False)}
    for chrom, cand_group in candidates.groupby("match_chr", sort=False):
        ints = by_chr.get(chrom)
        if ints is None or ints.empty:
            continue
        for cand in cand_group.itertuples(index=False):
            pos = int(cand.match_pos)
            subset = ints[(ints["start"] <= pos) & (ints["end"] >= pos)]
            if same_gene and "gene_symbol" in subset.columns:
                subset = subset[subset["gene_symbol"].astype(str).eq(str(cand.gene_symbol))]
            if subset.empty:
                continue
            left = pd.DataFrame([cand._asdict()] * len(subset)).reset_index(drop=True)
            right = subset.reset_index(drop=True)
            if "gene_symbol" in right.columns and "gene_symbol" in left.columns:
                right = right.rename(columns={"gene_symbol": "perturbed_gene_symbol"})
            hits.append(pd.concat([left, right], axis=1))
    if not hits:
        return pd.DataFrame()
    return pd.concat(hits, ignore_index=True).drop_duplicates()


def run_public_mpravardb(candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    mpra = pd.read_csv(MPRAVARDB_CSV, low_memory=False)
    mpra["chr_norm"] = mpra["chr"].map(normalize_chr_public)
    mpra["pos_int"] = pd.to_numeric(mpra["pos"], errors="coerce").astype("Int64")
    mpra["rsid"] = mpra["rsid"].astype(str)
    mpra["genome"] = mpra["genome"].astype(str)
    base_cols = ["chr", "pos", "ref", "alt", "genome", "rsid", "disease", "cellline", "Description", "log2FC", "pvalue", "fdr", "MPRA_study", "chr_norm", "pos_int"]
    mpra = mpra[base_cols].drop_duplicates()
    exact = candidates.merge(mpra, left_on="SNP", right_on="rsid", how="inner")
    exact["mpravardb_match_type"] = "exact_rsid_mpra_variant"
    coord = candidates.copy()
    coord["match_pos"] = pd.to_numeric(coord["match_pos"], errors="coerce").astype("Int64")
    coord = coord.merge(mpra, left_on=["match_chr", "match_pos"], right_on=["chr_norm", "pos_int"], how="inner")
    valid = (
        (coord["genome"].eq("hg19") & coord["coord_mode"].isin(["as_provided", "lift_hg38_to_hg19"]))
        | (coord["genome"].eq("hg38") & coord["coord_mode"].isin(["as_provided", "lift_hg19_to_hg38"]))
    )
    coord = coord[valid].copy()
    coord["mpravardb_match_type"] = "same_coordinate_mpra_variant"
    hits = pd.concat([exact, coord], ignore_index=True, sort=False).drop_duplicates()
    if hits.empty:
        hits.to_csv(PUBLIC_FUNCTIONAL_MPRAVARDB_HITS_OUT, sep="\t", index=False)
        return hits, pd.DataFrame()
    hits["candidate_rsid_equals_mpra_rsid"] = hits["SNP"].astype(str).eq(hits["rsid"].astype(str))
    hits = hits.sort_values(["candidate_rsid_equals_mpra_rsid", "priority_score"], ascending=[False, False])
    hits.to_csv(PUBLIC_FUNCTIONAL_MPRAVARDB_HITS_OUT, sep="\t", index=False)
    summary = (
        hits.groupby(["SNP", "gene_symbol", "position", "locus_500kb"], dropna=False)
        .agg(
            mpravardb_rows=("SNP", "size"),
            mpravardb_match_types=("mpravardb_match_type", lambda s: ";".join(sorted(set(map(str, s))))),
            mpravardb_rsids=("rsid", lambda s: first_unique(list(s), 10)),
            mpravardb_genomes=("genome", lambda s: first_unique(list(s), 4)),
            mpravardb_diseases=("disease", lambda s: first_unique(list(s), 10)),
            mpravardb_celllines=("cellline", lambda s: first_unique(list(s), 10)),
            mpravardb_studies=("MPRA_study", lambda s: first_unique(list(s), 10)),
            mpravardb_min_fdr=("fdr", lambda s: pd.to_numeric(s, errors="coerce").min()),
            mpravardb_max_abs_log2fc=("log2FC", lambda s: pd.to_numeric(s, errors="coerce").abs().max()),
            mpravardb_exact_rsid_match=("candidate_rsid_equals_mpra_rsid", "max"),
        )
        .reset_index()
    )
    return hits, summary


def load_public_mprabase_intervals() -> pd.DataFrame:
    query = """
        select
            ls.library_element_id,
            ls.genome_build,
            ls.element_coordinate,
            ls.library_element_name,
            s.sample_id,
            s.sample_name,
            s.Library_strategy,
            s.Organism,
            s.Cell_line_tissue,
            d.datasets_name,
            d.PMID,
            d.Reference,
            es.score
        from library_sequence ls
        join sample s on ls.library_id = s.library_id
        left join designed_library dl on ls.library_id = dl.library_id
        left join datasets d on dl.datasets_id = d.datasets_id
        left join element_score es
            on es.library_element_id = ls.library_element_id
            and es.sample_id = s.sample_id
        where s.Organism = 'Homo sapiens'
    """
    with sqlite3.connect(MPRABASE_DB) as conn:
        data = pd.read_sql_query(query, conn)
    rows: list[dict[str, object]] = []
    for row in data.itertuples(index=False):
        parsed = (None, None, None, None)
        for text in [row.element_coordinate, row.library_element_name]:
            parsed = parse_public_interval(text, row.genome_build)
            if parsed[0] is not None:
                break
        chrom, start, end, build = parsed
        if chrom is None:
            continue
        rows.append(
            {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "genome_build": build,
                "library_element_id": row.library_element_id,
                "library_element_name": row.library_element_name,
                "sample_name": row.sample_name,
                "library_strategy": row.Library_strategy,
                "cell_line_tissue": row.Cell_line_tissue,
                "datasets_name": row.datasets_name,
                "pmid": row.PMID,
                "reference": row.Reference,
                "mpra_score": row.score,
            }
        )
    out = pd.DataFrame(rows).drop_duplicates()
    return out[out["genome_build"].isin(["hg19", "hg38"])].copy()


def run_public_mprabase(candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    intervals = load_public_mprabase_intervals()
    hits = public_interval_hits(candidates, intervals, same_gene=False)
    if hits.empty:
        hits.to_csv(PUBLIC_FUNCTIONAL_MPRABASE_HITS_OUT, sep="\t", index=False, compression="gzip")
        return hits, pd.DataFrame()
    valid = (
        (hits["genome_build"].eq("hg19") & hits["coord_mode"].isin(["as_provided", "lift_hg38_to_hg19"]))
        | (hits["genome_build"].eq("hg38") & hits["coord_mode"].isin(["as_provided", "lift_hg19_to_hg38"]))
    )
    hits = hits[valid].copy()
    hits["mprabase_exact_rsid_in_element_name"] = [str(snp) in str(name) for snp, name in zip(hits["SNP"], hits["library_element_name"], strict=False)]
    hits = hits.sort_values(["mprabase_exact_rsid_in_element_name", "priority_score"], ascending=[False, False])
    hits.to_csv(PUBLIC_FUNCTIONAL_MPRABASE_HITS_OUT, sep="\t", index=False, compression="gzip")
    summary = (
        hits.groupby(["SNP", "gene_symbol", "position", "locus_500kb"], dropna=False)
        .agg(
            mprabase_rows=("SNP", "size"),
            mprabase_coordinate_modes=("coord_mode", lambda s: first_unique(list(s), 5)),
            mprabase_genome_builds=("genome_build", lambda s: first_unique(list(s), 5)),
            mprabase_unique_elements=("library_element_id", "nunique"),
            mprabase_samples=("sample_name", "nunique"),
            mprabase_studies=("datasets_name", "nunique"),
            mprabase_cell_line_tissues=("cell_line_tissue", lambda s: first_unique(list(s), 10)),
            mprabase_library_strategies=("library_strategy", lambda s: first_unique(list(s), 10)),
            mprabase_max_abs_score=("mpra_score", lambda s: pd.to_numeric(s, errors="coerce").abs().max()),
            mprabase_exact_rsid_in_element_name=("mprabase_exact_rsid_in_element_name", "max"),
        )
        .reset_index()
    )
    return hits, summary


def run_public_crispri_flowfish(candidates: pd.DataFrame) -> pd.DataFrame:
    data = pd.read_csv(CRISPRI_FLOWFISH_FILE, sep="\t")
    intervals = pd.DataFrame(
        {
            "chrom": data["chrPerturbationTarget"].map(normalize_chr_public),
            "start": data["startPerturbationTarget"].astype(int),
            "end": data["endPerturbationTarget"].astype(int),
            "gene_symbol": data["GeneSymbol"].astype(str),
            "crispri_effect_size": pd.to_numeric(data["EffectSize"], errors="coerce"),
            "crispri_significant": data["Significant"].astype(str),
            "crispri_regulated": data["Regulated"].astype(str),
            "crispri_include_in_model": data["IncludeInModel"].astype(str),
            "crispri_cell_type": data["CellType"].astype(str),
            "crispri_padj": pd.to_numeric(data["padj"], errors="coerce"),
            "crispri_reference": data["Reference"].astype(str),
        }
    )
    starts = intervals[["start", "end"]].min(axis=1)
    ends = intervals[["start", "end"]].max(axis=1)
    intervals["start"], intervals["end"] = starts, ends
    hits = public_interval_hits(candidates, intervals, same_gene=True)
    if not hits.empty:
        hits["direct_crispri_flowfish_support"] = hits["crispri_significant"].str.upper().eq("TRUE") & hits["crispri_regulated"].str.upper().eq("TRUE")
        hits = hits.sort_values(["direct_crispri_flowfish_support", "priority_score", "crispri_padj"], ascending=[False, False, True])
    hits.to_csv(PUBLIC_FUNCTIONAL_CRISPRI_HITS_OUT, sep="\t", index=False)
    return hits


def run_public_opentargets_crispr(candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    parts: list[pd.DataFrame] = []
    if OT_CRISPR_CACHE.exists():
        parts.append(pd.read_parquet(OT_CRISPR_CACHE).assign(ot_crispr_source="crispr"))
    if OT_CRISPR_SCREEN_CACHE.exists():
        parts.append(pd.read_parquet(OT_CRISPR_SCREEN_CACHE).assign(ot_crispr_source="crispr_screen"))
    if not parts:
        return pd.DataFrame(), pd.DataFrame()
    evidence = pd.concat(parts, ignore_index=True)
    evidence["target_symbol"] = evidence["targetFromSource"].astype(str)
    cols = [c for c in ["target_symbol", "diseaseFromSource", "diseaseId", "datatypeId", "datasourceId", "ot_crispr_source", "resourceScore", "score", "literature"] if c in evidence.columns]
    evidence = evidence[cols].copy()
    for col in evidence.columns:
        evidence[col] = evidence[col].map(lambda value: "|".join(map(str, value)) if isinstance(value, (list, tuple)) else str(value))
    hits = candidates.drop_duplicates(["SNP", "gene_symbol"]).merge(evidence.drop_duplicates(), left_on="gene_symbol", right_on="target_symbol", how="inner")
    hits = hits.sort_values(["priority_score", "score"], ascending=[False, False])
    hits.to_csv(PUBLIC_FUNCTIONAL_OT_CRISPR_HITS_OUT, sep="\t", index=False)
    summary = (
        hits.groupby(["SNP", "gene_symbol", "position", "locus_500kb"], dropna=False)
        .agg(
            ot_crispr_context_rows=("SNP", "size"),
            ot_crispr_max_score=("score", lambda s: pd.to_numeric(s, errors="coerce").max()),
            ot_crispr_diseases=("diseaseFromSource", lambda s: first_unique(list(s), 8)),
        )
        .reset_index()
    )
    return hits, summary


def run_public_functional_support_layer() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ensure_public_functional_cache()
    base = load_public_functional_candidates()
    coordinate_candidates = add_public_coordinate_modes(base)
    crispri_hits = run_public_crispri_flowfish(coordinate_candidates)
    _, mpravardb_summary = run_public_mpravardb(coordinate_candidates)
    _, mprabase_summary = run_public_mprabase(coordinate_candidates)
    _, ot_crispr_summary = run_public_opentargets_crispr(base)

    support = base.copy()
    for table in [mpravardb_summary, mprabase_summary, ot_crispr_summary]:
        if not table.empty:
            support = support.merge(table, on=["SNP", "gene_symbol", "position", "locus_500kb"], how="left")
    if not crispri_hits.empty:
        crispri_summary = (
            crispri_hits.groupby(["SNP", "gene_symbol", "position", "locus_500kb"], dropna=False)
            .agg(
                crispri_flowfish_rows=("SNP", "size"),
                direct_crispri_flowfish_support=("direct_crispri_flowfish_support", "max"),
                crispri_flowfish_min_padj=("crispri_padj", "min"),
                crispri_flowfish_cell_types=("crispri_cell_type", lambda s: first_unique(list(s), 8)),
            )
            .reset_index()
        )
        support = support.merge(crispri_summary, on=["SNP", "gene_symbol", "position", "locus_500kb"], how="left")

    for col in [
        "mpravardb_rows",
        "mprabase_rows",
        "mprabase_unique_elements",
        "ot_crispr_context_rows",
        "crispri_flowfish_rows",
    ]:
        values = support[col] if col in support.columns else pd.Series(0, index=support.index)
        support[col] = pd.to_numeric(values, errors="coerce").fillna(0).astype(int)
    support["mpravardb_exact_rsid_match"] = stable_bool_series(support.get("mpravardb_exact_rsid_match", pd.Series(False, index=support.index)))
    support["mpravardb_significant"] = pd.to_numeric(
        support["mpravardb_min_fdr"] if "mpravardb_min_fdr" in support.columns else pd.Series(np.nan, index=support.index),
        errors="coerce",
    ).lt(0.05)
    support["mprabase_element_overlap"] = support["mprabase_unique_elements"].gt(0)
    support["ot_crispr_gene_context"] = support["ot_crispr_context_rows"].gt(0)
    support["direct_crispri_flowfish_support"] = stable_bool_series(support.get("direct_crispri_flowfish_support", pd.Series(False, index=support.index)))
    support["any_public_functional_support"] = (
        support["mpravardb_exact_rsid_match"]
        | support["mpravardb_significant"]
        | support["mprabase_element_overlap"]
        | support["ot_crispr_gene_context"]
        | support["direct_crispri_flowfish_support"]
    )
    support.to_csv(MECH_PUBLIC_FUNCTIONAL_SUPPORT_OUT, sep="\t", index=False)

    rows: list[dict[str, object]] = []
    def add_summary(name: str, mask: pd.Series, interpretation: str) -> None:
        sub = support[mask]
        rows.append(
            {
                "result": name,
                "candidate_rows": int(len(sub)),
                "variants": int(sub["SNP"].nunique()) if not sub.empty else 0,
                "broad_loci_500kb": int(sub["locus_500kb"].nunique()) if not sub.empty else 0,
                "genes": int(sub["gene_symbol"].nunique()) if not sub.empty else 0,
                "interpretation": interpretation,
            }
        )
    add_summary("input_candidate_pairs", pd.Series(True, index=support.index), "Union of full-stack candidates and clean single-gene follow-up candidates.")
    add_summary("mpravardb_exact_rsid_variant", support["mpravardb_exact_rsid_match"], "Candidate rsID has public MPRA variant-effect data.")
    add_summary("mpravardb_significant_exact_rsid_variant", support["mpravardb_significant"], "Candidate rsID has FDR<0.05 public MPRA variant-effect evidence.")
    add_summary("mprabase_element_overlap", support["mprabase_element_overlap"], "Candidate coordinate falls inside a public MPRA-assayed element.")
    add_summary("opentargets_crispr_gene_context", support["ot_crispr_gene_context"], "Candidate gene has Open Targets CRISPR or CRISPR-screen disease/cell context.")
    add_summary("direct_crispri_flowfish_same_gene", support["direct_crispri_flowfish_support"], "Candidate coordinate overlaps a CRISPRi-FlowFISH perturbation interval for the same gene.")
    summary = pd.DataFrame(rows)
    summary.to_csv(MECH_PUBLIC_FUNCTIONAL_SUMMARY_OUT, sep="\t", index=False)

    ranked = support.copy()
    ranked["mechanism_priority_score"] = pd.to_numeric(ranked.get("priority_score"), errors="coerce").fillna(0.0)
    ranked["mechanism_priority_score"] += 20 * stable_bool_series(ranked.get("same_gene_ot_coloc_h4_ge_0_8", pd.Series(False, index=ranked.index))).astype(int)
    ranked["mechanism_priority_score"] += 14 * stable_bool_series(ranked.get("same_gene_eqtl_catalogue_pip_ge_0_5", pd.Series(False, index=ranked.index))).astype(int)
    ranked["mechanism_priority_score"] += 8 * stable_bool_series(ranked.get("same_gene_eqtl_catalogue_pip_ge_0_1", pd.Series(False, index=ranked.index))).astype(int)
    ranked["mechanism_priority_score"] += 10 * stable_bool_series(ranked.get("same_gene_abc_nasser_score_ge_0_05", pd.Series(False, index=ranked.index))).astype(int)
    ranked["mechanism_priority_score"] += 8 * stable_bool_series(ranked.get("same_gene_pqtl_h4_ge_0_8", pd.Series(False, index=ranked.index))).astype(int)
    ranked["mechanism_priority_score"] += 18 * ranked["mpravardb_significant"].astype(int)
    ranked["mechanism_priority_score"] += 6 * ranked["mpravardb_exact_rsid_match"].astype(int)
    ranked["mechanism_priority_score"] += 5 * ranked["mprabase_element_overlap"].astype(int)
    ranked["mechanism_priority_score"] += 4 * ranked["ot_crispr_gene_context"].astype(int)
    ranked["mechanism_priority_score"] += 12 * ranked["direct_crispri_flowfish_support"].astype(int)
    ranked["mechanism_priority_score"] += 5 * pd.to_numeric(ranked.get("candidate_gene_count"), errors="coerce").fillna(2).le(1).astype(int)
    ranked["mechanism_priority_score"] += pd.to_numeric(ranked.get("n_families"), errors="coerce").fillna(0).clip(upper=10)
    ranked["mechanism_priority_score"] += (pd.to_numeric(ranked.get("n_traits"), errors="coerce").fillna(0) / 5.0).clip(upper=10)
    ranked["mechanism_support_class"] = "multi_source_candidate"
    ranked.loc[ranked["mpravardb_significant"] | ranked["direct_crispri_flowfish_support"], "mechanism_support_class"] = "public_functional_support_candidate"
    ranked.loc[
        stable_bool_series(ranked.get("same_gene_ot_coloc_h4_ge_0_8", pd.Series(False, index=ranked.index)))
        & (
            stable_bool_series(ranked.get("same_gene_eqtl_catalogue_pip_ge_0_1", pd.Series(False, index=ranked.index)))
            | stable_bool_series(ranked.get("same_gene_abc_nasser", pd.Series(False, index=ranked.index)))
            | stable_bool_series(ranked.get("same_gene_pqtl_h4_ge_0_8", pd.Series(False, index=ranked.index)))
        ),
        "mechanism_support_class",
    ] = "triangulated_mechanism_anchor"
    ranked = ranked.sort_values(["mechanism_priority_score", "priority_score", "n_families", "n_traits"], ascending=[False, False, False, False])
    ranked.to_csv(MECH_PRIORITIZED_CANDIDATES_OUT, sep="\t", index=False)

    dossier_cols = [
        "locus_500kb",
        "SNP",
        "gene_symbol",
        "mechanism_priority_score",
        "mechanism_support_class",
        "primary_axis",
        "n_traits",
        "n_families",
        "candidate_gene_count",
        "same_gene_ot_coloc_h4_ge_0_8",
        "same_gene_eqtl_catalogue_pip_ge_0_1",
        "same_gene_eqtl_catalogue_pip_ge_0_5",
        "same_gene_abc_nasser",
        "same_gene_pqtl_h4_ge_0_8",
        "mpravardb_significant",
        "mpravardb_min_fdr",
        "mpravardb_max_abs_log2fc",
        "mprabase_element_overlap",
        "ot_crispr_gene_context",
    ]
    top = ranked[[col for col in dossier_cols if col in ranked.columns]].head(30)

    def display_label(value: object) -> str:
        labels = {
            "neural_enriched": "neural-enriched",
            "systemic_non_neural": "systemic non-neural",
            "systemic_other": "systemic",
            "mixed_trait_axis": "mixed trait axis",
            "neuro_systemic": "neuro-systemic",
            "triangulated_mechanism_anchor": "triangulated mechanism anchor",
            "multi_source_candidate": "multi-source candidate",
        }
        text = "" if value is None else str(value)
        return labels.get(text, text.replace("_", " "))

    lines = ["# Prioritized Candidate Mechanisms", "", "Ranked by atlas recurrence, same-gene external evidence, public functional support, and locus ambiguity.", ""]
    for i, row in enumerate(top.itertuples(index=False), start=1):
        data = row._asdict()
        evidence = []
        if data.get("same_gene_ot_coloc_h4_ge_0_8"):
            evidence.append("same-gene GWAS-molQTL colocalisation")
        if data.get("same_gene_eqtl_catalogue_pip_ge_0_5"):
            evidence.append("eQTL Catalogue PIP>=0.5")
        elif data.get("same_gene_eqtl_catalogue_pip_ge_0_1"):
            evidence.append("eQTL Catalogue PIP>=0.1")
        if data.get("same_gene_abc_nasser"):
            evidence.append("ABC enhancer-gene support")
        if data.get("same_gene_pqtl_h4_ge_0_8"):
            evidence.append("pQTL colocalisation")
        if data.get("mpravardb_significant"):
            evidence.append("FDR-significant MPRA variant")
        if data.get("mprabase_element_overlap"):
            evidence.append("MPRA element overlap")
        if data.get("ot_crispr_gene_context"):
            evidence.append("CRISPR gene context")
        lines.extend(
            [
                f"## {i}. {data.get('gene_symbol')} / {data.get('SNP')} / {data.get('locus_500kb')}",
                "",
                f"- Priority score: {float(data.get('mechanism_priority_score', 0)):.1f} ({display_label(data.get('mechanism_support_class'))})",
                f"- Atlas recurrence: {int(data.get('n_traits', 0))} traits, {int(data.get('n_families', 0))} families, axis {display_label(data.get('primary_axis'))}",
                f"- Candidate ambiguity: {int(float(data.get('candidate_gene_count', 0) or 0))} candidate gene(s)",
                f"- Evidence: {'; '.join(evidence) if evidence else 'multi-source candidate evidence'}",
                "",
            ]
        )
    MECH_CANDIDATE_DOSSIERS_OUT.write_text("\n".join(lines), encoding="utf-8")
    mpra_sig = summary[summary["result"].eq("mpravardb_significant_exact_rsid_variant")].iloc[0]
    mpra_any = summary[summary["result"].eq("mpravardb_exact_rsid_variant")].iloc[0]
    mprabase = summary[summary["result"].eq("mprabase_element_overlap")].iloc[0]
    crispr_gene = summary[summary["result"].eq("opentargets_crispr_gene_context")].iloc[0]
    functional_report = f"""

## Public Functional-Support Layer

| result | count |
|---|---:|
| candidate variant-gene rows tested | {len(support):,} |
| candidate variants with exact-rsID MPRAVarDB testing | {int(mpra_any.variants):,} |
| FDR<0.05 MPRAVarDB candidate rows | {int(mpra_sig.candidate_rows):,} |
| MPRAbase element-overlap candidate rows | {int(mprabase.candidate_rows):,} |
| Open Targets CRISPR gene-context candidate rows | {int(crispr_gene.candidate_rows):,} |
| ranked candidate mechanisms | {len(ranked):,} |

Public functional genomics resources add orthogonal prioritization evidence from MPRA, enhancer-element assays, CRISPR screens, and disease-cell-context annotations. Candidate rankings are generated directly from atlas recurrence and external evidence features.
"""
    if MECH_REPORT_OUT.exists():
        MECH_REPORT_OUT.write_text(MECH_REPORT_OUT.read_text() + functional_report)
    if REPORT_MD_OUT.exists():
        REPORT_MD_OUT.write_text(REPORT_MD_OUT.read_text() + functional_report)
    print(f"Public functional support: {display_path(MECH_PUBLIC_FUNCTIONAL_SUMMARY_OUT)}")
    print(f"Prioritized candidate mechanisms: {display_path(MECH_PRIORITIZED_CANDIDATES_OUT)}")
    return support, summary, ranked


def write_clean_main_outputs(
    locus_summary: pd.DataFrame,
    denominators: pd.DataFrame,
    ot_summary: pd.DataFrame,
    ot_gene_null: pd.DataFrame,
    coloc_summary: pd.DataFrame,
    coloc_gene_null: pd.DataFrame,
    coloc_metrics: pd.DataFrame,
    extended_summary: pd.DataFrame,
    extended_metrics: pd.DataFrame,
    eqtl_summary: pd.DataFrame,
    eqtl_gene_null: pd.DataFrame,
    eqtl_metrics: pd.DataFrame,
    abc_summary: pd.DataFrame,
    abc_gene_null: pd.DataFrame,
    abc_metrics: pd.DataFrame,
) -> None:
    def locus_row(name: str) -> pd.Series:
        return locus_summary[locus_summary["set"].eq(name)].iloc[0]

    def ot_row(name: str) -> pd.Series:
        return ot_summary[ot_summary["set"].eq(name)].iloc[0]

    def coloc_row(name: str) -> pd.Series:
        return coloc_summary[coloc_summary["set"].eq(name)].iloc[0]

    def extended_row(name: str) -> pd.Series:
        return extended_summary[extended_summary["set"].eq(name)].iloc[0]

    def eqtl_row(name: str) -> pd.Series:
        return eqtl_summary[eqtl_summary["set"].eq(name)].iloc[0]

    def abc_row(name: str) -> pd.Series:
        return abc_summary[abc_summary["set"].eq(name)].iloc[0]

    def null_row(table: pd.DataFrame, set_name: str, indicator: str) -> pd.Series:
        return table[(table["set"].eq(set_name)) & (table["indicator"].eq(indicator))].iloc[0]

    novelty_summary = pd.read_csv(OT_NOVELTY_SUMMARY_OUT, sep="\t")
    def novelty_row(name: str) -> pd.Series:
        return novelty_summary[novelty_summary["set"].eq(name)].iloc[0]

    main_rows = [
        {
            "result": "evidence_disjoint_locus_convergence_overall",
            "observed_loci": int(locus_row("overall_non_eqtl_tested").observed_concordant_loci_500kb),
            "null_loci": float(locus_row("overall_non_eqtl_tested").null_locus_mean),
            "fold": float(locus_row("overall_non_eqtl_tested").locus_fold_enrichment),
            "empirical_p": float(locus_row("overall_non_eqtl_tested").locus_empirical_p_upper),
        },
        {
            "result": "evidence_disjoint_locus_convergence_moderate",
            "observed_loci": int(locus_row("score_moderate_0.2_0.5").observed_concordant_loci_500kb),
            "null_loci": float(locus_row("score_moderate_0.2_0.5").null_locus_mean),
            "fold": float(locus_row("score_moderate_0.2_0.5").locus_fold_enrichment),
            "empirical_p": float(locus_row("score_moderate_0.2_0.5").locus_empirical_p_upper),
        },
        {
            "result": "opentargets_l2g_agreement_moderate",
            "observed_loci": int(null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_gene_agree_score_ge_0_05").observed_loci),
            "null_loci": float(null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_gene_agree_score_ge_0_05").null_locus_mean),
            "fold": float(null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_gene_agree_score_ge_0_05").fold_enrichment),
            "empirical_p": float(null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_gene_agree_score_ge_0_05").empirical_p_upper),
        },
        {
            "result": "opentargets_e2g_agreement_moderate",
            "observed_loci": int(null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_e2g_gene_agree").observed_loci),
            "null_loci": float(null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_e2g_gene_agree").null_locus_mean),
            "fold": float(null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_e2g_gene_agree").fold_enrichment),
            "empirical_p": float(null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_e2g_gene_agree").empirical_p_upper),
        },
        {
            "result": "opentargets_colocalisation_h4_ge_0_8_moderate",
            "observed_loci": int(null_row(coloc_gene_null, "moderate_exact_non_eqtl", "ot_coloc_h4_ge_0_8_gene_agree").observed_loci),
            "null_loci": float(null_row(coloc_gene_null, "moderate_exact_non_eqtl", "ot_coloc_h4_ge_0_8_gene_agree").null_locus_mean),
            "fold": float(null_row(coloc_gene_null, "moderate_exact_non_eqtl", "ot_coloc_h4_ge_0_8_gene_agree").fold_enrichment),
            "empirical_p": float(null_row(coloc_gene_null, "moderate_exact_non_eqtl", "ot_coloc_h4_ge_0_8_gene_agree").empirical_p_upper),
        },
        {
            "result": "opentargets_colocalisation_clpp_ge_0_01_moderate",
            "observed_loci": int(null_row(coloc_gene_null, "moderate_exact_non_eqtl", "ot_coloc_clpp_ge_0_01_gene_agree").observed_loci),
            "null_loci": float(null_row(coloc_gene_null, "moderate_exact_non_eqtl", "ot_coloc_clpp_ge_0_01_gene_agree").null_locus_mean),
            "fold": float(null_row(coloc_gene_null, "moderate_exact_non_eqtl", "ot_coloc_clpp_ge_0_01_gene_agree").fold_enrichment),
            "empirical_p": float(null_row(coloc_gene_null, "moderate_exact_non_eqtl", "ot_coloc_clpp_ge_0_01_gene_agree").empirical_p_upper),
        },
        {
            "result": "standalone_e2g_interval_support_moderate",
            "observed_loci": int(extended_row("moderate_exact_non_eqtl").standalone_e2g_gene_agree_loci),
            "null_loci": "",
            "fold": "",
            "empirical_p": "",
        },
        {
            "result": "eqtl_catalogue_pip_ge_0_1_same_gene_moderate",
            "observed_loci": int(null_row(eqtl_gene_null, "moderate_exact_non_eqtl", "eqtl_catalogue_pip_ge_0_1_same_gene").observed_loci),
            "null_loci": float(null_row(eqtl_gene_null, "moderate_exact_non_eqtl", "eqtl_catalogue_pip_ge_0_1_same_gene").null_locus_mean),
            "fold": float(null_row(eqtl_gene_null, "moderate_exact_non_eqtl", "eqtl_catalogue_pip_ge_0_1_same_gene").fold_enrichment),
            "empirical_p": float(null_row(eqtl_gene_null, "moderate_exact_non_eqtl", "eqtl_catalogue_pip_ge_0_1_same_gene").empirical_p_upper),
        },
        {
            "result": "eqtl_catalogue_pip_ge_0_5_same_gene_moderate",
            "observed_loci": int(null_row(eqtl_gene_null, "moderate_exact_non_eqtl", "eqtl_catalogue_pip_ge_0_5_same_gene").observed_loci),
            "null_loci": float(null_row(eqtl_gene_null, "moderate_exact_non_eqtl", "eqtl_catalogue_pip_ge_0_5_same_gene").null_locus_mean),
            "fold": float(null_row(eqtl_gene_null, "moderate_exact_non_eqtl", "eqtl_catalogue_pip_ge_0_5_same_gene").fold_enrichment),
            "empirical_p": float(null_row(eqtl_gene_null, "moderate_exact_non_eqtl", "eqtl_catalogue_pip_ge_0_5_same_gene").empirical_p_upper),
        },
        {
            "result": "abc_nasser_same_gene_moderate",
            "observed_loci": int(null_row(abc_gene_null, "moderate_exact_non_eqtl", "abc_nasser_any_same_gene").observed_loci),
            "null_loci": float(null_row(abc_gene_null, "moderate_exact_non_eqtl", "abc_nasser_any_same_gene").null_locus_mean),
            "fold": float(null_row(abc_gene_null, "moderate_exact_non_eqtl", "abc_nasser_any_same_gene").fold_enrichment),
            "empirical_p": float(null_row(abc_gene_null, "moderate_exact_non_eqtl", "abc_nasser_any_same_gene").empirical_p_upper),
        },
    ]
    pd.DataFrame(main_rows).to_csv(MAIN_RESULTS_OUT, sep="\t", index=False)

    observed = denominators[denominators["denominator"].eq("observed_exact_non_eqtl")].iloc[0]
    tiered = pd.DataFrame(
        [
            {
                "tier": "score_atlas",
                "definition": "all 1,587-GWAS score-atlas variants",
                "variants": 1_730_224,
                "broad_loci_500kb": "",
                "moderate_loci": "",
                "systemic_non_neural_loci": "",
                "same_gene_l2g_loci": "",
                "same_gene_e2g_loci": "",
                "same_gene_colocalisation_h4_ge_0_8_loci": "",
                "same_gene_colocalisation_clpp_ge_0_01_loci": "",
            },
            {
                "tier": "recurrent_atlas",
                "definition": "n_traits >= 3",
                "variants": 312_794,
                "broad_loci_500kb": "",
                "moderate_loci": "",
                "systemic_non_neural_loci": "",
                "same_gene_l2g_loci": "",
                "same_gene_e2g_loci": "",
                "same_gene_colocalisation_h4_ge_0_8_loci": "",
                "same_gene_colocalisation_clpp_ge_0_01_loci": "",
            },
            {
                "tier": "core_recurrent_atlas",
                "definition": "n_traits >= 10",
                "variants": 21_027,
                "broad_loci_500kb": "",
                "moderate_loci": "",
                "systemic_non_neural_loci": "",
                "same_gene_l2g_loci": "",
                "same_gene_e2g_loci": "",
                "same_gene_colocalisation_h4_ge_0_8_loci": "",
                "same_gene_colocalisation_clpp_ge_0_01_loci": "",
            },
            {
                "tier": "moderate_family_broad_layer",
                "definition": "moderate-score family-broad regulatory recurrence layer",
                "variants": 7_560,
                "broad_loci_500kb": 2_685,
                "moderate_loci": 2_685,
                "systemic_non_neural_loci": "",
                "same_gene_l2g_loci": "",
                "same_gene_e2g_loci": "",
                "same_gene_colocalisation_h4_ge_0_8_loci": "",
                "same_gene_colocalisation_clpp_ge_0_01_loci": "",
            },
            {
                "tier": "evidence_disjoint_exact",
                "definition": "GTEx QTL gene overlaps non-eQTL SCREEN-linked gene",
                "variants": int(observed.variants),
                "broad_loci_500kb": int(observed.loci_500kb),
                "moderate_loci": int(coloc_row("moderate_exact_non_eqtl").loci_500kb),
                "systemic_non_neural_loci": int(coloc_row("systemic_non_neural_exact").loci_500kb),
                "same_gene_l2g_loci": "",
                "same_gene_e2g_loci": "",
                "same_gene_colocalisation_h4_ge_0_8_loci": "",
                "same_gene_colocalisation_clpp_ge_0_01_loci": "",
            },
            {
                "tier": "external_target_gene_coherence",
                "definition": "evidence-disjoint exact candidates with same-gene Open Targets L2G/e2G agreement at GWAS credible-set coordinates",
                "variants": int(ot_row("observed_exact_non_eqtl").variants),
                "broad_loci_500kb": int(ot_row("observed_exact_non_eqtl").loci_500kb),
                "moderate_loci": int(ot_row("moderate_exact_non_eqtl").loci_500kb),
                "systemic_non_neural_loci": int(ot_row("systemic_non_neural_exact").loci_500kb),
                "same_gene_l2g_loci": int(ot_row("observed_exact_non_eqtl").ot_l2g_gene_agree_score_ge_0_05_loci),
                "same_gene_e2g_loci": int(ot_row("observed_exact_non_eqtl").ot_l2g_e2g_gene_agree_loci),
                "same_gene_colocalisation_h4_ge_0_8_loci": "",
                "same_gene_colocalisation_clpp_ge_0_01_loci": "",
            },
            {
                "tier": "external_colocalisation_support",
                "definition": "evidence-disjoint exact candidates with same-gene Open Targets GWAS-molecular-QTL colocalisation support",
                "variants": int(coloc_row("observed_exact_non_eqtl").variants),
                "broad_loci_500kb": int(coloc_row("observed_exact_non_eqtl").loci_500kb),
                "moderate_loci": int(coloc_row("moderate_exact_non_eqtl").loci_500kb),
                "systemic_non_neural_loci": int(coloc_row("systemic_non_neural_exact").loci_500kb),
                "same_gene_l2g_loci": "",
                "same_gene_e2g_loci": "",
                "same_gene_colocalisation_h4_ge_0_8_loci": int(coloc_row("observed_exact_non_eqtl").ot_coloc_h4_ge_0_8_gene_agree_loci),
                "same_gene_colocalisation_clpp_ge_0_01_loci": int(coloc_row("observed_exact_non_eqtl").ot_coloc_clpp_ge_0_01_gene_agree_loci),
            },
            {
                "tier": "orthogonal_public_database_support",
                "definition": "evidence-disjoint exact candidates with same-gene eQTL Catalogue fine-mapped QTL and ABC enhancer-gene map support",
                "variants": int(eqtl_row("observed_exact_non_eqtl").variants),
                "broad_loci_500kb": int(eqtl_row("observed_exact_non_eqtl").loci_500kb),
                "moderate_loci": int(eqtl_row("moderate_exact_non_eqtl").loci_500kb),
                "systemic_non_neural_loci": int(eqtl_row("systemic_non_neural_exact").loci_500kb),
                "same_gene_l2g_loci": "",
                "same_gene_e2g_loci": "",
                "same_gene_colocalisation_h4_ge_0_8_loci": "",
                "same_gene_colocalisation_clpp_ge_0_01_loci": "",
            },
        ]
    )
    for column in [
        "standalone_e2g_interval_loci",
        "same_gene_pqtl_h4_ge_0_8_loci",
        "eqtl_catalogue_pip_ge_0_1_loci",
        "eqtl_catalogue_pip_ge_0_5_loci",
        "abc_nasser_same_gene_loci",
        "abc_nasser_score_ge_0_05_loci",
        "not_exact_or_ld_proxy_r2_ge_0_8_loci",
        "no_external_gwas_credible_set_broad_locus_overlap_loci",
    ]:
        if column not in tiered.columns:
            tiered[column] = ""
    tiered.loc[tiered["tier"].eq("external_colocalisation_support"), "standalone_e2g_interval_loci"] = int(
        extended_row("observed_exact_non_eqtl").standalone_e2g_gene_agree_loci
    )
    tiered.loc[tiered["tier"].eq("external_colocalisation_support"), "same_gene_pqtl_h4_ge_0_8_loci"] = int(
        extended_row("observed_exact_non_eqtl").same_gene_pqtl_h4_ge_0_8_loci
    )
    tiered.loc[tiered["tier"].eq("orthogonal_public_database_support"), "eqtl_catalogue_pip_ge_0_1_loci"] = int(
        eqtl_row("observed_exact_non_eqtl").eqtl_catalogue_pip_ge_0_1_same_gene_loci
    )
    tiered.loc[tiered["tier"].eq("orthogonal_public_database_support"), "eqtl_catalogue_pip_ge_0_5_loci"] = int(
        eqtl_row("observed_exact_non_eqtl").eqtl_catalogue_pip_ge_0_5_same_gene_loci
    )
    tiered.loc[tiered["tier"].eq("orthogonal_public_database_support"), "abc_nasser_same_gene_loci"] = int(
        abc_row("observed_exact_non_eqtl").abc_nasser_any_same_gene_loci
    )
    tiered.loc[tiered["tier"].eq("orthogonal_public_database_support"), "abc_nasser_score_ge_0_05_loci"] = int(
        abc_row("observed_exact_non_eqtl").abc_nasser_score_ge_0_05_same_gene_loci
    )
    tiered.loc[tiered["tier"].eq("orthogonal_public_database_support"), "not_exact_or_ld_proxy_r2_ge_0_8_loci"] = int(
        novelty_row("observed_exact_non_eqtl").not_exact_or_ld_proxy_r2_ge_0_8_loci
    )
    tiered.loc[tiered["tier"].eq("orthogonal_public_database_support"), "no_external_gwas_credible_set_broad_locus_overlap_loci"] = int(
        novelty_row("observed_exact_non_eqtl").no_external_gwas_credible_set_broad_locus_overlap_loci
    )
    tiered.to_csv(TIERED_STACK_OUT, sep="\t", index=False)

    ot_variant = pd.read_csv(OT_VARIANT_METRICS_OUT, sep="\t")
    variant_evidence = ot_variant.merge(coloc_metrics.drop(columns=["locus_500kb"]), on="variant_key", how="left")
    variant_evidence = variant_evidence.merge(
        extended_metrics.drop(columns=["SNP", "locus_500kb"], errors="ignore"),
        on="variant_key",
        how="left",
    )
    variant_evidence = variant_evidence.merge(
        eqtl_metrics.drop(columns=["locus_500kb"], errors="ignore"),
        on="variant_key",
        how="left",
    )
    variant_evidence = variant_evidence.merge(
        abc_metrics.drop(columns=["locus_500kb"], errors="ignore"),
        on="variant_key",
        how="left",
    )
    for column in [
        "standalone_e2g_gene_agree",
        "same_gene_pqtl_h4_ge_0_8",
        "same_gene_pqtl_clpp_ge_0_01",
        "eqtl_catalogue_any_same_gene",
        "eqtl_catalogue_pip_ge_0_1_same_gene",
        "eqtl_catalogue_pip_ge_0_5_same_gene",
        "abc_nasser_any_same_gene",
        "abc_nasser_score_ge_0_05_same_gene",
        "abc_nasser_score_ge_0_10_same_gene",
    ]:
        if column in variant_evidence.columns:
            variant_evidence[column] = (
                variant_evidence[column].astype(object).where(variant_evidence[column].notna(), False).astype(bool)
            )
    variant_evidence.to_csv(VARIANT_EVIDENCE_OUT, sep="\t", index=False)

    raw = pd.read_csv(
        INPUT,
        sep="\t",
        usecols=[
            "SNP",
            "position",
            "max_score",
            "overlap_non_eqtl_screen_genes",
            "overlap_non_eqtl_screen_count",
        ],
        low_memory=False,
    )
    raw = normalize_coordinates(raw)
    raw["variant_key"] = raw["CHR"].astype(str) + ":" + raw["BP"].astype(str)
    raw["candidate_gene_set"] = raw["overlap_non_eqtl_screen_genes"].map(gene_set)
    evidence = variant_evidence.merge(
        raw[
            [
                "SNP",
                "variant_key",
                "position",
                "max_score",
                "overlap_non_eqtl_screen_genes",
                "overlap_non_eqtl_screen_count",
                "candidate_gene_set",
            ]
        ],
        on=["SNP", "variant_key"],
        how="left",
    )
    coloc_hits = pd.read_csv(OT_COLOC_HITS_OUT, sep="\t")
    coloc_h4 = coloc_hits[coloc_hits["h4_ge_0_8"].astype(bool)]
    coloc_clpp = coloc_hits[coloc_hits["clpp_ge_0_01"].astype(bool)]

    h4_by_locus_gene = coloc_h4.groupby(["locus_500kb", "gene_symbol"]).agg(
        max_coloc_h4=("h4", "max"),
        max_coloc_clpp=("clpp", "max"),
        n_h4_coloc_rows=("rightStudyLocusId", "nunique"),
        n_h4_coloc_types=("rightStudyType", "nunique"),
    )
    clpp_pairs = set(zip(coloc_clpp["locus_500kb"].astype(str), coloc_clpp["gene_symbol"].astype(str)))
    h4_pairs = set(h4_by_locus_gene.index)

    candidate_rows: list[dict[str, object]] = []
    for row in evidence.itertuples(index=False):
        candidate_genes = set(row.candidate_gene_set) if isinstance(row.candidate_gene_set, (set, frozenset)) else set()
        l2g_genes = gene_set(row.ot_l2g_agree_symbols)
        locus_id = str(row.locus_500kb)
        coloc_genes = {gene for locus, gene in h4_pairs if locus == locus_id}
        full_stack_genes = sorted(candidate_genes & l2g_genes & coloc_genes)
        e2g_interval_genes = gene_set(getattr(row, "standalone_e2g_gene_symbols", ""))
        pqtl_h4_genes = gene_set(getattr(row, "same_gene_pqtl_h4_genes", ""))
        eqtl01_genes = gene_set(getattr(row, "eqtl_catalogue_pip_ge_0_1_genes", ""))
        eqtl05_genes = gene_set(getattr(row, "eqtl_catalogue_pip_ge_0_5_genes", ""))
        abc_any_genes = gene_set(getattr(row, "abc_nasser_any_genes", ""))
        abc05_genes = gene_set(getattr(row, "abc_nasser_score_ge_0_05_genes", ""))
        for gene in full_stack_genes:
            h4_stats = h4_by_locus_gene.loc[(locus_id, gene)]
            has_clpp = (locus_id, gene) in clpp_pairs
            has_standalone_e2g = gene in e2g_interval_genes
            has_pqtl_h4 = gene in pqtl_h4_genes
            has_eqtl01 = gene in eqtl01_genes
            has_eqtl05 = gene in eqtl05_genes
            has_abc_any = gene in abc_any_genes
            has_abc05 = gene in abc05_genes
            evidence_score = (
                100
                + 25 * bool(row.ot_l2g_e2g_gene_agree)
                + 20 * bool(row.ot_l2g_top_gene_agree)
                + 20 * bool(row.ot_l2g_coloc_gene_agree)
                + 15 * has_clpp
                + 12 * has_standalone_e2g
                + 15 * has_eqtl05
                + 10 * has_eqtl01
                + 10 * has_abc05
                + 6 * has_abc_any
                + 10 * has_pqtl_h4
                + min(float(row.n_families), 20.0)
                + min(float(row.n_traits) / 10.0, 20.0)
                + min(float(row.ot_l2g_max_score_for_candidate_gene) * 10.0, 10.0)
                + min(float(h4_stats.max_coloc_h4) * 10.0, 10.0)
            )
            candidate_rows.append(
                {
                    "locus_500kb": locus_id,
                    "gene_symbol": gene,
                    "SNP": row.SNP,
                    "position": row.position,
                    "source_set": row.source_set,
                    "score_band": row.score_band,
                    "primary_axis": row.primary_axis,
                    "n_traits": int(row.n_traits),
                    "n_families": int(row.n_families),
                    "max_score": float(row.max_score),
                    "candidate_gene_count": int(row.overlap_non_eqtl_screen_count) if pd.notna(row.overlap_non_eqtl_screen_count) else len(candidate_genes),
                    "ot_gwas_cs_max_pip": float(row.ot_gwas_cs_max_pip),
                    "ot_l2g_max_score_for_candidate_gene": float(row.ot_l2g_max_score_for_candidate_gene),
                    "same_gene_l2g": bool(row.ot_l2g_gene_agree_score_ge_0_05),
                    "same_gene_top_l2g": bool(row.ot_l2g_top_gene_agree),
                    "same_gene_e2g": bool(row.ot_l2g_e2g_gene_agree),
                    "same_gene_l2g_coloc_feature": bool(row.ot_l2g_coloc_gene_agree),
                    "same_gene_ot_coloc_h4_ge_0_8": True,
                    "same_gene_ot_coloc_clpp_ge_0_01": bool(has_clpp),
                    "same_gene_standalone_e2g_interval": bool(has_standalone_e2g),
                    "standalone_e2g_max_score": float(getattr(row, "standalone_e2g_max_score", 0.0) or 0.0),
                    "same_gene_eqtl_catalogue_pip_ge_0_1": bool(has_eqtl01),
                    "same_gene_eqtl_catalogue_pip_ge_0_5": bool(has_eqtl05),
                    "same_gene_abc_nasser": bool(has_abc_any),
                    "same_gene_abc_nasser_score_ge_0_05": bool(has_abc05),
                    "abc_nasser_max_score": float(getattr(row, "abc_nasser_max_score", 0.0) or 0.0),
                    "abc_nasser_biosample_count": int(getattr(row, "abc_nasser_biosample_count", 0) or 0),
                    "same_gene_pqtl_h4_ge_0_8": bool(has_pqtl_h4),
                    "max_coloc_h4": float(h4_stats.max_coloc_h4),
                    "max_coloc_clpp": float(h4_stats.max_coloc_clpp),
                    "n_h4_coloc_study_loci": int(h4_stats.n_h4_coloc_rows),
                    "n_h4_coloc_qtl_types": int(h4_stats.n_h4_coloc_types),
                    "priority_score": float(evidence_score),
                    "candidate_class": "low_ambiguity_candidate" if int(row.overlap_non_eqtl_screen_count) <= 2 else "multi_gene_candidate",
                }
            )

    prioritized_candidates = pd.DataFrame(candidate_rows)
    if not prioritized_candidates.empty:
        prioritized_candidates = prioritized_candidates.sort_values(
            [
                "candidate_class",
                "priority_score",
                "n_h4_coloc_study_loci",
                "n_families",
                "n_traits",
            ],
            ascending=[True, False, False, False, False],
        )
    prioritized_candidates.to_csv(PRIORITIZED_CANDIDATES_OUT, sep="\t", index=False)

    moderate_coloc = null_row(coloc_gene_null, "moderate_exact_non_eqtl", "ot_coloc_h4_ge_0_8_gene_agree")
    systemic_coloc = null_row(coloc_gene_null, "systemic_non_neural_exact", "ot_coloc_h4_ge_0_8_gene_agree")
    moderate_eqtl01 = null_row(eqtl_gene_null, "moderate_exact_non_eqtl", "eqtl_catalogue_pip_ge_0_1_same_gene")
    moderate_eqtl05 = null_row(eqtl_gene_null, "moderate_exact_non_eqtl", "eqtl_catalogue_pip_ge_0_5_same_gene")
    moderate_abc = null_row(abc_gene_null, "moderate_exact_non_eqtl", "abc_nasser_any_same_gene")
    moderate_novelty = novelty_row("moderate_exact_non_eqtl")
    report = f"""# Atlas Exploration Report

## Main Claim
Atlas-scale G-Atlas recurrence reveals a cross-trait regulatory convergence layer whose candidate genes repeatedly align across evidence-disjoint QTL/chromatin links, GWAS-molecular-QTL colocalisation, fine-mapped molecular-QTL, and enhancer-gene evidence.

## Main Results
| analysis | set | observed_loci | null_loci | fold | empirical_p |
|---|---:|---:|---:|---:|---:|
| evidence-disjoint locus convergence | cross-trait layer | {int(locus_row("score_moderate_0.2_0.5").observed_concordant_loci_500kb)} | {locus_row("score_moderate_0.2_0.5").null_locus_mean:.1f} | {locus_row("score_moderate_0.2_0.5").locus_fold_enrichment:.2f} | {locus_row("score_moderate_0.2_0.5").locus_empirical_p_upper:.4g} |
| Open Targets L2G agreement | cross-trait layer | {int(null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_gene_agree_score_ge_0_05").observed_loci)} | {null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_gene_agree_score_ge_0_05").null_locus_mean:.1f} | {null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_gene_agree_score_ge_0_05").fold_enrichment:.2f} | {null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_gene_agree_score_ge_0_05").empirical_p_upper:.4g} |
| Open Targets e2G agreement | cross-trait layer | {int(null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_e2g_gene_agree").observed_loci)} | {null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_e2g_gene_agree").null_locus_mean:.1f} | {null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_e2g_gene_agree").fold_enrichment:.2f} | {null_row(ot_gene_null, "moderate_exact_non_eqtl", "ot_l2g_e2g_gene_agree").empirical_p_upper:.4g} |
| Open Targets GWAS-molQTL colocalisation | cross-trait layer | {int(moderate_coloc.observed_loci)} | {moderate_coloc.null_locus_mean:.1f} | {moderate_coloc.fold_enrichment:.2f} | {moderate_coloc.empirical_p_upper:.4g} |
| Open Targets GWAS-molQTL colocalisation | systemic non-neural | {int(systemic_coloc.observed_loci)} | {systemic_coloc.null_locus_mean:.1f} | {systemic_coloc.fold_enrichment:.2f} | {systemic_coloc.empirical_p_upper:.4g} |
| eQTL Catalogue SuSiE support, PIP >= 0.1 | cross-trait layer | {int(moderate_eqtl01.observed_loci)} | {moderate_eqtl01.null_locus_mean:.1f} | {moderate_eqtl01.fold_enrichment:.2f} | {moderate_eqtl01.empirical_p_upper:.4g} |
| eQTL Catalogue SuSiE support, PIP >= 0.5 | cross-trait layer | {int(moderate_eqtl05.observed_loci)} | {moderate_eqtl05.null_locus_mean:.1f} | {moderate_eqtl05.fold_enrichment:.2f} | {moderate_eqtl05.empirical_p_upper:.4g} |
| ABC enhancer-gene support | cross-trait layer | {int(moderate_abc.observed_loci)} | {moderate_abc.null_locus_mean:.1f} | {moderate_abc.fold_enrichment:.2f} | {moderate_abc.empirical_p_upper:.4g} |

## Supplementary Context
| analysis | set | count |
|---|---:|---:|
| Open Targets standalone E2G interval support | moderate loci | {int(extended_row("moderate_exact_non_eqtl").standalone_e2g_gene_agree_loci)} |
| same-gene pQTL colocalisation, h4 >= 0.8 | moderate loci | {int(extended_row("moderate_exact_non_eqtl").same_gene_pqtl_h4_ge_0_8_loci)} |
| not exact coordinate or r2 >= 0.8 LD proxy of external GWAS credible-set variant | moderate loci | {int(moderate_novelty.not_exact_or_ld_proxy_r2_ge_0_8_loci)} |
| no broad 500 kb external GWAS credible-set overlap | moderate loci | {int(moderate_novelty.no_external_gwas_credible_set_broad_locus_overlap_loci)} |
| prioritized locus-gene candidate rows | full same-gene stack | {len(prioritized_candidates)} |

## Reporting Scope
| result type | role |
|---|---|
| regulatory target-gene convergence map | primary supported output |
| mechanism-ready candidate loci | prioritized follow-up layer |
| external colocalisation, molecular-QTL, enhancer-gene, and public functional evidence | support layers for candidate ranking |
| locus novelty, therapeutic directionality, and experimental validation | reported in dedicated outputs where directly supported |
"""
    REPORT_MD_OUT.write_text(report)


def main() -> int:
    total_steps = 12
    print("Atlas exploration: downstream score-atlas analysis")
    print(f"Package: {display_path(ROOT)}")
    print(f"Raw atlas input: {display_path(RAW_ATLAS_INPUT)}")
    print(f"Public-data cache: {display_path(CACHE_ROOT)}")
    reset_outputs()
    print(f"Outputs reset: {display_path(OUT)}")
    ensure_public_cache_release()

    section(1, total_steps, "Rebuild candidate table from raw atlas", "INPUT")
    build_candidate_table_from_raw_atlas()

    section(2, total_steps, "Evidence-disjoint locus convergence", "MAIN")
    summary = run_analysis()
    denominators = pd.read_csv(DENOMINATORS_OUT, sep="\t")

    print("Denominators")
    for row in denominators.itertuples(index=False):
        print(f"- {row.denominator}: {row.variants:,} variants; {row.loci_500kb:,} 500kb loci")

    print("Primary locus-level results")
    for key in [
        "overall_non_eqtl_tested",
        "score_moderate_0.2_0.5",
        "source_lower_recurrence_exact",
        "axis_systemic_non_neural",
    ]:
        row = summary[summary["set"].eq(key)].iloc[0]
        print(
            f"- {key}: {int(row.observed_concordant_loci_500kb):,} observed loci vs "
            f"{row.null_locus_mean:.1f} null; {row.locus_fold_enrichment:.2f}x; "
            f"p={row.locus_empirical_p_upper:.4g}"
        )

    section(3, total_steps, "Open Targets L2G/e2G target-gene coherence", "MAIN")
    ot_summary, _, ot_gene_null = run_opentargets_anchor()
    exact_row = ot_summary[ot_summary["set"].eq("observed_exact_non_eqtl")].iloc[0]
    moderate_row = ot_summary[ot_summary["set"].eq("moderate_exact_non_eqtl")].iloc[0]
    print(
        "- observed_exact_non_eqtl: "
        f"{int(exact_row.ot_gwas_cs_position_any_loci):,} GWAS credible-set coordinate loci; "
        f"{int(exact_row.ot_l2g_gene_agree_score_ge_0_05_loci):,} L2G-agreeing loci; "
        f"{int(exact_row.ot_l2g_e2g_gene_agree_loci):,} e2G-agreeing loci"
    )
    print(
        "- moderate_exact_non_eqtl: "
        f"{int(moderate_row.ot_gwas_cs_position_any_loci):,} GWAS credible-set coordinate loci; "
        f"{int(moderate_row.ot_l2g_gene_agree_score_ge_0_05_loci):,} L2G-agreeing loci; "
        f"{int(moderate_row.ot_l2g_e2g_gene_agree_loci):,} e2G-agreeing loci"
    )
    key = ot_gene_null[
        ot_gene_null["set"].eq("moderate_exact_non_eqtl")
        & ot_gene_null["indicator"].eq("ot_l2g_e2g_gene_agree")
    ].iloc[0]
    print(
        "- moderate_exact_non_eqtl/e2G gene-label null: "
        f"{int(key.observed_loci):,} observed loci vs {key.null_locus_mean:.1f} null; "
        f"{key.fold_enrichment:.2f}x; p={key.empirical_p_upper:.4g}"
    )

    section(4, total_steps, "Open Targets same-gene GWAS-molQTL colocalisation", "MAIN")
    coloc_summary, coloc_gene_null, coloc_metrics = run_opentargets_colocalisation()
    coloc_key = coloc_gene_null[
        coloc_gene_null["set"].eq("moderate_exact_non_eqtl")
        & coloc_gene_null["indicator"].eq("ot_coloc_h4_ge_0_8_gene_agree")
    ].iloc[0]
    print(
        "- moderate_exact_non_eqtl/h4>=0.8 same-gene colocalisation: "
        f"{int(coloc_key.observed_loci):,} observed loci vs {coloc_key.null_locus_mean:.1f} null; "
        f"{coloc_key.fold_enrichment:.2f}x; p={coloc_key.empirical_p_upper:.4g}"
    )

    section(5, total_steps, "Open Targets E2G and pQTL context", "SUPPLEMENT")
    extended_summary, extended_metrics, _ = run_opentargets_extended_layers()
    ext_row = extended_summary[extended_summary["set"].eq("moderate_exact_non_eqtl")].iloc[0]
    print(
        "- moderate_exact_non_eqtl supplementary context: "
        f"{int(ext_row.standalone_e2g_gene_agree_loci):,} standalone E2G loci; "
        f"{int(ext_row.same_gene_pqtl_h4_ge_0_8_loci):,} pQTL h4>=0.8 loci"
    )

    section(6, total_steps, "eQTL Catalogue SuSiE molecular-QTL support", "MAIN")
    eqtl_summary, eqtl_gene_null, eqtl_metrics = run_eqtl_catalogue_anchor()
    eqtl_key = eqtl_gene_null[
        eqtl_gene_null["set"].eq("moderate_exact_non_eqtl")
        & eqtl_gene_null["indicator"].eq("eqtl_catalogue_pip_ge_0_1_same_gene")
    ].iloc[0]
    print(
        "- moderate_exact_non_eqtl/eQTL Catalogue PIP>=0.1 same-gene: "
        f"{int(eqtl_key.observed_loci):,} observed loci vs {eqtl_key.null_locus_mean:.1f} null; "
        f"{eqtl_key.fold_enrichment:.2f}x; p={eqtl_key.empirical_p_upper:.4g}"
    )

    section(7, total_steps, "ABC Nasser 2021 enhancer-gene support", "MAIN")
    abc_summary, abc_gene_null, abc_metrics = run_abc_nasser_anchor()
    abc_key = abc_gene_null[
        abc_gene_null["set"].eq("moderate_exact_non_eqtl")
        & abc_gene_null["indicator"].eq("abc_nasser_any_same_gene")
    ].iloc[0]
    print(
        "- moderate_exact_non_eqtl/ABC same-gene: "
        f"{int(abc_key.observed_loci):,} observed loci vs {abc_key.null_locus_mean:.1f} null; "
        f"{abc_key.fold_enrichment:.2f}x; p={abc_key.empirical_p_upper:.4g}"
    )

    section(8, total_steps, "Write main and supplementary outputs", "OUTPUT")
    write_clean_main_outputs(
        summary,
        denominators,
        ot_summary,
        ot_gene_null,
        coloc_summary,
        coloc_gene_null,
        coloc_metrics,
        extended_summary,
        extended_metrics,
        eqtl_summary,
        eqtl_gene_null,
        eqtl_metrics,
        abc_summary,
        abc_gene_null,
        abc_metrics,
    )
    print(f"Report: {display_path(REPORT_MD_OUT)}")
    print(f"Main results: {display_path(MAIN_RESULTS_OUT)}")
    print(f"Tiered stack: {display_path(TIERED_STACK_OUT)}")
    print(f"Variant evidence: {display_path(VARIANT_EVIDENCE_OUT)}")
    print(f"Prioritized candidates: {display_path(PRIORITIZED_CANDIDATES_OUT)}")

    section(9, total_steps, "Mechanism-readiness candidate layer", "MIDDLE")
    _, _, _, mechanism_ready, experiment_ready = run_mechanism_readiness_layer()
    print(
        "- mechanism-readiness layer: "
        f"{len(mechanism_ready):,} clean mechanism-ready loci; "
        f"{len(experiment_ready):,} experiment-ready candidate variants"
    )

    section(10, total_steps, "Public functional-support and candidate prioritization layer", "MIDDLE")
    public_support, public_summary, prioritized_mechanisms = run_public_functional_support_layer()
    public_key = public_summary[public_summary["result"].eq("mpravardb_significant_exact_rsid_variant")].iloc[0]
    print(
        "- public functional support: "
        f"{int(public_key.candidate_rows):,} FDR<0.05 MPRA candidate rows; "
        f"{int(public_key.broad_loci_500kb):,} broad loci"
    )
    print(
        "- prioritized candidate mechanisms: "
        f"{len(prioritized_mechanisms):,} ranked variant-gene rows; "
        f"{prioritized_mechanisms['locus_500kb'].nunique():,} broad loci"
    )

    section(11, total_steps, "Conditional recurrence contribution model", "AUDIT")
    _, conditional_yields, conditional_models = run_conditional_recurrence_contribution_model()
    conditional_key = conditional_models[
        conditional_models["model"].eq("external_same_gene_colocalisation_h4")
        & conditional_models["term"].eq("moderate_locus")
    ]
    if not conditional_key.empty:
        row = conditional_key.iloc[0]
        print(
            "- external h4 colocalisation/moderate_locus adjusted audit: "
            f"OR={row.odds_ratio:.2f}; 95% CI {row.ci95_low:.2f}-{row.ci95_high:.2f}; "
            f"p={row.wald_p:.4g}"
        )
    yield_key = conditional_yields[conditional_yields["set"].eq("moderate_exact_loci")].iloc[0]
    print(
        "- moderate_exact_loci descriptive yield: "
        f"{int(yield_key.full_stack_loci):,} full-stack loci; "
        f"{int(yield_key.ot_coloc_h4_loci):,} h4>=0.8 colocalised loci"
    )

    section(12, total_steps, "Leave-one-resource-out target-gene prediction", "SUPPLEMENT")
    loo_summary, _ = run_leave_one_resource_out_target_gene_analysis()
    loo_key = loo_summary[
        loo_summary["set"].eq("moderate_exact_non_eqtl")
        & loo_summary["holdout"].eq("ot_coloc_h4_ge_0_8")
        & loo_summary["mode"].eq("top1_vs_uniform_local_gene_null")
    ].iloc[0]
    print(
        "- moderate_exact_non_eqtl/leave-one-out h4>=0.8 colocalisation: "
        f"{int(loo_key.observed_loci):,} observed loci vs {float(loo_key.null_locus_mean):.1f} null; "
        f"{float(loo_key.fold_enrichment):.2f}x; p={float(loo_key.empirical_p_upper):.4g}"
    )

    print("Clean outputs")
    print(f"- {display_path(REPORT_MD_OUT)}")
    print(f"- {display_path(MAIN_RESULTS_OUT)}")
    print(f"- {display_path(TIERED_STACK_OUT)}")
    print(f"- {display_path(VARIANT_EVIDENCE_OUT)}")
    print(f"- {display_path(PRIORITIZED_CANDIDATES_OUT)}")
    print(f"- {display_path(MECH)}")
    print(f"- {display_path(MECH_PRIORITIZED_CANDIDATES_OUT)}")
    print(f"- {display_path(CONDITIONAL_MODEL_RESULTS_OUT)}")
    print(f"- {display_path(LEAVE_ONE_OUT_SUMMARY_OUT)}")
    print(f"- {display_path(SUPP)}")
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
