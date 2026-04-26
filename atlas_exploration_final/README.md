# Atlas Exploration

Self-contained downstream analysis of the `1,587`-GWAS G-Atlas score atlas.

## Entry Point

```bash
python3 -u run_atlas_exploration.py 2>&1 | tee logs/run_atlas_exploration.log
```

## Inputs

| file | role |
|---|---|
| `inputs/no_window_gatlas_1587_gwas.parquet` | fixed G-Atlas score-atlas table: 1,730,224 variants aggregated across 1,587 GWAS |

The input parquet is the handoff from model inference to downstream atlas analysis. Each row is a variant with `CHR`, `BP`, `snp`, `n_traits`, `traits`, `scores`, `max_score`, and `mean_score`. The `traits` and `scores` list columns record the GWAS accessions and G-Atlas scores in which that variant was observed. The file name preserves the original run label; the analysis treats it as the frozen `1,587`-GWAS score atlas.

Public annotation resources are cached inside this package under `cache_v3/`. The run checks the manifest, downloads missing or size-mismatched public cache archives from `https://gatlas.szags.uk/`, extracts them, and then regenerates the outputs.

| cache | resource |
|---|---|
| `cache_v3/gtex_v8` | GTEx v8 significant eQTL and sQTL tar archives |
| `cache_v3/screen_ccre` | SCREEN cCRE BED files and Human-Gene-Links |
| `cache_v3/gwas_catalog` | GWAS Catalog association table |
| `cache_v3/causaldb` | CAUSALdb credible-set resources used for exact-rsID flags |
| `cache_v3/trait_family_map` | GCST-to-trait-family annotation used to compute recurrence breadth |
| `cache_v3/public_functional_support` | public MPRA, CRISPRi-FlowFISH, MPRAbase, MPRAVarDB, and liftOver resources |
| `cache_v3/opentargets_platform_26_03` | Open Targets 26.03 credible sets, L2G/e2G, colocalisation, E2G, pQTL context |
| `cache_v3/eqtl_catalogue_r8_beta_susie` | eQTL Catalogue r8 beta SuSiE molecular-QTL credible sets |
| `cache_v3/abc_nasser2021` | Nasser 2021 ABC enhancer-gene predictions |

## Main Outputs

| file | content |
|---|---|
| `outputs/atlas_report.md` | concise report with main results, supplementary context, and reporting scope |
| `outputs/atlas_main_results.tsv` | main numerical results for Figure 4-style reporting |
| `outputs/atlas_tiered_mechanism_stack.tsv` | tiered atlas-to-target-gene evidence stack |
| `outputs/atlas_variant_evidence.tsv` | per-variant public evidence table used by the final atlas report |
| `outputs/atlas_prioritized_locus_gene_candidates.tsv` | ranked locus-gene candidates with same-gene GTEx/non-eQTL SCREEN, L2G, and colocalisation support |
| `logs/run_atlas_exploration.log` | timestamped section-by-section run log |

## Mechanism-Readiness Outputs

These tables form the practical follow-up layer between the main evidence stack and supplementary audits.

| file | content |
|---|---|
| `outputs/mechanism_readiness/atlas_mechanism_readiness_report.md` | concise mechanism-readiness summary |
| `outputs/mechanism_readiness/mechanism_discovery_taxonomy.tsv` | candidate counts by ambiguity, score band, source, and axis |
| `outputs/mechanism_readiness/target_gene_recurrence.tsv` | recurrent target genes across loci and trait axes |
| `outputs/mechanism_readiness/locus_ambiguity.tsv` | broad-locus target ambiguity and support summary |
| `outputs/mechanism_readiness/mechanism_ready_loci.tsv` | clean single-gene, 3D-supported loci |
| `outputs/mechanism_readiness/experiment_ready_candidates.tsv` | ranked variants for targeted validation |
| `outputs/mechanism_readiness/public_functional_support.tsv` | public MPRA, MPRA-element, CRISPRi, and CRISPR gene-context support joined to candidate pairs |
| `outputs/mechanism_readiness/public_functional_support_summary.tsv` | compact counts for public functional-support layers |
| `outputs/mechanism_readiness/prioritized_candidate_mechanisms.tsv` | naturally ranked candidate mechanisms from atlas recurrence plus external evidence |
| `outputs/mechanism_readiness/candidate_mechanism_dossiers.md` | concise ranked candidate summaries generated from the scoring table |

## Main Analyses

| step | analysis |
|---:|---|
| 1 | rebuild candidate regulatory layer from the raw atlas parquet |
| 2 | evidence-disjoint GTEx/non-eQTL SCREEN convergence after 500 kb locus collapse |
| 3 | Open Targets L2G/e2G target-gene agreement at GWAS credible-set coordinates |
| 4 | Open Targets same-gene GWAS-molecular-QTL colocalisation |
| 5 | eQTL Catalogue SuSiE same-variant, same-gene molecular-QTL support |
| 6 | Nasser 2021 ABC same-variant, same-gene enhancer-gene support |
| 7 | mechanism-readiness layer for clean candidate regulatory variant-gene hypotheses |
| 8 | public functional-support and candidate-prioritization layer |

## Supplementary Outputs

Supplementary tables are written under `outputs/supplementary/`.

| group | role |
|---|---|
| `locus_convergence_*` | denominator, null draws, and locus-level convergence summaries |
| `opentargets26_*` | Open Targets coordinate, L2G/e2G, colocalisation, E2G, pQTL, and target-context tables |
| `eqtl_catalogue_*` | eQTL Catalogue same-gene molecular-QTL support and nulls |
| `abc_nasser2021_*` | ABC enhancer-gene support and nulls |
| `conditional_recurrence_*` | adjusted audit of recurrence contribution after broad covariates |
| `leave_one_resource_out_*` | held-out evidence prediction tables and null draws |

## Reporting Scope

| result type | role |
|---|---|
| regulatory target-gene convergence map | primary supported output |
| mechanism-ready candidate loci | prioritized follow-up layer |
| external colocalisation, molecular-QTL, enhancer-gene, and public functional evidence | support layers for candidate ranking |
| locus novelty, therapeutic directionality, and experimental validation | reported in dedicated outputs where directly supported |
