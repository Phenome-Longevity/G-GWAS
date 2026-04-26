# Atlas Exploration Report

## Main Claim
Atlas-scale G-Atlas recurrence reveals a cross-trait regulatory convergence layer whose candidate genes repeatedly align across evidence-disjoint QTL/chromatin links, GWAS-molecular-QTL colocalisation, fine-mapped molecular-QTL, and enhancer-gene evidence.

## Main Results
| analysis | set | observed_loci | null_loci | fold | empirical_p |
|---|---:|---:|---:|---:|---:|
| evidence-disjoint locus convergence | cross-trait layer | 738 | 219.7 | 3.36 | 0.0004998 |
| Open Targets L2G agreement | cross-trait layer | 288 | 16.3 | 17.67 | 0.0004998 |
| Open Targets e2G agreement | cross-trait layer | 251 | 15.0 | 16.76 | 0.0004998 |
| Open Targets GWAS-molQTL colocalisation | cross-trait layer | 246 | 15.1 | 16.30 | 0.0004998 |
| Open Targets GWAS-molQTL colocalisation | systemic non-neural | 127 | 6.0 | 21.21 | 0.0004998 |
| eQTL Catalogue SuSiE support, PIP >= 0.1 | cross-trait layer | 55 | 3.3 | 16.55 | 0.0004998 |
| eQTL Catalogue SuSiE support, PIP >= 0.5 | cross-trait layer | 18 | 0.9 | 19.31 | 0.0004998 |
| ABC enhancer-gene support | cross-trait layer | 49 | 1.0 | 47.76 | 0.0004998 |

## Supplementary Context
| analysis | set | count |
|---|---:|---:|
| Open Targets standalone E2G interval support | moderate loci | 239 |
| same-gene pQTL colocalisation, h4 >= 0.8 | moderate loci | 18 |
| not exact coordinate or r2 >= 0.8 LD proxy of external GWAS credible-set variant | moderate loci | 482 |
| no broad 500 kb external GWAS credible-set overlap | moderate loci | 1 |
| prioritized locus-gene candidate rows | full same-gene stack | 448 |

## Reporting Scope
| result type | role |
|---|---|
| regulatory target-gene convergence map | primary supported output |
| mechanism-ready candidate loci | prioritized follow-up layer |
| external colocalisation, molecular-QTL, enhancer-gene, and public functional evidence | support layers for candidate ranking |
| locus novelty, therapeutic directionality, and experimental validation | reported in dedicated outputs where directly supported |


## Mechanism-Ready Candidate Layer
| result | count |
|---|---:|
| single-gene candidate variants | 873 |
| clean single-gene loci | 466 |
| mechanism-ready clean single-gene loci | 224 |
| experiment-ready candidate variants | 430 |

These counts define the practical follow-up layer: low-ambiguity candidate regulatory variant-gene hypotheses suitable for targeted colocalisation, MPRA, CRISPRi, base editing, and disease-context follow-up.


## Public Functional-Support Layer

| result | count |
|---|---:|
| candidate variant-gene rows tested | 796 |
| candidate variants with exact-rsID MPRAVarDB testing | 25 |
| FDR<0.05 MPRAVarDB candidate rows | 3 |
| MPRAbase element-overlap candidate rows | 1 |
| Open Targets CRISPR gene-context candidate rows | 22 |
| ranked candidate mechanisms | 796 |

Public functional genomics resources add orthogonal prioritization evidence from MPRA, enhancer-element assays, CRISPR screens, and disease-cell-context annotations. Candidate rankings are generated directly from atlas recurrence and external evidence features.
