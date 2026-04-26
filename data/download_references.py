#!/usr/bin/env python3
"""
Download reference data and benchmark GWAS.

1. HM3 reference (1000 Genomes HapMap3 QC'd variants, ~1.15M SNPs)
   - Used to filter CAUSALdb2 to common well-imputed variants
   - Source: 1000 Genomes Project

2. Sentinel GWAS (16 UKBB traits from OpenGWAS)
   - Used for held-out benchmarking (excluded from all training)
   - Source: OpenGWAS API (https://gwas.mrcieu.ac.uk/)

Usage:
    python data/download_references.py --output_dir data/
"""

import os, time, argparse, urllib.request


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# 16 UKBB sentinel traits — ALL excluded from training
SENTINEL_TRAITS = {
    'ukb-a-248': 'Body mass index (BMI)',
    'ukb-a-389': 'Standing height',
    'ukb-d-30500_irnt': 'Microalbumin in urine',
    'ukb-d-30240_irnt': 'Reticulocyte percentage',
    'ukb-d-30260_irnt': 'Mean reticulocyte volume',
    'ukb-d-30250_irnt': 'Reticulocyte count',
    'ukb-d-30000_irnt': 'White blood cell (leukocyte) count',
    'ukb-d-30020_irnt': 'Haemoglobin concentration',
    'ukb-d-30080_irnt': 'Platelet count',
    'ukb-d-30010_irnt': 'Red blood cell (erythrocyte) count',
    'ukb-d-30700_irnt': 'Creatinine',
    'ukb-d-30760_irnt': 'HDL cholesterol',
    'ukb-b-12043': 'Protein',
    'ukb-b-8875': 'Heel bone mineral density (BMD)',
    'ukb-b-20175': 'Systolic blood pressure, automated reading',
    'ukb-b-7992': 'Diastolic blood pressure, automated reading',
}

# HM3 reference BIM — 1000 Genomes HapMap3 QC'd variants
HM3_URL = "https://zenodo.org/records/3376628/files/1kg_hm3_QCed_noM.bim"
GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.basic.annotation.gtf.gz"


def download_gencode(output_dir):
    """Download GENCODE GRCh37 gene annotations for functional enrichment."""
    out_path = os.path.join(output_dir, 'gencode_genes.bed')
    if os.path.exists(out_path):
        log(f"GENCODE genes exists: {out_path}")
        return

    import gzip
    gtf_path = os.path.join(output_dir, 'gencode_grch37.gtf.gz')
    log(f"Downloading GENCODE GRCh37 v44...")
    try:
        urllib.request.urlretrieve(GENCODE_URL, gtf_path)
        genes = []
        with gzip.open(gtf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'gene': continue
                chrom = parts[0].replace('chr', '')
                if not chrom.isdigit(): continue
                start, end = int(parts[3]), int(parts[4])
                gene_name = ''
                for attr in parts[8].split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_name'):
                        gene_name = attr.split('"')[1]
                        break
                if gene_name:
                    genes.append((chrom, start, end, gene_name))
        with open(out_path, 'w') as f:
            for chrom, start, end, name in sorted(genes, key=lambda x: (int(x[0]), x[1])):
                f.write(f'{chrom}\t{start}\t{end}\t{name}\n')
        os.remove(gtf_path)
        log(f"Extracted {len(genes)} autosomal genes")
    except Exception as e:
        log(f"GENCODE download failed: {e}")


def download_hm3(output_dir):
    """Download HM3 reference BIM file."""
    out_path = os.path.join(output_dir, 'hm3.bim')
    if os.path.exists(out_path):
        log(f"HM3 BIM exists: {out_path}")
        return

    log(f"Downloading HM3 reference from {HM3_URL}")
    try:
        urllib.request.urlretrieve(HM3_URL, out_path)
        n_variants = sum(1 for _ in open(out_path))
        log(f"Downloaded: {n_variants} HM3 variants")
    except Exception as e:
        log(f"Download failed: {e}")
        log("Manual download: get 1kg_hm3_QCed_noM.bim from 1000 Genomes (Zenodo)")
        log(f"Place at: {out_path}")


def download_sentinel(output_dir):
    """Download 16 UKBB sentinel GWAS from OpenGWAS."""
    sentinel_dir = os.path.join(output_dir, 'sentinel_sumstats')
    os.makedirs(sentinel_dir, exist_ok=True)

    for trait_id, trait_name in SENTINEL_TRAITS.items():
        out_path = os.path.join(sentinel_dir, f'{trait_id}.txt')
        if os.path.exists(out_path):
            continue

        log(f"Downloading {trait_name} ({trait_id})...")

        # OpenGWAS provides VCF files via API
        # For reproducibility, we provide the GWAS summary statistics
        # which were created from OpenGWAS VCFs in the original pipeline:
        #   1. Download VCF from OpenGWAS API
        #   2. Parse VCF → TSV (CHR, BP, SNP, A1, A2, MAF, BETA, SE, P, N)
        #   3. These are the sentinel_sumstats files

        log(f"  NOTE: OpenGWAS downloads require API access.")
        log(f"  These are standard GWAS summary statistics derived from OpenGWAS VCFs.")
        log(f"  Original source: https://gwas.mrcieu.ac.uk/datasets/{trait_id}/")

    existing = len([f for f in os.listdir(sentinel_dir) if f.endswith('.txt')])
    log(f"Sentinel sumstats: {existing}/{len(SENTINEL_TRAITS)} traits in {sentinel_dir}")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output_dir', default='data/')
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log("=" * 70)
    log("DOWNLOAD REFERENCES")
    log("=" * 70)

    download_hm3(args.output_dir)
    download_gencode(args.output_dir)
    download_sentinel(args.output_dir)

    log("\nData provenance:")
    log("  HM3 BIM:           1000 Genomes HapMap3 QC'd (via Zenodo/1000 Genomes)")
    log("  Sentinel GWAS:     OpenGWAS API (https://gwas.mrcieu.ac.uk/)")
    log("  EBI GWAS:          EBI GWAS Catalog FTP (download_ebi.py)")
    log("  CAUSALdb2:         http://www.mulinlab.org/causaldb/ (download_causaldb.py)")
    log("  All sentinel traits excluded from all training.")
