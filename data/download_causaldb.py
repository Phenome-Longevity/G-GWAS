#!/usr/bin/env python3
"""
Step 0c: Download and prepare CAUSALdb2 fine-mapping data.

CAUSALdb2 contains real fine-mapping PIPs from 7 methods (SuSiE, FINEMAP,
PAINTOR, CAVIARBF, PolyFun-FINEMAP, PolyFun-SuSiE, ABF) across ~15K GWAS.

Source: http://www.mulinlab.org/causaldb/
Data: credible_set.txt (36M variants, 6.6 GB)

This script downloads the data and filters to HM3 variants for training.

Usage:
    python data/download_causaldb.py --output_dir /path/to/causaldb/
"""

import os, time, argparse
import pandas as pd
import numpy as np


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# CAUSALdb2 API
CAUSALDB_API = "http://www.mulinlab.org/causaldb/block/readFile"
# Fields: t0=ABF, t1=FINEMAP, t2=PAINTOR, t3=CAVIARBF, t4=SUSIE, t5=POLYFUN_FINEMAP, t6=POLYFUN_SUSIE


def download_causaldb(output_dir):
    """Download the full CAUSALdb2 credible set file."""
    log("NOTE: CAUSALdb2 download requires manual steps:")
    log("  1. Go to http://www.mulinlab.org/causaldb/")
    log("  2. Download credible_set.txt (v2.1, ~6.6 GB)")
    log("  3. Place in output directory")
    log("")
    log("Or use the API for individual GWAS:")
    log(f"  POST {CAUSALDB_API}")
    log("  Params: fileName, type, chr, posStart, posEnd, id, popu")
    log("  Fields: t0=ABF, t1=FINEMAP, t2=PAINTOR, t3=CAVIARBF,")
    log("          t4=SUSIE, t5=POLYFUN_FINEMAP, t6=POLYFUN_SUSIE")


def filter_to_hm3(credible_set_path, hm3_bim_path, output_path):
    """Filter CAUSALdb2 to HM3 variants and save as parquet."""
    log(f"Loading CAUSALdb2: {credible_set_path}")
    # The credible set file is large — read in chunks
    chunks = []
    for chunk in pd.read_csv(credible_set_path, sep='\t', chunksize=1_000_000):
        # Standardize columns
        chunk = chunk.rename(columns=lambda c: c.lower().strip())
        chunks.append(chunk)
        if len(chunks) % 10 == 0:
            log(f"  Read {len(chunks)}M rows")

    df = pd.concat(chunks, ignore_index=True)
    log(f"Total: {len(df)} variants, {df['meta_id'].nunique() if 'meta_id' in df.columns else '?'} GWAS")

    # Load HM3 positions
    if hm3_bim_path and os.path.exists(hm3_bim_path):
        log(f"Filtering to HM3: {hm3_bim_path}")
        hm3 = pd.read_csv(hm3_bim_path, sep='\t', header=None,
                           names=['chr', 'rsid', 'cm', 'bp', 'a1', 'a2'])
        hm3_pos = set(zip(hm3['chr'], hm3['bp']))
        before = len(df)
        df = df[df.apply(lambda r: (int(r['chr']), int(r['bp'])) in hm3_pos, axis=1)]
        log(f"  {before} → {len(df)} variants (HM3 filter)")

    # Save as parquet
    df.to_parquet(output_path, index=False)
    log(f"Saved: {output_path} ({os.path.getsize(output_path)/1e9:.1f} GB)")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--output_dir', required=True)
    p.add_argument('--credible_set', default=None, help='Path to credible_set.txt if already downloaded')
    p.add_argument('--hm3_bim', default=None, help='HM3 .bim file for filtering')
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.credible_set:
        output = os.path.join(args.output_dir, 'susie_pips.parquet')
        filter_to_hm3(args.credible_set, args.hm3_bim, output)
    else:
        download_causaldb(args.output_dir)
