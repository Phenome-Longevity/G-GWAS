#!/usr/bin/env python3
"""
Step 0b: Preprocess raw GWAS TSVs → binary .npz blocks for encoder pretraining.

For each GWAS:
  1. Read harmonised TSV (chromosome, base_pair_location, beta, SE, p, MAF)
  2. Split into 500kb genomic blocks
  3. Per block: build 6-feature vector per variant
  4. Save as .npz (fast reload for pretraining)

Features per variant (6):
  0: beta (effect size)
  1: standard_error
  2: -log10(p_value)
  3: effect_allele_frequency (MAF)
  4: relative position within 500kb block (0-1)
  5: chromosome / 23

Usage:
    python data/preprocess.py --input_dir /path/to/raw_gwas/ --output_dir /path/to/npz/ --workers 7
"""

import sys, os, time, argparse
import numpy as np
import pandas as pd
from pathlib import Path
from multiprocessing import Pool

BLOCK_SIZE = 500_000       # 500kb genomic blocks
MAX_VARIANTS_PER_BLOCK = 500
MAX_BLOCKS_PER_GWAS = 200  # subsample if more blocks


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def preprocess_one(args):
    gwas_path, out_path = args
    if os.path.exists(out_path):
        return ('skip', gwas_path, 0)

    try:
        df = pd.read_csv(gwas_path, sep='\t', comment='#', low_memory=False,
                         usecols=lambda c: c in [
                             'chromosome', 'base_pair_location', 'beta',
                             'standard_error', 'p_value', 'effect_allele_frequency',
                             'odds_ratio'])

        if 'chromosome' not in df.columns or 'base_pair_location' not in df.columns:
            return ('no_cols', gwas_path, 0)

        df = df.dropna(subset=['chromosome', 'base_pair_location'])
        df['chromosome'] = df['chromosome'].astype(str).str.replace('chr', '')
        df = df[df['chromosome'].isin([str(i) for i in range(1, 23)])]
        if len(df) < 100:
            return ('small', gwas_path, len(df))

        df['chromosome'] = df['chromosome'].astype(int)
        df['base_pair_location'] = df['base_pair_location'].astype(int)
        df['block'] = df['base_pair_location'] // BLOCK_SIZE

        # Beta
        if 'beta' in df.columns:
            df['beta'] = pd.to_numeric(df['beta'], errors='coerce')
        if 'beta' not in df.columns or df['beta'].isna().all():
            if 'odds_ratio' in df.columns:
                df['beta'] = np.log(pd.to_numeric(df['odds_ratio'], errors='coerce'))
            if 'beta' not in df.columns or df['beta'].isna().all():
                return ('no_beta', gwas_path, 0)

        df['beta'] = df['beta'].fillna(0).astype(np.float32)
        df['se'] = pd.to_numeric(df.get('standard_error', pd.Series(dtype=float)),
                                  errors='coerce').fillna(1).astype(np.float32)
        pvals = pd.to_numeric(df['p_value'], errors='coerce').clip(lower=1e-300)
        df['neg_log_p'] = -np.log10(pvals).astype(np.float32)
        df['maf'] = pd.to_numeric(df.get('effect_allele_frequency', pd.Series(dtype=float)),
                                   errors='coerce').fillna(0.5).astype(np.float32)

        # Build blocks
        block_keys = sorted(df.groupby(['chromosome', 'block']).groups.keys())
        if len(block_keys) > MAX_BLOCKS_PER_GWAS:
            rng = np.random.default_rng(hash(str(gwas_path)) % (2**32))
            idxs = rng.choice(len(block_keys), MAX_BLOCKS_PER_GWAS, replace=False)
            block_keys = [block_keys[i] for i in sorted(idxs)]

        blocks = []
        total_variants = 0
        for chrom, blk in block_keys:
            sub = df[(df['chromosome'] == chrom) & (df['block'] == blk)]
            if len(sub) < 10:
                continue
            if len(sub) > MAX_VARIANTS_PER_BLOCK:
                sub = sub.sample(MAX_VARIANTS_PER_BLOCK, random_state=42)
            sub = sub.sort_values('base_pair_location')

            n = len(sub)
            feat = np.zeros((n, 6), dtype=np.float32)
            feat[:, 0] = sub['beta'].values
            feat[:, 1] = sub['se'].values
            feat[:, 2] = sub['neg_log_p'].values
            feat[:, 3] = sub['maf'].values
            feat[:, 4] = (sub['base_pair_location'].values - blk * BLOCK_SIZE) / BLOCK_SIZE
            feat[:, 5] = chrom / 23.0

            if np.isfinite(feat).all():
                blocks.append(feat)
                total_variants += n

        if len(blocks) < 5:
            return ('few_blocks', gwas_path, len(blocks))

        save_dict = {f'block_{i}': blocks[i] for i in range(len(blocks))}
        save_dict['n_blocks'] = np.array([len(blocks)])
        save_dict['n_variants'] = np.array([total_variants])
        np.savez_compressed(out_path, **save_dict)
        return ('ok', gwas_path, total_variants)

    except Exception as e:
        return ('error', gwas_path, str(e)[:80])


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--input_dir', required=True, help='Directory of raw GWAS .tsv.gz files')
    p.add_argument('--output_dir', required=True, help='Directory for preprocessed .npz files')
    p.add_argument('--workers', type=int, default=4)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log("=" * 70)
    log("PREPROCESS GWAS → .npz blocks")
    log(f"  Block size: {BLOCK_SIZE} bp")
    log(f"  Max variants/block: {MAX_VARIANTS_PER_BLOCK}")
    log(f"  Max blocks/GWAS: {MAX_BLOCKS_PER_GWAS}")
    log(f"  Features: [beta, SE, -log10p, MAF, rel_pos, chrom/23]")
    log("=" * 70)

    gwas_files = sorted(Path(args.input_dir).glob('*.tsv.gz'))
    tasks = [(str(gf), os.path.join(args.output_dir, gf.stem.split('.')[0] + '.npz'))
             for gf in gwas_files]
    log(f"Found {len(tasks)} GWAS files, using {args.workers} workers")

    stats = {'ok': 0, 'skip': 0, 'error': 0, 'no_beta': 0, 'few_blocks': 0, 'no_cols': 0, 'small': 0}
    t0 = time.time()

    with Pool(args.workers) as pool:
        for i, result in enumerate(pool.imap_unordered(preprocess_one, tasks)):
            status = result[0]
            stats[status] = stats.get(status, 0) + 1
            if (i + 1) % 500 == 0:
                elapsed = time.time() - t0
                rate = (i + 1) / elapsed
                log(f"  {i+1}/{len(tasks)} ({rate:.1f}/s) | ok={stats['ok']} skip={stats['skip']}")

    log(f"\nDone in {(time.time()-t0)/60:.1f} min. Stats: {stats}")
