#!/usr/bin/env python3
"""
G-Atlas inference: GWAS summary statistics → causal variant scores.

Usage:
    python predict.py input.tsv output.scores.tsv
    python predict.py input.tsv output.scores.tsv --region chr6:25000000-35000000
    python predict.py input.tsv output.scores.tsv --device cuda

Input: GWAS summary statistics (TSV, optionally gzipped).
Output: TSV with columns: SNP, CHR, BP, SCORE
    SCORE = G-Atlas causal prioritization score (0-1).
    Note: these are neural network predictions trained to approximate
    SuSiE posteriors, not Bayesian posterior inclusion probabilities.

GPU-accelerated: uses CUDA for LD features + model inference when available.
Falls back to CPU automatically.
"""

import sys, os, time, argparse
import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.model import GAtlasEncoder, PIPHead
from src.features import compute_ld_auto, compute_ld_features, BLOCK_SIZE, N_LD_FEATURES

MAX_VARIANTS_PER_BLOCK = 500
MODEL_DIR = os.path.join(os.path.dirname(__file__), '..', 'models')

# Try polars for fast CSV reading, fall back to pandas
try:
    import polars as pl
    HAS_POLARS = True
except ImportError:
    HAS_POLARS = False
import pandas as pd


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_gwas(input_path, region=None):
    """Read GWAS file. Uses polars if available (~50× faster)."""
    t0 = time.time()

    if HAS_POLARS:
        df = pl.read_csv(input_path, separator='\t',
                         infer_schema_length=10000, ignore_errors=True)
        # Standardize column names — handle EBI harmonized (hm_*) and standard formats
        # Priority: hm_ prefixed > standard (hm_ are harmonized, more reliable)
        col_map = {}
        cols_lower = {c: c.lower() for c in df.columns}
        mapped = set()  # track which target columns are already mapped

        # First pass: hm_ prefixed columns (highest priority, only if they have real data)
        def has_data(col_name):
            try:
                s = df[col_name].cast(str)
                bad = s.is_null() | s.str.to_uppercase().is_in(["NA", "NULL", ""])
                return int(bad.sum()) < len(df)
            except Exception:
                return False
        for col, cl in cols_lower.items():
            if cl == 'hm_chrom' and 'CHR' not in mapped and has_data(col): col_map[col] = 'CHR'; mapped.add('CHR')
            elif cl == 'hm_pos' and 'BP' not in mapped and has_data(col): col_map[col] = 'BP'; mapped.add('BP')
            elif cl == 'hm_rsid' and 'SNP' not in mapped and has_data(col): col_map[col] = 'SNP'; mapped.add('SNP')
            elif cl == 'hm_beta' and 'BETA' not in mapped and has_data(col): col_map[col] = 'BETA'; mapped.add('BETA')
            elif cl in ('hm_se', 'hm_standard_error') and 'SE' not in mapped and has_data(col): col_map[col] = 'SE'; mapped.add('SE')
            elif cl in ('hm_p_value', 'hm_pval') and 'P' not in mapped and has_data(col): col_map[col] = 'P'; mapped.add('P')
            elif cl in ('hm_effect_allele_frequency', 'hm_eaf') and 'MAF' not in mapped and has_data(col): col_map[col] = 'MAF'; mapped.add('MAF')

        # Second pass: standard columns (only if not already mapped)
        for col, cl in cols_lower.items():
            if cl in ('chromosome', 'chr', 'chrom', '#chrom') and 'CHR' not in mapped: col_map[col] = 'CHR'; mapped.add('CHR')
            elif cl in ('base_pair_location', 'bp', 'pos', 'position') and 'BP' not in mapped: col_map[col] = 'BP'; mapped.add('BP')
            elif cl in ('rsid', 'snp', 'rsids') and 'SNP' not in mapped: col_map[col] = 'SNP'; mapped.add('SNP')
            elif cl in ('variant_id', 'id') and 'SNP' not in mapped: col_map[col] = 'SNP'; mapped.add('SNP')
            elif cl in ('beta', 'effect') and 'BETA' not in mapped: col_map[col] = 'BETA'; mapped.add('BETA')
            elif cl in ('standard_error', 'se', 'sebeta') and 'SE' not in mapped: col_map[col] = 'SE'; mapped.add('SE')
            elif cl in ('p_value', 'p', 'pval', 'pvalue') and 'P' not in mapped: col_map[col] = 'P'; mapped.add('P')
            elif cl in ('effect_allele_frequency', 'maf', 'af', 'eaf') and 'MAF' not in mapped: col_map[col] = 'MAF'; mapped.add('MAF')

        # Compute P from neg_log10_p_value or log(P) if no P column
        if 'P' not in mapped:
            for col, cl in cols_lower.items():
                if cl in ('neg_log10_p_value', 'log(p)', 'log_p'):
                    df = df.with_columns((10.0 ** (-pl.col(col).cast(pl.Float64).abs())).alias('P'))
                    mapped.add('P')
                    break

        df = df.rename(col_map)
        df = df.drop_nulls(subset=['CHR', 'BP', 'BETA', 'SE', 'P'])
        if 'SNP' not in df.columns:
            df = df.with_columns(
                (pl.lit('chr') + pl.col('CHR').cast(str) + pl.lit(':') + pl.col('BP').cast(str)).alias('SNP')
            )
        df = df.unique('SNP')
        df = df.to_pandas()
    else:
        df = pd.read_csv(input_path, sep='\t', comment='#', low_memory=False)
        col_map = {}
        mapped = set()
        cols_lower = {c: c.lower() for c in df.columns}
        def has_data(col_name):
            s = df[col_name]
            bad = s.isna() | s.astype(str).str.upper().isin(["NA", "NULL", ""])
            return int(bad.sum()) < len(df)
        # hm_ prefixed first
        for col, cl in cols_lower.items():
            if cl == 'hm_chrom' and 'CHR' not in mapped: col_map[col] = 'CHR'; mapped.add('CHR')
            elif cl == 'hm_pos' and 'BP' not in mapped: col_map[col] = 'BP'; mapped.add('BP')
            elif cl == 'hm_rsid' and 'SNP' not in mapped: col_map[col] = 'SNP'; mapped.add('SNP')
            elif cl == 'hm_beta' and 'BETA' not in mapped: col_map[col] = 'BETA'; mapped.add('BETA')
            elif cl in ('hm_se', 'hm_standard_error') and 'SE' not in mapped and has_data(col): col_map[col] = 'SE'; mapped.add('SE')
            elif cl in ('hm_p_value', 'hm_pval') and 'P' not in mapped and has_data(col): col_map[col] = 'P'; mapped.add('P')
            elif cl in ('hm_effect_allele_frequency', 'hm_eaf') and 'MAF' not in mapped and has_data(col): col_map[col] = 'MAF'; mapped.add('MAF')
        # Standard columns
        for col, cl in cols_lower.items():
            if cl in ('chromosome', 'chr', 'chrom', '#chrom') and 'CHR' not in mapped: col_map[col] = 'CHR'; mapped.add('CHR')
            elif cl in ('base_pair_location', 'bp', 'pos', 'position') and 'BP' not in mapped: col_map[col] = 'BP'; mapped.add('BP')
            elif cl in ('rsid', 'snp', 'rsids') and 'SNP' not in mapped: col_map[col] = 'SNP'; mapped.add('SNP')
            elif cl in ('variant_id', 'id') and 'SNP' not in mapped: col_map[col] = 'SNP'; mapped.add('SNP')
            elif cl in ('beta', 'effect') and 'BETA' not in mapped: col_map[col] = 'BETA'; mapped.add('BETA')
            elif cl in ('standard_error', 'se', 'sebeta') and 'SE' not in mapped: col_map[col] = 'SE'; mapped.add('SE')
            elif cl in ('p_value', 'p', 'pval', 'pvalue') and 'P' not in mapped: col_map[col] = 'P'; mapped.add('P')
            elif cl in ('effect_allele_frequency', 'maf', 'af', 'eaf') and 'MAF' not in mapped: col_map[col] = 'MAF'; mapped.add('MAF')
        # Compute P from log if needed
        if 'P' not in mapped:
            for col, cl in cols_lower.items():
                if cl in ('neg_log10_p_value', 'log(p)', 'log_p'):
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                    df['P'] = 10.0 ** (-df[col].abs())
                    mapped.add('P')
                    break
        df = df.rename(columns=col_map)
        df = df.dropna(subset=['CHR', 'BP', 'BETA', 'SE', 'P'])

    # Filter to autosomes, handle CHR formatting
    df['CHR'] = df['CHR'].astype(str).str.replace('chr', '')
    df = df[df['CHR'].isin([str(i) for i in range(1, 23)])]
    df['CHR'] = df['CHR'].astype(int)
    df['BP'] = df['BP'].astype(int)

    # Drop rows with missing critical values instead of filling
    df['BETA'] = pd.to_numeric(df['BETA'], errors='coerce')
    df['SE'] = pd.to_numeric(df['SE'], errors='coerce')
    df = df.dropna(subset=['BETA', 'SE'])
    df['BETA'] = df['BETA'].astype(np.float32)
    df['SE'] = df['SE'].astype(np.float32)
    df['P'] = pd.to_numeric(df['P'], errors='coerce').clip(lower=1e-300)
    if 'MAF' not in df.columns:
        df['MAF'] = 0.5
    df['MAF'] = pd.to_numeric(df['MAF'], errors='coerce').fillna(0.5).astype(np.float32)

    # Create SNP column BEFORE dedup if missing
    if 'SNP' not in df.columns:
        df['SNP'] = 'chr' + df['CHR'].astype(str) + ':' + df['BP'].astype(str)
    df = df.drop_duplicates('SNP')

    # Filter by region if specified
    if region:
        chrom, positions = region.replace('chr', '').split(':')
        start, end = positions.split('-')
        df = df[(df['CHR'] == int(chrom)) & (df['BP'] >= int(start)) & (df['BP'] <= int(end))]

    elapsed = time.time() - t0
    reader = 'polars' if HAS_POLARS else 'pandas'
    log(f"Read ({reader}): {len(df):,} variants in {elapsed:.1f}s")
    return df


def predict_gwas(input_path, output_path, device='cpu', region=None):
    """Run G-Atlas on a GWAS file."""

    encoder = GAtlasEncoder().to(device)
    encoder.load_state_dict(torch.load(
        os.path.join(MODEL_DIR, 'encoder.pt'), weights_only=True, map_location=device))
    if device == 'cuda':
        encoder = encoder.half()
    encoder.eval()

    pip_head = PIPHead().to(device)
    pip_head.load_state_dict(torch.load(
        os.path.join(MODEL_DIR, 'pip_head.pt'), weights_only=True, map_location=device))
    if device == 'cuda':
        pip_head = pip_head.half()
    pip_head.eval()

    df = read_gwas(input_path, region)
    if len(df) == 0:
        log("No variants after filtering")
        return

    # Build features as numpy, sort by (chr, block, bp)
    t0 = time.time()
    chrs = df['CHR'].values.astype(np.int32)
    bps = df['BP'].values.astype(np.int64)
    block_ids = bps // BLOCK_SIZE
    sort_idx = np.lexsort((bps, block_ids, chrs))

    N = len(sort_idx)
    all_feat = np.zeros((N, 6), dtype=np.float32)
    all_feat[:, 0] = df['BETA'].values.astype(np.float32)[sort_idx]
    all_feat[:, 1] = df['SE'].values.astype(np.float32)[sort_idx]
    all_feat[:, 2] = -np.log10(np.clip(df['P'].values.astype(np.float64)[sort_idx], 1e-300, 1.0)).astype(np.float32)
    all_feat[:, 3] = df['MAF'].values.astype(np.float32)[sort_idx]
    all_feat[:, 4] = (bps[sort_idx] - block_ids[sort_idx] * BLOCK_SIZE) / BLOCK_SIZE
    all_feat[:, 5] = chrs[sort_idx] / 23.0
    all_feat = np.nan_to_num(all_feat).astype(np.float32)

    # Find block boundaries
    combined = chrs[sort_idx].astype(np.int64) * 1_000_000 + block_ids[sort_idx]
    changes = np.where(np.diff(combined) != 0)[0] + 1
    bounds = np.concatenate([[0], changes, [N]])

    # Build chunks (500-variant windows)
    chunks = []
    for i in range(len(bounds) - 1):
        s, e = bounds[i], bounds[i+1]
        if e - s < 5: continue
        for wi in range(s, e, MAX_VARIANTS_PER_BLOCK):
            ce = min(wi + MAX_VARIANTS_PER_BLOCK, e)
            if ce - wi < 3: continue
            chunks.append((wi, ce))

    elapsed_prep = time.time() - t0
    n_total = sum(e - s for s, e in chunks)
    log(f"Prep: {len(chunks)} chunks, {n_total:,} variants in {elapsed_prep:.1f}s")

    # GPU inference
    t0 = time.time()
    all_scores = np.full(N, np.nan, dtype=np.float32)
    BATCH = 256 if device == 'cuda' else 32

    use_amp = device == 'cuda'
    ctx = torch.amp.autocast('cuda') if use_amp else torch.no_grad()

    with torch.no_grad(), ctx:
        for bi in range(0, len(chunks), BATCH):
            bc = chunks[bi:bi+BATCH]
            B = len(bc); ml = max(e - s for s, e in bc)

            padded = torch.zeros(B, ml, 6, device=device)
            mask = torch.ones(B, ml, dtype=torch.bool, device=device)
            for k, (s, e) in enumerate(bc):
                padded[k, :e-s] = torch.from_numpy(all_feat[s:e])
                mask[k, :e-s] = False

            # GPU LD features
            ld = compute_ld_auto(padded.float() if device == 'cuda' else padded, 100)

            # Model
            inp = padded.half() if device == 'cuda' else padded
            ld_inp = ld.half() if device == 'cuda' else ld
            h = encoder(inp, mask)
            logits = pip_head(h, ld_inp, mask)
            scores = torch.sigmoid(logits)

            for k, (s, e) in enumerate(bc):
                all_scores[s:e] = scores[k, :e-s].float().cpu().numpy()

    if device == 'cuda':
        torch.cuda.synchronize()
    elapsed_gpu = time.time() - t0
    elapsed_total = elapsed_prep + elapsed_gpu
    log(f"GPU: {n_total:,} variants in {elapsed_gpu:.1f}s = {n_total/elapsed_gpu:,.0f}/sec")

    # Build output in original order
    snps = df['SNP'].values[sort_idx]
    out_chrs = chrs[sort_idx]
    out_bps = bps[sort_idx]

    result = pd.DataFrame({
        'SNP': snps,
        'CHR': out_chrs,
        'BP': out_bps,
        'SCORE': all_scores,
    })
    result = result.sort_values(['CHR', 'BP'])

    if output_path.endswith('.gz'):
        result.to_csv(output_path, sep='\t', index=False, compression='gzip')
    else:
        result.to_csv(output_path, sep='\t', index=False)

    n_scored = int(np.isfinite(all_scores).sum())
    n_unscored = N - n_scored
    log(f"Done: {n_scored:,} of {N:,} variants scored ({n_unscored} in blocks too small to score), {elapsed_total:.1f}s total")
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='G-Atlas: GWAS → causal variant scores')
    parser.add_argument('input', help='GWAS summary statistics file (TSV)')
    parser.add_argument('output', help='Output scores file (TSV)')
    parser.add_argument('--device', default='cuda' if torch.cuda.is_available() else 'cpu')
    parser.add_argument('--region', default=None, help='Genomic region (e.g., chr6:25000000-35000000)')
    args = parser.parse_args()

    predict_gwas(args.input, args.output, args.device, args.region)
