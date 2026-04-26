#!/usr/bin/env python3
"""
Step 2: Train PIP head on real SuSiE PIPs from CAUSALdb2.
Encoder is FROZEN. Only PIP head trains.

Usage:
    python train/train_pip.py --encoder models/encoder.pt --causaldb /path/to/hm3_credible_sets.parquet --output models/pip_head.pt

Input: pretrained encoder + CAUSALdb2 HM3 credible sets parquet.
Output: trained PIP head weights.
"""

import sys, os, time, argparse
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.model import GAtlasEncoder, PIPHead
from src.features import compute_ld_features, BLOCK_SIZE

BATCH_SIZE = 64
MAX_TRAIN_BLOCKS = 50000
MAX_TEST_BLOCKS = 10000


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def build_blocks(df, gwas_ids, max_blocks):
    """Build training blocks from CAUSALdb2 data."""
    blocks = []
    for gwas_id in gwas_ids:
        if len(blocks) >= max_blocks:
            break
        sub = df[df['meta_id'] == gwas_id]
        if len(sub) < 20:
            continue
        for (chrom, blk), bdf in sub.groupby([sub['chr'], sub['bp'] // BLOCK_SIZE]):
            if len(bdf) < 5 or len(blocks) >= max_blocks:
                continue
            if len(bdf) > 500:
                bdf = bdf.sample(500, random_state=42)
            bdf = bdf.sort_values('bp')
            n = len(bdf)
            feat = np.zeros((n, 6), dtype=np.float32)
            feat[:, 0] = bdf['beta'].values
            feat[:, 1] = bdf['se'].values
            feat[:, 2] = -np.log10(np.clip(bdf['p'].values, 1e-300, 1.0))
            feat[:, 3] = bdf['maf'].values
            feat[:, 4] = (bdf['bp'].values - blk * BLOCK_SIZE) / BLOCK_SIZE
            feat[:, 5] = chrom / 23.0
            feat = np.nan_to_num(feat).astype(np.float32)
            ld = compute_ld_features(feat)
            pip = bdf['susie'].values.astype(np.float32) if 'susie' in bdf.columns else np.zeros(n, np.float32)
            blocks.append({'features': feat, 'ld_feat': ld, 'pip': pip})
    return blocks


def train(args):
    device = args.device
    log("=" * 70)
    log(f"TRAIN PIP HEAD ON REAL SuSiE PIPs ({device})")
    log("=" * 70)

    # Load CAUSALdb2
    log(f"Loading CAUSALdb2: {args.causaldb}")
    cdb = pd.read_parquet(args.causaldb)
    gwas_ids = sorted(cdb['meta_id'].unique())
    log(f"CAUSALdb2: {len(cdb)} variants, {len(gwas_ids)} GWAS")

    # Verify no sentinel UKBB traits in training data
    SENTINEL_PREFIXES = ('ukb-a-', 'ukb-b-', 'ukb-d-')
    leaked = [g for g in gwas_ids if any(g.startswith(p) for p in SENTINEL_PREFIXES)]
    assert len(leaked) == 0, f"LEAKAGE: {len(leaked)} sentinel UKBB traits found in CAUSALdb2: {leaked[:5]}"
    log(f"Leakage check: 0 sentinel UKBB traits in training data")

    # Split
    rng = np.random.default_rng(42)
    rng.shuffle(gwas_ids)
    n_train = int(len(gwas_ids) * 0.8)
    train_ids = gwas_ids[:n_train]
    test_ids = gwas_ids[n_train:]
    log(f"Train: {len(train_ids)}, Test: {len(test_ids)}")

    # Build blocks
    log("Building training blocks...")
    train_blocks = build_blocks(cdb, train_ids, MAX_TRAIN_BLOCKS)
    log(f"Training blocks: {len(train_blocks)}")
    log("Building test blocks...")
    test_blocks = build_blocks(cdb, test_ids, MAX_TEST_BLOCKS)
    log(f"Test blocks: {len(test_blocks)}")

    # Load frozen encoder
    encoder = GAtlasEncoder().to(device)
    encoder.load_state_dict(torch.load(args.encoder, weights_only=True, map_location=device))
    encoder.eval()
    for p in encoder.parameters():
        p.requires_grad = False

    # PIP head
    pip_head = PIPHead().to(device)
    optimizer = torch.optim.Adam(pip_head.parameters(), lr=1e-3)

    best_auc = 0
    for epoch in range(args.epochs):
        pip_head.train()
        rng2 = np.random.default_rng(epoch)
        rng2.shuffle(train_blocks)
        total_loss = 0
        n_batches = 0

        for bi in range(0, len(train_blocks), BATCH_SIZE):
            batch = train_blocks[bi:bi+BATCH_SIZE]
            bs = len(batch)
            ml = max(len(b['features']) for b in batch)
            fp = torch.zeros(bs, ml, 6, device=device)
            lf = torch.zeros(bs, ml, 10, device=device)
            target = torch.zeros(bs, ml, device=device)
            mask = torch.ones(bs, ml, dtype=torch.bool, device=device)

            for k, b in enumerate(batch):
                n = len(b['features'])
                fp[k, :n] = torch.tensor(b['features'])
                lf[k, :n] = torch.tensor(b['ld_feat'])
                target[k, :n] = torch.tensor(b['pip'])
                mask[k, :n] = False

            with torch.no_grad():
                h = encoder(fp, mask)
            logits = pip_head(h, lf, mask)
            pred = torch.sigmoid(logits)
            loss = F.binary_cross_entropy(pred[~mask], target[~mask].clamp(0, 1))

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
            n_batches += 1

            if (bi // BATCH_SIZE + 1) % 200 == 0:
                log(f"    batch {bi//BATCH_SIZE+1}/{len(train_blocks)//BATCH_SIZE} loss={loss.item():.4f}")

        # Eval
        pip_head.eval()
        all_pred, all_true = [], []
        with torch.no_grad():
            for bi in range(0, len(test_blocks), BATCH_SIZE):
                batch = test_blocks[bi:bi+BATCH_SIZE]
                bs = len(batch)
                ml = max(len(b['features']) for b in batch)
                fp = torch.zeros(bs, ml, 6, device=device)
                lf = torch.zeros(bs, ml, 10, device=device)
                target = torch.zeros(bs, ml, device=device)
                mask = torch.ones(bs, ml, dtype=torch.bool, device=device)
                for k, b in enumerate(batch):
                    n = len(b['features'])
                    fp[k, :n] = torch.tensor(b['features'])
                    lf[k, :n] = torch.tensor(b['ld_feat'])
                    target[k, :n] = torch.tensor(b['pip'])
                    mask[k, :n] = False
                h = encoder(fp, mask)
                logits = pip_head(h, lf, mask)
                pred = torch.sigmoid(logits)
                valid = ~mask
                all_pred.extend(pred[valid].cpu().numpy())
                all_true.extend(target[valid].cpu().numpy())

        binary_true = (np.array(all_true) > 0.5).astype(int)
        if binary_true.sum() > 0 and binary_true.sum() < len(binary_true):
            auc = roc_auc_score(binary_true, all_pred)
        else:
            auc = 0
        spearman = np.corrcoef(np.argsort(all_pred).argsort(), np.argsort(all_true).argsort())[0, 1]

        log(f"  Epoch {epoch+1}: loss={total_loss/max(n_batches,1):.4f} AUC={auc:.4f} Sp={spearman:+.4f}")

        if auc > best_auc:
            best_auc = auc
            torch.save(pip_head.state_dict(), args.output)
            log(f"  → saved best (AUC={auc:.4f})")

    log(f"\nBest AUC: {best_auc:.4f}")
    log(f"Saved: {args.output}")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--encoder', required=True, help='Pretrained encoder .pt')
    p.add_argument('--causaldb', required=True, help='CAUSALdb2 HM3 parquet')
    p.add_argument('--output', default='models/pip_head.pt')
    p.add_argument('--device', default='cuda' if torch.cuda.is_available() else 'cpu')
    p.add_argument('--epochs', type=int, default=5)
    train(p.parse_args())
