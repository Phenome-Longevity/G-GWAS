#!/usr/bin/env python3
"""
Step 1: Self-supervised pretraining on GWAS summary statistics.
Objectives: masked z-score reconstruction + noise-aware denoising.

Usage:
    python train/pretrain.py --gwas_dir /path/to/npz/files --output models/encoder.pt

Input: directory of .npz files, each containing blocks of (N_variants, 6) features.
Output: pretrained encoder weights.
"""

import sys, os, time, glob, argparse
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.model import GAtlasEncoder

MASK_RATIO = 0.15
SAVE_EVERY = 5000


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


class ReconstructionHead(nn.Module):
    def __init__(self, d_model=128):
        super().__init__()
        self.head = nn.Sequential(nn.Linear(d_model, 64), nn.GELU(), nn.Linear(64, 3))
    def forward(self, x):
        return self.head(x)


class DenoisingHead(nn.Module):
    def __init__(self, d_model=128):
        super().__init__()
        self.head = nn.Sequential(nn.Linear(d_model, 64), nn.GELU(), nn.Linear(64, 1))
    def forward(self, x):
        return self.head(x).squeeze(-1)


def train(args):
    device = args.device
    log("=" * 70)
    log(f"PRETRAIN ENCODER ({device})")
    log("=" * 70)

    # Preload all blocks into RAM
    npz_files = sorted(glob.glob(os.path.join(args.gwas_dir, '*.npz')))
    # Verify no UKBB files in pretraining (all should be EBI GCST* or similar)
    ukbb_files = [f for f in npz_files if 'ukb-' in os.path.basename(f).lower() or 'ukb_' in os.path.basename(f).lower()]
    assert len(ukbb_files) == 0, f"LEAKAGE: {len(ukbb_files)} UKBB files in pretrain dir: {ukbb_files[:5]}"
    log(f"Found {len(npz_files)} .npz files (0 UKBB, leakage check passed). Preloading...")

    all_blocks = []
    for i, f in enumerate(npz_files):
        try:
            data = np.load(f, allow_pickle=True)
            n_blocks = int(data.get('n_blocks', [0])[0]) if 'n_blocks' in data else 0
            for j in range(n_blocks):
                key = f'block_{j}'
                if key in data:
                    b = data[key]
                    if b.shape[0] >= 20 and b.shape[1] >= 6:
                        all_blocks.append(b.astype(np.float32))
        except:
            pass
        if (i + 1) % 2000 == 0:
            log(f"  {i+1}/{len(npz_files)}, {len(all_blocks):,} blocks")

    log(f"Preloaded {len(all_blocks):,} blocks ({sum(b.nbytes for b in all_blocks)/1e9:.1f} GB)")

    # Model
    encoder = GAtlasEncoder().to(device)
    if args.checkpoint and os.path.exists(args.checkpoint):
        encoder.load_state_dict(torch.load(args.checkpoint, weights_only=True, map_location=device))
        log(f"Loaded checkpoint: {args.checkpoint}")

    recon_head = ReconstructionHead().to(device)
    denoise_head = DenoisingHead().to(device)

    params = list(encoder.parameters()) + list(recon_head.parameters()) + list(denoise_head.parameters())
    optimizer = torch.optim.AdamW(params, lr=args.lr, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.steps)

    log(f"Params: {sum(p.numel() for p in params):,}")
    log(f"Steps: {args.steps}, batch: {args.batch_size}, lr: {args.lr}\n")

    rng = np.random.default_rng(42)
    n_blocks = len(all_blocks)
    t0 = time.time()
    nan_count = 0

    try:
        for step in range(args.steps):
            encoder.train(); recon_head.train(); denoise_head.train()

            idx = rng.integers(0, n_blocks, size=args.batch_size)
            blocks = [torch.tensor(all_blocks[i]) for i in idx]
            max_len = max(len(b) for b in blocks)
            feat = torch.zeros(len(blocks), max_len, 6, device=device)
            mask = torch.ones(len(blocks), max_len, dtype=torch.bool, device=device)
            for i, b in enumerate(blocks):
                feat[i, :len(b)] = b.to(device)
                mask[i, :len(b)] = False
            valid = ~mask

            # Masked reconstruction
            rand_mask = (torch.rand_like(feat[:,:,0]) < MASK_RATIO) & valid
            feat_masked = feat.clone()
            feat_masked[:,:,0][rand_mask] = 0
            feat_masked[:,:,1][rand_mask] = 0
            feat_masked[:,:,2][rand_mask] = 0

            hidden = encoder(feat_masked, mask)
            recon = recon_head(hidden)
            loss_recon = (
                F.mse_loss(recon[:,:,0][rand_mask], feat[:,:,0][rand_mask]) +
                F.mse_loss(recon[:,:,1][rand_mask], feat[:,:,1][rand_mask]) +
                F.mse_loss(recon[:,:,2][rand_mask], feat[:,:,2][rand_mask])
            ) / 3

            # Denoising
            se_safe = feat[:,:,1].clamp(min=1e-6)
            clean_z = feat[:,:,0] / se_safe
            clean_z = clean_z.clamp(-30, 30)
            noise_factor = 1 + torch.rand_like(se_safe) * 2
            noisy_se = feat[:,:,1] * noise_factor
            noisy_beta = feat[:,:,0] + torch.randn_like(feat[:,:,0]) * noisy_se * 0.5
            feat_noisy = feat.clone()
            feat_noisy[:,:,0] = noisy_beta
            feat_noisy[:,:,1] = noisy_se

            hidden_noisy = encoder(feat_noisy, mask)
            pred_clean_z = denoise_head(hidden_noisy)
            loss_denoise = F.mse_loss(pred_clean_z[valid], clean_z[valid])

            loss = loss_recon + 0.5 * loss_denoise

            if not torch.isfinite(loss):
                nan_count += 1
                if nan_count > 100:
                    log(f"Too many NaN losses ({nan_count}), stopping early")
                    break
                continue

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(params, 0.5)
            optimizer.step()
            scheduler.step()

            if (step + 1) % 100 == 0:
                elapsed = time.time() - t0
                rate = (step + 1) / elapsed
                eta_min = (args.steps - step - 1) / rate / 60
                log(f"Step {step+1}/{args.steps} | loss={loss.item():.4f} "
                    f"(recon={loss_recon.item():.4f} denoise={loss_denoise.item():.4f}) "
                    f"| {rate:.1f}/s ETA {eta_min:.0f}min")

            if (step + 1) % SAVE_EVERY == 0:
                torch.save(encoder.state_dict(), args.output.replace('.pt', f'_step{step+1}.pt'))
                torch.save(encoder.state_dict(), args.output)
                log(f"  Checkpoint saved: step {step+1}")

    except Exception as e:
        log(f"ERROR during training: {e}")
    finally:
        torch.save(encoder.state_dict(), args.output)
        log(f"Saved: {args.output} (step ~{step+1})")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--gwas_dir', required=True, help='Directory of .npz GWAS blocks')
    p.add_argument('--output', default='models/encoder.pt')
    p.add_argument('--checkpoint', default=None, help='Resume from checkpoint')
    p.add_argument('--device', default='cuda' if torch.cuda.is_available() else 'cpu')
    p.add_argument('--steps', type=int, default=50000)
    p.add_argument('--batch_size', type=int, default=32)
    p.add_argument('--lr', type=float, default=1e-4)
    train(p.parse_args())
