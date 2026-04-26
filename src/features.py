"""
Feature computation for G-Atlas.

Base features (6 per variant, from GWAS summary statistics):
  0: beta (effect size)
  1: SE (standard error)
  2: -log10(p-value)
  3: MAF (effect allele frequency)
  4: relative position within 500kb block (0-1)
  5: chromosome / 23

LD features (10 per variant, computed from base features):
  0: nearby variant count (normalized)
  1: local variant density
  2: local LD score (sum chi² in window)
  3: max -log10p in window (normalized)
  4: mean |beta| in window
  5: count of significant variants in window (normalized)
  6: chi² (normalized)
  7: z-score (normalized)
  8: effective sample size (normalized)
  9: MAF

GPU acceleration: auto-detects CUDA → uses PyTorch GPU ops.
Falls back to CPU numpy if no GPU available.
"""

import numpy as np

try:
    import torch
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False

# CuPy CUDA kernel for maximum speed
HAS_CUPY = False
_cupy_ld_kernel = None
try:
    import cupy as cp
    _cupy_ld_kernel = cp.RawKernel(r'''
    extern "C" __global__
    void ld_features(
        const float* feat, float* ld_out,
        const int B, const int N, const float window_kb
    ) {
        int b = blockIdx.x;
        int i = threadIdx.x + blockDim.x * blockIdx.y;
        if (b >= B || i >= N) return;
        int base = b * N * 6;
        float beta_i = feat[base+i*6+0];
        float se_i = fmaxf(feat[base+i*6+1], 1e-8f);
        float nlp_i = feat[base+i*6+2];
        float maf_i = feat[base+i*6+3];
        float pos_i = feat[base+i*6+4] * window_kb * 10.0f;
        float z_i = beta_i / se_i;
        float chi2_i = z_i * z_i;
        float nearby=0, sum_chi2=0, max_nlp=-1e9f, sum_ab=0, cnt_sig=0;
        float mc2=0, maz=0, mnl=0;
        for (int j=0; j<N; j++) {
            float zj = feat[base+j*6+0] / fmaxf(feat[base+j*6+1], 1e-8f);
            float c2=zj*zj; float nj=feat[base+j*6+2];
            if (c2>mc2) mc2=c2; if (fabsf(zj)>maz) maz=fabsf(zj); if (nj>mnl) mnl=nj;
        }
        mc2=fmaxf(mc2,1e-8f); maz=fmaxf(maz,1e-8f); mnl=fmaxf(mnl,1e-8f);
        for (int j=0; j<N; j++) {
            float pj = feat[base+j*6+4] * window_kb * 10.0f;
            if (fabsf(pos_i-pj) < window_kb) {
                nearby+=1; float zj=feat[base+j*6+0]/fmaxf(feat[base+j*6+1],1e-8f);
                sum_chi2+=zj*zj; float nj=feat[base+j*6+2];
                if (nj>max_nlp) max_nlp=nj;
                sum_ab+=fabsf(feat[base+j*6+0]); if (nj>5.0f) cnt_sig+=1;
            }
        }
        int ob = b*N*10;
        ld_out[ob+i*10+0]=nearby/fmaxf((float)N,1.0f);
        ld_out[ob+i*10+1]=nearby/fmaxf(window_kb,1.0f);
        ld_out[ob+i*10+2]=sum_chi2/fmaxf((float)N,1.0f);
        ld_out[ob+i*10+3]=(max_nlp>-1e8f?max_nlp:0.0f)/mnl;
        ld_out[ob+i*10+4]=(nearby>0)?sum_ab/nearby:0.0f;
        ld_out[ob+i*10+5]=cnt_sig/fmaxf((float)N,1.0f);
        ld_out[ob+i*10+6]=chi2_i/mc2;
        ld_out[ob+i*10+7]=z_i/maz;
        float neff=1.0f/(se_i*se_i*2.0f*fmaxf(maf_i*(1.0f-maf_i),1e-8f));
        ld_out[ob+i*10+8]=log1pf(neff)/15.0f;
        ld_out[ob+i*10+9]=maf_i;
    }
    ''', 'ld_features')
    HAS_CUPY = True
except Exception:
    pass

BLOCK_SIZE = 500_000
N_LD_FEATURES = 10


def build_features(df_block, chrom, block_start):
    """Build 6-feature array from a DataFrame of variants in one block."""
    n = len(df_block)
    feat = np.zeros((n, 6), dtype=np.float32)
    feat[:, 0] = df_block['BETA'].values if 'BETA' in df_block else df_block['beta'].values
    feat[:, 1] = df_block['SE'].values if 'SE' in df_block else df_block['standard_error'].values

    p_col = 'P' if 'P' in df_block else 'p_value'
    feat[:, 2] = -np.log10(np.clip(df_block[p_col].values.astype(float), 1e-300, 1.0))

    maf_col = 'MAF' if 'MAF' in df_block else 'effect_allele_frequency'
    feat[:, 3] = df_block[maf_col].values if maf_col in df_block else 0.5

    bp_col = 'BP' if 'BP' in df_block else 'base_pair_location'
    feat[:, 4] = (df_block[bp_col].values - block_start) / BLOCK_SIZE
    feat[:, 5] = chrom / 23.0

    return np.nan_to_num(feat).astype(np.float32)


def compute_ld_features(base_features, window_kb=100):
    """Compute 10 LD-aware features. CPU numpy version (fallback)."""
    n = len(base_features)
    ld = np.zeros((n, N_LD_FEATURES), dtype=np.float32)

    beta = base_features[:, 0]
    se = np.clip(base_features[:, 1], 1e-8, None)
    z = beta / se
    chi2 = z ** 2
    neg_log_p = base_features[:, 2]
    maf = base_features[:, 3]
    positions = base_features[:, 4] * window_kb * 10

    max_chi2 = max(chi2.max(), 1e-8)
    max_absz = max(np.abs(z).max(), 1e-8)
    max_nlp = max(neg_log_p.max(), 1e-8)

    for i in range(n):
        dists = np.abs(positions - positions[i])
        in_window = dists < window_kb
        nearby = np.sum(in_window)

        ld[i, 0] = nearby / max(n, 1)
        ld[i, 1] = nearby / max(window_kb, 1)
        ld[i, 2] = np.sum(chi2[in_window]) / max(n, 1)
        ld[i, 3] = np.max(neg_log_p[in_window]) / max_nlp
        ld[i, 4] = np.mean(np.abs(beta[in_window]))
        ld[i, 5] = np.sum((neg_log_p > 5) & in_window) / max(n, 1)
        ld[i, 6] = chi2[i] / max_chi2
        ld[i, 7] = z[i] / max_absz
        n_eff = 1.0 / (se[i]**2 * 2 * max(maf[i] * (1-maf[i]), 1e-8))
        ld[i, 8] = np.log1p(n_eff) / 15.0
        ld[i, 9] = maf[i]

    return np.nan_to_num(ld).astype(np.float32)


def compute_ld_gpu(feat_batch, window_kb=100):
    """GPU-accelerated LD features using PyTorch. (B, N, 6) → (B, N, 10)."""
    B, N, _ = feat_batch.shape
    beta = feat_batch[:, :, 0]
    se = feat_batch[:, :, 1].clamp(min=1e-8)
    z = beta / se
    chi2 = z ** 2
    nlp = feat_batch[:, :, 2]
    maf = feat_batch[:, :, 3]
    pos = feat_batch[:, :, 4] * window_kb * 10

    dists = (pos.unsqueeze(2) - pos.unsqueeze(1)).abs()
    iw = dists < window_kb
    nc = iw.float().sum(dim=2)

    ld = torch.zeros(B, N, 10, device=feat_batch.device)
    ld[:, :, 0] = nc / max(N, 1)
    ld[:, :, 1] = nc / max(window_kb, 1)
    ld[:, :, 2] = (iw.float() * chi2.unsqueeze(1)).sum(2) / max(N, 1)
    mnlp = nlp.max(1, keepdim=True).values.clamp(min=1e-8)
    ld[:, :, 3] = torch.where(iw, nlp.unsqueeze(1).expand_as(iw),
                               torch.tensor(-1e9, device=feat_batch.device)).max(2).values / mnlp
    ld[:, :, 4] = torch.where(iw, beta.abs().unsqueeze(1).expand_as(iw),
                               torch.zeros_like(iw, dtype=torch.float)).sum(2) / nc.clamp(min=1)
    ld[:, :, 5] = (iw & (nlp > 5).unsqueeze(1)).float().sum(2) / max(N, 1)
    ld[:, :, 6] = chi2 / chi2.max(1, keepdim=True).values.clamp(min=1e-8)
    ld[:, :, 7] = z / z.abs().max(1, keepdim=True).values.clamp(min=1e-8)
    neff = 1.0 / (se**2 * 2 * (maf * (1 - maf)).clamp(min=1e-8))
    ld[:, :, 8] = torch.log1p(neff) / 15.0
    ld[:, :, 9] = maf
    return torch.nan_to_num(ld)


def compute_ld_cupy(feat_tensor, window_kb=100.0):
    """CuPy CUDA kernel — fastest. (B, N, 6) torch CUDA → (B, N, 10) torch CUDA."""
    B, N, _ = feat_tensor.shape
    feat_cp = cp.asarray(feat_tensor.contiguous().float())
    out_cp = cp.zeros((B, N, 10), dtype=cp.float32)
    threads = min(N, 512)
    blocks_y = (N + threads - 1) // threads
    _cupy_ld_kernel((B, blocks_y), (threads,),
                    (feat_cp, out_cp, np.int32(B), np.int32(N), np.float32(window_kb)))
    return torch.as_tensor(out_cp, device=feat_tensor.device)


def compute_ld_auto(feat_batch, window_kb=100):
    """Auto-detect best backend: CuPy CUDA → PyTorch GPU → CPU numpy."""
    if HAS_CUPY and HAS_TORCH and isinstance(feat_batch, torch.Tensor) and feat_batch.is_cuda:
        return compute_ld_cupy(feat_batch, window_kb)
    elif HAS_TORCH and isinstance(feat_batch, torch.Tensor) and feat_batch.is_cuda:
        return compute_ld_gpu(feat_batch, window_kb)
    elif HAS_TORCH and isinstance(feat_batch, torch.Tensor):
        return compute_ld_gpu(feat_batch, window_kb)
    else:
        return compute_ld_features(feat_batch, window_kb)
