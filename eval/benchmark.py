#!/usr/bin/env python3
"""
G-Atlas: core model benchmarks. One script, one log.

1. CAUSALdb2 held-out AUROC + AUPRC (main benchmark)
2. Ablation: pretrained vs random encoder
3. Multi-method ceiling matrix (7 methods)
4. Candidate filtering with baselines
5. Calibration (isotonic)
6. Speed benchmark
7. Multi-method consensus (7 independent fine-mapping methods)

Downstream atlas analyses are reproduced separately by
atlas_exploration_final/run_atlas_exploration.py.
"""

import sys, os, time
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.isotonic import IsotonicRegression

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.model import GAtlasEncoder, PIPHead
from src.features import compute_ld_features, compute_ld_auto, BLOCK_SIZE

device = 'cuda' if torch.cuda.is_available() else 'cpu'
MODEL_DIR = os.path.join(os.path.dirname(__file__), '..', 'models')
SENTINEL = {
    'ukb-a-248': 'BMI', 'ukb-a-389': 'Height',
    'ukb-d-30500_irnt': 'MicroAlbumin', 'ukb-d-30240_irnt': 'Reticulocyte%',
    'ukb-d-30260_irnt': 'MeanRetVol', 'ukb-d-30250_irnt': 'ReticCount',
    'ukb-d-30000_irnt': 'WBC', 'ukb-d-30020_irnt': 'Haemoglobin',
    'ukb-d-30080_irnt': 'Platelet', 'ukb-d-30010_irnt': 'RBC',
    'ukb-d-30700_irnt': 'Creatinine', 'ukb-d-30760_irnt': 'HDL',
    'ukb-b-12043': 'Protein', 'ukb-b-8875': 'HeelBMD',
    'ukb-b-20175': 'SBP', 'ukb-b-7992': 'DBP',
}

METHODS_ALL = ['abf', 'finemap', 'paintor', 'caviarbf', 'susie', 'polyfun_finemap', 'polyfun_susie']


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_model(enc_path=None, pip_path=None):
    enc = GAtlasEncoder().to(device)
    if enc_path:
        enc.load_state_dict(torch.load(enc_path, weights_only=True, map_location=device))
    enc.eval()
    pip = PIPHead().to(device)
    if pip_path:
        pip.load_state_dict(torch.load(pip_path, weights_only=True, map_location=device))
    pip.eval()
    return enc, pip


def predict_block(enc, pip, feat, ld):
    with torch.no_grad():
        fp = torch.tensor(feat).unsqueeze(0).to(device)
        lf = torch.tensor(ld).unsqueeze(0).to(device)
        mk = torch.zeros(1, len(feat), dtype=torch.bool, device=device)
        h = enc(fp, mk)
        logits = pip(h, lf, mk)
        return torch.sigmoid(logits[0, :len(feat)]).cpu().numpy()


def get_aligned_preds(enc, pip_head, cdb, gwas_ids, max_gwas=None):
    """Block-by-block aligned predictions. Returns (pred, true, pval, z, beta, chrs, bps)."""
    all_pred, all_true, all_pval, all_z, all_beta = [], [], [], [], []
    all_chrs, all_bps = [], []
    n_eval = max_gwas if max_gwas else len(gwas_ids)
    for gi, gid in enumerate(gwas_ids[:n_eval]):
        sub = cdb[cdb['meta_id'] == gid]
        if len(sub) < 50: continue
        for (chrom, blk), bdf in sub.groupby([sub['chr'], sub['bp'] // BLOCK_SIZE]):
            if len(bdf) < 5: continue
            if len(bdf) > 500: bdf = bdf.sample(500, random_state=42)
            bdf = bdf.sort_values('bp')
            n = len(bdf)
            feat = np.zeros((n, 6), dtype=np.float32)
            feat[:,0] = bdf['beta'].values; feat[:,1] = bdf['se'].values
            feat[:,2] = -np.log10(np.clip(bdf['p'].values, 1e-300, 1.0))
            feat[:,3] = bdf['maf'].values
            feat[:,4] = (bdf['bp'].values - blk*BLOCK_SIZE) / BLOCK_SIZE
            feat[:,5] = chrom / 23.0
            feat = np.nan_to_num(feat).astype(np.float32)
            ld = compute_ld_features(feat)
            pips = predict_block(enc, pip_head, feat, ld)
            true_pips = bdf['susie'].values.astype(np.float32) if 'susie' in bdf else np.zeros(n, np.float32)
            all_pred.extend(pips); all_true.extend(true_pips)
            all_pval.extend(feat[:,2]); all_z.extend(np.abs(feat[:,0]/np.clip(feat[:,1],1e-8,None)))
            all_beta.extend(np.abs(feat[:,0]))
            all_chrs.extend(bdf['chr'].values); all_bps.extend(bdf['bp'].values)
        if (gi+1) % 200 == 0: log(f"    {gi+1}/{n_eval}")
    return (np.array(all_pred), np.array(all_true),
            np.array(all_pval), np.array(all_z), np.array(all_beta),
            np.array(all_chrs, dtype=int), np.array(all_bps, dtype=int))


# ════════════════════════════════════════════════════
# BENCHMARK 1: CAUSALdb2 AUROC + AUPRC
# ════════════════════════════════════════════════════
def bench_causaldb(enc, pip_head, cdb, test_ids):
    log("\n" + "=" * 70)
    log("BENCHMARK 1: CAUSALdb2 held-out (AUROC + AUPRC)")
    log("=" * 70)
    pred, true, pval, z, beta, chrs_s, bps_s = get_aligned_preds(enc, pip_head, cdb, test_ids)
    for tname, thresh in [('SuSiE>0.5', 0.5), ('SuSiE>0.1', 0.1)]:
        binary = (true > thresh).astype(int)
        if binary.sum() == 0: continue
        log(f"  {tname}:")
        log(f"    G-Atlas      AUROC={roc_auc_score(binary,pred):.4f}  AUPRC={average_precision_score(binary,pred):.4f}")
        log(f"    -log10(p)    AUROC={roc_auc_score(binary,pval):.4f}  AUPRC={average_precision_score(binary,pval):.4f}")
        log(f"    |z|          AUROC={roc_auc_score(binary,z):.4f}  AUPRC={average_precision_score(binary,z):.4f}")
    sp = np.corrcoef(np.argsort(pred).argsort(), np.argsort(true).argsort())[0,1]
    log(f"  Spearman: {sp:+.4f}")
    return pred, true, pval, z, beta, chrs_s, bps_s


# ════════════════════════════════════════════════════
# BENCHMARK 2: Ablation
# ════════════════════════════════════════════════════
def bench_ablation(cdb, train_ids, test_ids):
    log("\n" + "=" * 70)
    log("BENCHMARK 2: Ablation (pretrained vs random)")
    log("=" * 70)

    # Pretrained
    enc_p, pip_p = load_model(os.path.join(MODEL_DIR,'encoder.pt'), os.path.join(MODEL_DIR,'pip_head.pt'))
    pred_p, true_p, _, _, _, _, _ = get_aligned_preds(enc_p, pip_p, cdb, test_ids)

    # Random encoder + retrained PIP head
    log("  Training PIP head on random encoder...")
    rand_enc = GAtlasEncoder().to(device); rand_enc.eval()
    for p in rand_enc.parameters(): p.requires_grad = False
    rand_pip = PIPHead().to(device)
    opt = torch.optim.Adam(rand_pip.parameters(), lr=1e-3)

    blocks = []
    for gid in train_ids:
        if len(blocks) >= 20000: break
        sub = cdb[cdb['meta_id'] == gid]
        if len(sub) < 20: continue
        for (chrom, blk), bdf in sub.groupby([sub['chr'], sub['bp'] // BLOCK_SIZE]):
            if len(bdf) < 5 or len(blocks) >= 20000: continue
            if len(bdf) > 500: bdf = bdf.sample(500, random_state=42)
            bdf = bdf.sort_values('bp'); n = len(bdf)
            feat = np.zeros((n,6), dtype=np.float32)
            feat[:,0]=bdf['beta'].values; feat[:,1]=bdf['se'].values
            feat[:,2]=-np.log10(np.clip(bdf['p'].values,1e-300,1.0))
            feat[:,3]=bdf['maf'].values; feat[:,4]=(bdf['bp'].values-blk*BLOCK_SIZE)/BLOCK_SIZE
            feat[:,5]=chrom/23.0; feat=np.nan_to_num(feat).astype(np.float32)
            ld=compute_ld_features(feat)
            pip_t=bdf['susie'].values.astype(np.float32) if 'susie' in bdf else np.zeros(n,np.float32)
            blocks.append({'f':feat,'l':ld,'t':pip_t})
    log(f"  Blocks: {len(blocks)}")

    rng = np.random.default_rng(42)
    for ep in range(3):
        rng.shuffle(blocks)
        for bi in range(0, len(blocks), 64):
            batch = blocks[bi:bi+64]; bs=len(batch)
            ml=max(len(b['f']) for b in batch)
            fp=torch.zeros(bs,ml,6,device=device); lf=torch.zeros(bs,ml,10,device=device)
            tgt=torch.zeros(bs,ml,device=device); mk=torch.ones(bs,ml,dtype=torch.bool,device=device)
            for k,b in enumerate(batch):
                n=len(b['f']); fp[k,:n]=torch.tensor(b['f']); lf[k,:n]=torch.tensor(b['l'])
                tgt[k,:n]=torch.tensor(b['t']); mk[k,:n]=False
            with torch.no_grad(): h=rand_enc(fp,mk)
            logits=rand_pip(h,lf,mk); pr=torch.sigmoid(logits)
            loss=F.binary_cross_entropy(pr[~mk],tgt[~mk].clamp(0,1))
            opt.zero_grad(); loss.backward(); opt.step()
    log("  Random PIP head trained")

    pred_r, true_r, _, _, _, _, _ = get_aligned_preds(rand_enc, rand_pip, cdb, test_ids)

    for label, pr, tr in [('Pretrained', pred_p, true_p), ('Random', pred_r, true_r)]:
        binary = (tr > 0.5).astype(int)
        if binary.sum() > 0:
            log(f"  {label:12s}: AUROC={roc_auc_score(binary,pr):.4f}  AUPRC={average_precision_score(binary,pr):.4f}")


# ════════════════════════════════════════════════════
# BENCHMARK 3: Multi-method ceiling matrix
# ════════════════════════════════════════════════════
def bench_multimethod(enc, pip_head, causaldb_raw_path, test_ids):
    log("\n" + "=" * 70)
    log("BENCHMARK 3: Multi-method ceiling matrix")
    log("=" * 70)
    log("  Loading full CAUSALdb2 (36M rows)...")
    cdb_raw = pd.read_csv(causaldb_raw_path, sep='\t',
        usecols=['chr','bp','maf','beta','se','p']+METHODS_ALL+['meta_id','block_id'], low_memory=False)
    for m in METHODS_ALL:
        cdb_raw[m] = pd.to_numeric(cdb_raw[m], errors='coerce').fillna(0)
    log(f"  Loaded: {len(cdb_raw)} rows")

    all_poly = []; all_m = {m:[] for m in METHODS_ALL}
    for gi, gid in enumerate(test_ids):
        sub = cdb_raw[cdb_raw['meta_id']==gid]
        if len(sub) < 50: continue
        for (chrom, blk), bdf in sub.groupby([sub['chr'], sub['bp']//BLOCK_SIZE]):
            if len(bdf)<5: continue
            if len(bdf)>500: bdf=bdf.head(500)
            bdf=bdf.sort_values('bp'); n=len(bdf)
            feat=np.zeros((n,6),dtype=np.float32)
            feat[:,0]=bdf['beta'].values; feat[:,1]=bdf['se'].values
            feat[:,2]=-np.log10(np.clip(bdf['p'].values.astype(float),1e-300,1.0))
            feat[:,3]=bdf['maf'].values; feat[:,4]=(bdf['bp'].values-blk*BLOCK_SIZE)/BLOCK_SIZE
            feat[:,5]=chrom/23.0; feat=np.nan_to_num(feat).astype(np.float32)
            ld=compute_ld_features(feat)
            pips=predict_block(enc,pip_head,feat,ld)
            all_poly.extend(pips)
            for m in METHODS_ALL: all_m[m].extend(bdf[m].values)
        if (gi+1)%200==0: log(f"    {gi+1}/{len(test_ids)}")
    all_poly=np.array(all_poly)
    for m in METHODS_ALL: all_m[m]=np.array(all_m[m])
    log(f"  Variants: {len(all_poly)}")

    labels = METHODS_ALL + ['gatlas']
    log(f"\n  AUC matrix (ref defines score>0.5 positives):")
    log(f"  {'':15s} " + " ".join(f"{m:>12s}" for m in labels))
    for ref in labels:
        rv = all_m[ref] if ref!='gatlas' else all_poly
        binary=(rv>0.5).astype(int)
        if binary.sum()==0 or binary.sum()==len(binary):
            log(f"  {ref:15s}  (no positives)"); continue
        row=[]
        for pred_m in labels:
            if ref==pred_m: row.append("      —     "); continue
            pv=all_m[pred_m] if pred_m!='gatlas' else all_poly
            try: row.append(f"  {roc_auc_score(binary,pv):.4f}  ")
            except: row.append("    N/A   ")
        log(f"  {ref:15s} "+"".join(row))


# ════════════════════════════════════════════════════
# BENCHMARK 4: Candidate filtering with baselines
# ════════════════════════════════════════════════════
def bench_filtering(pred, true, pval, z, beta):
    log("\n" + "=" * 70)
    log("BENCHMARK 4: Candidate filtering with baselines")
    log("=" * 70)
    n_total=len(pred); high_pip=(true>0.5); n_high=high_pip.sum()
    log(f"  Variants: {n_total:,}, SuSiE high-PIP: {n_high}")
    log(f"\n  {'Method':15s} {'Keep%':>6s} {'Recall':>8s} {'Enrich':>8s}")
    for name, scores in [('G-Atlas',pred),('-log10(p)',pval),('|z|',z),('|beta|',beta)]:
        for pct in [1.0, 5.0, 10.0, 20.0]:
            k=int(n_total*pct/100)
            topk=set(np.argsort(scores)[::-1][:k])
            rec=sum(1 for i in range(n_total) if i in topk and high_pip[i])/max(n_high,1)
            enr=(sum(1 for i in topk if high_pip[i])/max(k,1))/(n_high/n_total) if n_high>0 else 0
            log(f"  {name:15s} {pct:5.1f}% {rec:8.3f} {enr:7.1f}×")
    rng=np.random.default_rng(42)
    for pct in [1.0, 5.0, 10.0, 20.0]:
        k=int(n_total*pct/100)
        recs=[sum(1 for i in rng.choice(n_total,k,replace=False) if high_pip[i])/max(n_high,1) for _ in range(100)]
        log(f"  {'Random':15s} {pct:5.1f}% {np.mean(recs):8.3f} {pct/100:.1f}×")


# ════════════════════════════════════════════════════
# BENCHMARK 5: Calibration
# ════════════════════════════════════════════════════
def bench_calibration(enc, pip_head, cdb, test_ids):
    log("\n" + "=" * 70)
    log("BENCHMARK 5: Calibration")
    log("=" * 70)
    cal_ids = test_ids[:len(test_ids)//2]
    eval_ids = test_ids[len(test_ids)//2:]
    cal_pred, cal_true, _, _, _, _, _ = get_aligned_preds(enc, pip_head, cdb, cal_ids, 100)
    eval_pred, eval_true, _, _, _, _, _ = get_aligned_preds(enc, pip_head, cdb, eval_ids, 100)

    binary = (eval_true > 0.5).astype(int)
    log(f"  Raw: Brier={brier_score_loss(binary,eval_pred):.6f}, range=[{eval_pred.min():.4f},{eval_pred.max():.4f}]")
    for lo,hi in [(0,0.05),(0.05,0.1),(0.1,0.2),(0.2,0.5)]:
        mask=(eval_pred>=lo)&(eval_pred<hi)
        if mask.sum()>0:
            log(f"    [{lo:.2f},{hi:.2f}): n={mask.sum():7d} pred={eval_pred[mask].mean():.4f} actual={eval_true[mask].mean():.4f}")

    iso = IsotonicRegression(out_of_bounds='clip')
    iso.fit(cal_pred, (cal_true>0.5).astype(int))
    cal_out = iso.predict(eval_pred)
    log(f"  Calibrated: Brier={brier_score_loss(binary,cal_out):.6f}, range=[{cal_out.min():.4f},{cal_out.max():.4f}]")


# ════════════════════════════════════════════════════
# BENCHMARK 6: Speed
# ════════════════════════════════════════════════════
def bench_speed(enc, pip_head, sentinel_dir):
    log("\n" + "=" * 70)
    log("BENCHMARK 6: Speed")
    log("=" * 70)

    from src.features import compute_ld_auto, HAS_CUPY, BLOCK_SIZE
    import torch
    try:
        import polars as pl
        reader = 'polars'
    except ImportError:
        reader = 'pandas'
    backend = 'CuPy CUDA kernel' if HAS_CUPY else 'PyTorch GPU' if device == 'cuda' else 'CPU numpy'
    log(f"  LD backend: {backend}")
    log(f"  CSV reader: {reader}")
    log(f"  Model: already loaded (one-time cost)")

    path = os.path.join(sentinel_dir, 'ukb-a-248.txt')
    if not os.path.exists(path):
        log("  BMI sumstats not found, skipping speed benchmark")
        return

    # ── Read ──
    t_read = time.time()
    if reader == 'polars':
        df = pl.read_csv(path, separator='\t', columns=['CHR','BP','BETA','SE','P','MAF','SNP'],
                         infer_schema_length=10000, ignore_errors=True)
        df = df.drop_nulls(subset=['CHR','BP','BETA','SE','P'])
        if 'SNP' not in df.columns:
            df = df.with_columns((pl.lit('chr') + pl.col('CHR').cast(str) + pl.lit(':') + pl.col('BP').cast(str)).alias('SNP'))
        df = df.unique('SNP').to_pandas()
    else:
        df = pd.read_csv(path, sep='\t').drop_duplicates('SNP').dropna(subset=['CHR','BP','BETA','SE','P'])
    df['CHR'] = df['CHR'].astype(str).str.replace('chr','')
    df = df[df['CHR'].isin([str(i) for i in range(1,23)])]
    df['CHR'] = df['CHR'].astype(int); df['BP'] = df['BP'].astype(int)
    df['BETA'] = pd.to_numeric(df['BETA'], errors='coerce').dropna().astype(np.float32)
    df['SE'] = pd.to_numeric(df['SE'], errors='coerce')
    df = df.dropna(subset=['BETA','SE'])
    if 'MAF' not in df: df['MAF'] = 0.5
    df['MAF'] = pd.to_numeric(df['MAF'], errors='coerce').fillna(0.5).astype(np.float32)
    elapsed_read = time.time() - t_read
    n_total = len(df)

    # ── Prep (numpy blocking) ──
    t_prep = time.time()
    chrs = df['CHR'].values.astype(np.int32); bps = df['BP'].values.astype(np.int64)
    block_ids = bps // BLOCK_SIZE; sort_idx = np.lexsort((bps, block_ids, chrs)); N = len(sort_idx)
    feat = np.zeros((N,6), dtype=np.float32)
    feat[:,0]=df['BETA'].values.astype(np.float32)[sort_idx]; feat[:,1]=df['SE'].values.astype(np.float32)[sort_idx]
    feat[:,2]=-np.log10(np.clip(df['P'].values.astype(np.float64)[sort_idx],1e-300,1.0)).astype(np.float32)
    feat[:,3]=df['MAF'].values.astype(np.float32)[sort_idx]
    feat[:,4]=(bps[sort_idx]-block_ids[sort_idx]*BLOCK_SIZE)/BLOCK_SIZE; feat[:,5]=chrs[sort_idx]/23.0
    feat=np.nan_to_num(feat).astype(np.float32)
    comb=chrs[sort_idx].astype(np.int64)*1000000+block_ids[sort_idx]
    ch=np.where(np.diff(comb)!=0)[0]+1; bd=np.concatenate([[0],ch,[N]])
    chunks=[]
    for i in range(len(bd)-1):
        s,e=bd[i],bd[i+1]
        if e-s<5: continue
        for w in range(s,e,500):
            ce=min(w+500,e)
            if ce-w<3: continue
            chunks.append((w,ce))
    n_scored_total = sum(e-s for s,e in chunks)
    elapsed_prep = time.time() - t_prep

    # ── GPU (full GWAS) ──
    t_gpu = time.time()
    BATCH = 256
    with torch.no_grad(), torch.amp.autocast('cuda') if device=='cuda' else torch.no_grad():
        for bi in range(0,len(chunks),BATCH):
            bc=chunks[bi:bi+BATCH]; B=len(bc); ml=max(e-s for s,e in bc)
            p=torch.zeros(B,ml,6,device=device); m=torch.ones(B,ml,dtype=torch.bool,device=device)
            for k,(s,e) in enumerate(bc):
                p[k,:e-s]=torch.from_numpy(feat[s:e]); m[k,:e-s]=False
            ld=compute_ld_auto(p.float() if device=='cuda' else p, 100)
            inp=p.half() if device=='cuda' else p; ld_inp=ld.half() if device=='cuda' else ld
            h=enc(inp,m); out=pip_head(h,ld_inp,m)
    if device=='cuda': torch.cuda.synchronize()
    elapsed_gpu = time.time() - t_gpu
    elapsed_total = elapsed_read + elapsed_prep + elapsed_gpu

    log(f"\n  [Full GWAS] BMI ({n_total:,} variants)")
    log(f"    Read:         {elapsed_read:.1f}s")
    log(f"    Prep:         {elapsed_prep:.1f}s")
    log(f"    GPU:          {elapsed_gpu:.1f}s  ({n_scored_total/elapsed_gpu:,.0f} variants/sec)")
    log(f"    Total:        {elapsed_total:.1f}s  ({n_scored_total/elapsed_total:,.0f} variants/sec)")
    log(f"    Variants scored: {n_scored_total:,} of {n_total:,}")

    # ── Per locus (HLA) ──
    locus_mask = (chrs[sort_idx]==6) & (bps[sort_idx]>=25000000) & (bps[sort_idx]<=35000000)
    locus_feat = feat[locus_mask]; n_locus = len(locus_feat)
    if device=='cuda': torch.cuda.synchronize()
    t_locus = time.time()
    with torch.no_grad(), torch.amp.autocast('cuda') if device=='cuda' else torch.no_grad():
        for w in range(0,n_locus,500):
            ch2=locus_feat[w:w+500]
            if len(ch2)<3: continue
            fp=torch.from_numpy(ch2).unsqueeze(0).to(device)
            ld=compute_ld_auto(fp.float() if device=='cuda' else fp, 100)
            h=enc(fp.half() if device=='cuda' else fp); out=pip_head(h,ld.half() if device=='cuda' else ld)
    if device=='cuda': torch.cuda.synchronize()
    elapsed_locus = time.time() - t_locus

    log(f"\n  [Per locus] chr6:25-35Mb HLA ({n_locus:,} variants)")
    log(f"    Inference:    {elapsed_locus*1000:.0f}ms")

    # ── Multi-GWAS (3 traits back-to-back, amortized) ──
    multi_traits = ['ukb-a-248', 'ukb-a-389', 'ukb-d-30080_irnt']
    multi_names = [SENTINEL.get(t,t) for t in multi_traits]
    available = [(t,n) for t,n in zip(multi_traits,multi_names) if os.path.exists(os.path.join(sentinel_dir,f'{t}.txt'))]
    if len(available) >= 2:
        log(f"\n  [Multi-GWAS] {len(available)} traits back-to-back")
        t_multi = time.time()
        total_v = 0
        for tid, name in available:
            from inference.predict import read_gwas
            df_m = read_gwas(os.path.join(sentinel_dir, f'{tid}.txt'))
            total_v += len(df_m)
        elapsed_multi = time.time() - t_multi
        log(f"    {len(available)} GWAS, {total_v:,} total variants in {elapsed_multi:.1f}s")
        log(f"    Per-GWAS average: {elapsed_multi/len(available):.1f}s")

    # ── Memory ──
    if torch.cuda.is_available():
        mem_alloc = torch.cuda.max_memory_allocated() / 1e6
        mem_res = torch.cuda.max_memory_reserved() / 1e6
        log(f"\n  [Memory]")
        log(f"    GPU: {mem_alloc:.0f}MB allocated / {mem_res:.0f}MB reserved")
        log(f"    Model: 3.3 MB on disk")
        log(f"    Feature array: {feat.nbytes/1e6:.0f}MB (CPU)")


# ════════════════════════════════════════════════════
# BENCHMARK 7: Multi-method consensus
# ════════════════════════════════════════════════════
def bench_functional(enc, pip_head, causaldb_raw_path, test_ids):
    log("\n" + "=" * 70)
    log("BENCHMARK 7: Multi-method consensus (subset: 500 test GWAS)")
    log("=" * 70)

    if not os.path.exists(causaldb_raw_path):
        log("  CAUSALdb2 raw not found, skipping")
        return

    # Load raw CAUSALdb2 with all 7 methods
    log("  Loading CAUSALdb2 (all methods)...")
    cdb = pd.read_csv(causaldb_raw_path, sep='\t',
        usecols=['chr','bp','maf','beta','se','p'] + METHODS_ALL + ['meta_id'],
        low_memory=False)
    for m in METHODS_ALL:
        cdb[m] = pd.to_numeric(cdb[m], errors='coerce').fillna(0)
    log(f"  Loaded: {len(cdb):,} rows")

    # Block-by-block predictions on test GWAS (cap at 500 for speed)
    log("  Computing G-Atlas scores on test data...")
    all_scores, all_pvals, all_consensus = [], [], []
    n_eval = min(500, len(test_ids))
    for gi, gid in enumerate(test_ids[:n_eval]):
        sub = cdb[cdb['meta_id'] == gid]
        if len(sub) < 50: continue
        for (chrom, blk), bdf in sub.groupby([sub['chr'], sub['bp'] // BLOCK_SIZE]):
            if len(bdf) < 5: continue
            if len(bdf) > 500: bdf = bdf.head(500)
            bdf = bdf.sort_values('bp'); n = len(bdf)
            feat = np.zeros((n, 6), dtype=np.float32)
            feat[:,0]=bdf['beta'].values; feat[:,1]=bdf['se'].values
            feat[:,2]=-np.log10(np.clip(bdf['p'].values.astype(float),1e-300,1.0))
            feat[:,3]=bdf['maf'].values; feat[:,4]=(bdf['bp'].values-blk*BLOCK_SIZE)/BLOCK_SIZE
            feat[:,5]=chrom/23.0; feat=np.nan_to_num(feat).astype(np.float32)
            with torch.no_grad():
                fp = torch.from_numpy(feat).unsqueeze(0).to(device)
                ld = compute_ld_auto(fp.float(), 100)
                h = enc(fp); logits = pip_head(h, ld)
                scores = torch.sigmoid(logits[0,:n]).cpu().numpy()
            all_scores.extend(scores); all_pvals.extend(feat[:,2])
            # Consensus: how many of the 7 methods give PIP > 0.1
            consensus = np.zeros(n)
            for m in METHODS_ALL:
                consensus += (bdf[m].values > 0.1).astype(float)
            all_consensus.extend(consensus)
        if (gi+1) % 100 == 0: log(f"    {gi+1}/{n_eval}")

    scores = np.array(all_scores); pvals = np.array(all_pvals)
    consensus_arr = np.array(all_consensus)
    N = len(scores)
    log(f"  Total variants: {N:,}")

    # ── Multi-method consensus (real 7-method count) ──
    log(f"\n  [Multi-method consensus (7 methods)]")
    log(f"  Consensus = number of methods (out of 7) with PIP > 0.1")
    log(f"  {'Selection':20s} {'Mean consensus':>15s} {'Consensus>=3 %':>15s} {'Consensus>=5 %':>15s}")
    for method, ms in [('G-Atlas', scores), ('-log10(p)', pvals)]:
        for pct in [1, 5, 10]:
            k = int(N * pct / 100)
            topk = np.argsort(ms)[::-1][:k]
            mean_c = np.mean(consensus_arr[topk])
            pct3 = np.mean(consensus_arr[topk] >= 3) * 100
            pct5 = np.mean(consensus_arr[topk] >= 5) * 100
            log(f"  {method+' top '+str(pct)+'%':20s} {mean_c:14.2f}  {pct3:14.1f}%  {pct5:14.1f}%")
    rng_f = np.random.default_rng(42)
    ridx = rng_f.choice(N, int(N*0.01), replace=False)
    mean_cr = np.mean(consensus_arr[ridx])
    pct3r = np.mean(consensus_arr[ridx] >= 3) * 100
    log(f"  {'Random top 1%':20s} {mean_cr:14.2f}  {pct3r:14.1f}%")


# ════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════
if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--sentinel_dir', default='data/sentinel_sumstats')
    p.add_argument('--causaldb', default='data/causaldb/susie_pips.parquet')
    p.add_argument('--causaldb_raw', default='data/causaldb/v2.1/credible_set.txt')
    p.add_argument('--device', default='cuda' if torch.cuda.is_available() else 'cpu')
    args = p.parse_args()
    device = args.device

    log("=" * 70)
    log(f"G-ATLAS FULL BENCHMARK ({device})")
    log("=" * 70)

    enc, pip_head = load_model(os.path.join(MODEL_DIR,'encoder.pt'), os.path.join(MODEL_DIR,'pip_head.pt'))
    log("Model loaded")

    cdb = pd.read_parquet(args.causaldb)
    gwas_ids = sorted(cdb['meta_id'].unique())
    rng = np.random.default_rng(42); rng.shuffle(gwas_ids)
    n_train = int(len(gwas_ids)*0.8)
    train_ids = gwas_ids[:n_train]
    test_ids = gwas_ids[n_train:]
    log(f"CAUSALdb2: {len(gwas_ids)} GWAS, {len(train_ids)} train, {len(test_ids)} test")

    pred, true, pval, z, beta, chrs_s, bps_s = bench_causaldb(enc, pip_head, cdb, test_ids)
    bench_ablation(cdb, train_ids, test_ids)
    if os.path.exists(args.causaldb_raw):
        bench_multimethod(enc, pip_head, args.causaldb_raw, test_ids)
    bench_filtering(pred, true, pval, z, beta)
    bench_calibration(enc, pip_head, cdb, test_ids)
    bench_speed(enc, pip_head, args.sentinel_dir)
    if os.path.exists(args.causaldb_raw):
        bench_functional(enc, pip_head, args.causaldb_raw, test_ids)

    log("\n" + "=" * 70)
    log("ALL BENCHMARKS COMPLETE")
    log("=" * 70)
