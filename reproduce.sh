#!/bin/bash
# ============================================================================
# G-Atlas: EVERYTHING. Download → preprocess → train → benchmark.
# One script. No flags. Sentinel GWAS require OpenGWAS API access.
# ============================================================================

set -e
cd "$(dirname "$0")"
DEVICE="${DEVICE:-cuda}"
WORKERS="${WORKERS:-4}"
STEPS="${STEPS:-100000}"  # 100K steps to match original 2-stage training

echo "============================================"
echo "G-Atlas — full pipeline"
echo "============================================"

mkdir -p models results data/raw_gwas data/gwas_31k data/causaldb

# ── REFERENCES (HM3, sentinel GWAS) ──
python3 -u data/download_references.py --output_dir data/

# ── GWAS DATA ──
NPZ_COUNT=$(ls data/gwas_31k/*.npz 2>/dev/null | wc -l)
if [ "$NPZ_COUNT" -lt 1000 ]; then
    RAW_COUNT=$(ls data/raw_gwas/*.tsv.gz 2>/dev/null | wc -l)
    if [ "$RAW_COUNT" -lt 1000 ]; then
        echo "[0a] Downloading raw GWAS from EBI..."
        python3 -u data/download_ebi.py --output_dir data/raw_gwas
    fi
    echo "[0b] Preprocessing raw GWAS → .npz blocks..."
    python3 -u data/preprocess.py --input_dir data/raw_gwas --output_dir data/gwas_31k --workers "$WORKERS"
fi
echo "GWAS .npz: $(ls data/gwas_31k/*.npz | wc -l) files ✓"

# ── CAUSALDB ──
if [ ! -f data/causaldb/susie_pips.parquet ]; then
    if [ ! -f data/causaldb/credible_set.txt ] && [ ! -f data/causaldb/v2.1/credible_set.txt ]; then
        echo "[0c] Downloading CAUSALdb2..."
        python3 -u data/download_causaldb.py --output_dir data/causaldb
    fi
    RAW_TXT="data/causaldb/credible_set.txt"
    [ -f data/causaldb/v2.1/credible_set.txt ] && RAW_TXT="data/causaldb/v2.1/credible_set.txt"
    echo "[0c] Filtering CAUSALdb2 to HM3 → parquet..."
    python3 -u data/download_causaldb.py --output_dir data/causaldb --credible_set "$RAW_TXT" --hm3_bim data/hm3.bim
fi
echo "CAUSALdb2: ✓"

CAUSALDB_RAW="data/causaldb/v2.1/credible_set.txt"
if [ ! -f "$CAUSALDB_RAW" ]; then
    CAUSALDB_RAW="data/causaldb/credible_set.txt"
fi
if [ ! -f "$CAUSALDB_RAW" ]; then
    echo "ERROR: CAUSALdb2 credible-set table missing. Expected data/causaldb/v2.1/credible_set.txt or data/causaldb/credible_set.txt." >&2
    exit 1
fi

# ── TRAIN ENCODER ──
echo ""
ENCODER_RETRAINED=0
if [ -f models/encoder.pt ]; then
    echo "[1/3] encoder.pt exists — skipping pretrain. Delete to retrain."
else
    echo "[1/3] Pretraining encoder ($STEPS steps)..."
    python3 -u train/pretrain.py \
        --gwas_dir data/gwas_31k \
        --output models/encoder.pt \
        --device "$DEVICE" \
        --steps "$STEPS" \
        2>&1 | tee results/pretrain.log
    ENCODER_RETRAINED=1
fi

# ── TRAIN PIP HEAD ──
echo ""
if [ "$ENCODER_RETRAINED" -eq 0 ] && [ -f models/pip_head.pt ] && [ ! models/encoder.pt -nt models/pip_head.pt ]; then
    echo "[2/3] pip_head.pt exists — skipping PIP training. Delete to retrain."
else
    echo "[2/3] Training PIP head..."
    python3 -u train/train_pip.py \
        --encoder models/encoder.pt \
        --causaldb data/causaldb/susie_pips.parquet \
        --output models/pip_head.pt \
        --device "$DEVICE" \
        2>&1 | tee results/train_pip.log
fi

# ── BENCHMARK ──
echo ""
echo "[3/3] Benchmarks..."
python3 -u eval/benchmark.py \
    --causaldb data/causaldb/susie_pips.parquet \
    --causaldb_raw "$CAUSALDB_RAW" \
    --device "$DEVICE" \
    2>&1 | tee results/benchmark.log

echo ""
echo "============================================"
echo "Done. Results: results/benchmark.log"
echo "Inference: python3 inference/predict.py INPUT.tsv OUTPUT.scores.tsv"
echo "============================================"
