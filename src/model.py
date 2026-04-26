"""
G-Atlas model architecture.
Pretrained GWAS encoder + fine-mapping PIP head.

Encoder: 4-layer transformer, 128d, 4 heads, 811K params.
  - Pretrained on 31,200 EBI GWAS via masked z-score reconstruction + denoising.
  - Input: (batch, seq, 6) variant features [beta, SE, -log10p, MAF, rel_pos, chrom/23]
  - Output: (batch, seq, 128) contextual variant representations.

PIP Head: Linear(128,64) → GELU → Dropout → Linear(64,1) → Sigmoid.
  - Trained on 8,232 CAUSALdb2 GWAS with real SuSiE PIPs.
  - Uses LD adapter to incorporate local LD context.
  - Output: per-variant causal prioritization score (trained to approximate SuSiE posteriors).
"""

import torch
import torch.nn as nn


class GAtlasEncoder(nn.Module):
    """Pretrained GWAS block encoder."""

    def __init__(self, n_input=6, d_model=128, n_heads=4, n_layers=4, dropout=0.1):
        super().__init__()
        self.d_model = d_model
        self.input_proj = nn.Sequential(
            nn.Linear(n_input, d_model),
            nn.GELU(),
            nn.LayerNorm(d_model),
            nn.Linear(d_model, d_model),
        )
        self.pos_proj = nn.Linear(1, d_model)
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model, nhead=n_heads,
            dim_feedforward=d_model * 4,
            dropout=dropout, batch_first=True,
            activation='gelu',
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)

    def forward(self, features, padding_mask=None):
        """
        features: (batch, seq, 6)
        padding_mask: (batch, seq) bool, True = padded
        Returns: (batch, seq, 128)
        """
        x = self.input_proj(features)
        x = x + self.pos_proj(features[:, :, 4:5])
        if padding_mask is not None:
            x = self.transformer(x, src_key_padding_mask=padding_mask)
        else:
            x = self.transformer(x)
        return x


class LDAdapter(nn.Module):
    """Adapts LD-aware features into encoder space with learnable gate."""

    def __init__(self, n_ld=10, d_model=128, dropout=0.1):
        super().__init__()
        self.proj = nn.Sequential(
            nn.Linear(n_ld, d_model // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model // 2, d_model),
        )
        self.gate = nn.Parameter(torch.zeros(1))

    def forward(self, ld_features):
        return self.proj(ld_features) * torch.sigmoid(self.gate)


class PIPHead(nn.Module):
    """Predicts per-variant causal prioritization score (approximates SuSiE PIP)."""

    def __init__(self, d_model=128, n_ld=10):
        super().__init__()
        self.ld_adapter = LDAdapter(n_ld, d_model, dropout=0.1)
        self.head = nn.Sequential(
            nn.Linear(d_model, 64),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(64, 1),
        )

    def forward(self, h_encoder, ld_features, padding_mask=None):
        """
        h_encoder: (batch, seq, 128) from GAtlasEncoder
        ld_features: (batch, seq, 10) local LD context
        Returns: (batch, seq) logits. Apply sigmoid for PIP.
        """
        h = h_encoder + self.ld_adapter(ld_features)
        return self.head(h).squeeze(-1)
