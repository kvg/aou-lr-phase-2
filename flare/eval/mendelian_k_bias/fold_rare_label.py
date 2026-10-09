#!/usr/bin/env python3
"""Cost of keeping a rare ancestry label: one VCF scored with and without folding it away.

The uniform-label simulation (``simulate_k_bias.py``) gives +9-11% per extra label.
A real five-panel recipe is not uniform: SAS is ~5% of ancestry. This scores the SAME
synthetic pin-like VCF twice with the real scorer, as-is (5 labels) and with SAS folded
into EUR (4 labels), at three per-call error rates and three seeds each, and prints how
much higher the five-label rate is. It sets ``LABEL_COUNT_BIAS_PER_LABEL`` in
``scripts/flare_lai_exp.py``.

    python3 flare/eval/mendelian_k_bias/fold_rare_label.py
"""
from __future__ import annotations

import statistics as st
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
import test_flare_alphabet_pipeline as T  # noqa: E402


def main() -> int:
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        T._shared_files(tmp)
        to_pin = lambda p, rng: T.PIN_CODE[p]  # noqa: E731
        for eps in (0.005, 0.02, 0.05):
            raw, folded = [], []
            for seed in (11, 12, 13):
                gz = T._write(tmp, f"pin_{eps}_{seed}", to_pin, eps, range(5), T.PIN_HDR, seed)
                raw.append(T._score(tmp, f"raw_{eps}_{seed}", gz)["violations_per_informative_locus"])
                folded.append(T._score(
                    tmp, f"fold_{eps}_{seed}", gz, "--project-labels",
                    str(tmp / "pin_sas_eur.labels.tsv"))["violations_per_informative_locus"])
            r, f = st.mean(raw), st.mean(folded)
            print(f"eps={eps:<6} five labels={r:.4f}  SAS->EUR (four)={f:.4f}  "
                  f"five-label rate is +{(r / f - 1) * 100:.1f}% higher")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
