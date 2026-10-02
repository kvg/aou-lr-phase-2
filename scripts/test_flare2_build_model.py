#!/usr/bin/env python3
"""Unit tests for scripts/flare2_build_model.py.

Uses a stub in place of upstream ``create_model_file.py`` that writes the same
model layout (GMM-ordered ``anc_i`` labels + correlation comment lines), so the
relabel / reorder / gate logic runs without scikit-learn. Set
``FLARE2_CREATE_MODEL_SCRIPT`` to the real upstream script to also run an
end-to-end clustering check on a synthetic ``.panels`` file.
"""

from __future__ import annotations

import json
import os
import random
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from flare2_build_model import main  # noqa: E402
from flare_model import parse_model  # noqa: E402

# Upstream-shaped output: anc_0 is EUR-like, anc_1 AFR-like, anc_2 EAS-like.
STUB = '''
import sys
nanc, panels, prefix = int(sys.argv[1]), sys.argv[2], sys.argv[3]
p_rows = [
    [0.05, 0.85, 0.05, 0.05],
    [0.90, 0.04, 0.03, 0.03],
    [0.02, 0.08, 0.80, 0.10],
]
with open(prefix + ".model", "w") as fh:
    fh.write("# model is output from cluster_model_singlewin.py\\n")
    fh.write("# lag-1-window correlations are {}\\n")
    fh.write("# max cross-correlation is 0.01\\n")
    fh.write("# min auto-correlation is MINAC\\n")
    fh.write("anc_0\\tanc_1\\tanc_2\\n\\n")
    fh.write("afr\\teur\\teas\\tsas\\n\\n")
    fh.write("10\\n\\n")
    fh.write("0.2\\t0.5\\t0.3\\n\\n")
    for row in p_rows:
        fh.write("\\t".join(str(x) for x in row) + "\\n")
    fh.write("\\n")
    for i in range(3):
        fh.write("\\t".join([str(0.001 * (i + 1))] * 4) + "\\n")
    fh.write("\\n0.1\\t0.2\\t0.3\\n")
open(prefix + ".log", "w").write("stub\\n")
'''


def _write_stub(td: Path, min_ac: float) -> Path:
    stub = td / "create_model_file.py"
    stub.write_text(STUB.replace("MINAC", str(min_ac)))
    return stub


def test_relabel_and_reorder():
    with tempfile.TemporaryDirectory() as tdn:
        td = Path(tdn)
        stub = _write_stub(td, 0.6)
        (td / "x.panels").write_text("")
        prefix = str(td / "out")
        assert main([
            "--create-model-script", str(stub), "--panels", str(td / "x.panels"),
            "--nanc", "3", "--out-prefix", prefix,
        ]) == 0
        m = parse_model(prefix + ".model")
        # Canonical order: sorted by dominant panel name.
        assert m.ancestries == ["anc0_afr", "anc1_eas", "anc2_eur"], m.ancestries
        # Blocks permuted with the ancestries (old anc_1, anc_2, anc_0).
        assert m.props == [0.5, 0.3, 0.2]
        assert m.panel_weights[0][0] == 0.9
        assert m.panel_weights[2][1] == 0.85
        assert abs(m.miscopy[0][0] - 0.002) < 1e-12
        assert m.ibd_rates == [0.2, 0.3, 0.1]
        summary = json.loads(Path(prefix + ".summary.json").read_text())
        assert summary["passed"] is True
        assert summary["min_autocorr"] == 0.6
        assert [r["upstream_label"] for r in summary["ancestries"]] == ["anc_1", "anc_2", "anc_0"]
        labels = Path(prefix + ".labels.tsv").read_text().splitlines()
        assert labels[1].startswith("0\tanc0_afr\tanc_1\tafr")


def test_autocorr_gate_fails():
    with tempfile.TemporaryDirectory() as tdn:
        td = Path(tdn)
        stub = _write_stub(td, 0.1)
        (td / "x.panels").write_text("")
        prefix = str(td / "out")
        try:
            main([
                "--create-model-script", str(stub), "--panels", str(td / "x.panels"),
                "--nanc", "3", "--min-autocorr", "0.25", "--out-prefix", prefix,
            ])
        except SystemExit as exc:
            assert "autocorrelation" in str(exc)
        else:
            raise AssertionError("expected gate failure")
        summary = json.loads(Path(prefix + ".summary.json").read_text())
        assert summary["passed"] is False


def _synthetic_panels(path: Path, *, n_haps: int = 60, n_loc: int = 40, seed: int = 3) -> None:
    """Three ancestries with long tracts; copying probs near one-hot + noise."""
    rng = random.Random(seed)
    panels = ["afr", "eur", "eas"]
    lines = ["#nrefhaps=300", "#loc\tmarker\thap\trho\t" + "\t".join(panels)]
    labels = [[rng.randrange(3)] for _ in range(n_haps)]
    for h in range(n_haps):
        for _ in range(1, n_loc):
            prev = labels[h][-1]
            labels[h].append(prev if rng.random() > 0.05 else rng.randrange(3))
    for loc in range(n_loc):
        for h in range(n_haps):
            probs = [rng.uniform(0.0, 0.08) for _ in panels]
            probs[labels[h][loc]] += 0.9
            s = sum(probs)
            probs = [x / s for x in probs]
            lines.append(
                f"{loc}\tchr20:{1000 * (loc + 1)}\t{h}\t0.5\t" + "\t".join(f"{x:.5f}" for x in probs)
            )
    path.write_text("\n".join(lines) + "\n")


def test_upstream_end_to_end():
    script = os.environ.get("FLARE2_CREATE_MODEL_SCRIPT", "")
    if not script:
        print("skip test_upstream_end_to_end (FLARE2_CREATE_MODEL_SCRIPT unset)")
        return
    with tempfile.TemporaryDirectory() as tdn:
        td = Path(tdn)
        _synthetic_panels(td / "syn.panels")
        prefix = str(td / "syn")
        outs = []
        for _ in range(2):
            assert main([
                "--create-model-script", script, "--panels", str(td / "syn.panels"),
                "--nanc", "3", "--seed", "7", "--out-prefix", prefix,
            ]) == 0
            outs.append(Path(prefix + ".model").read_text())
        assert outs[0] == outs[1], "seeded runs should be identical"
        m = parse_model(prefix + ".model")
        assert m.ancestries == ["anc0_afr", "anc1_eas", "anc2_eur"], m.ancestries


if __name__ == "__main__":
    test_relabel_and_reorder()
    test_autocorr_gate_fails()
    test_upstream_end_to_end()
    print("ok")
