#!/usr/bin/env python3
"""CLI: python -m manuscript_figures  [fig1 ...]"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow `python manuscript_figures` and `python -m manuscript_figures`.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from manuscript_figures import FIGURES, render_all


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("figures", nargs="*", help="Optional subset, e.g. fig1 fig3")
    args = p.parse_args()
    which = args.figures or None
    if which:
        known = {n for n, _ in FIGURES}
        unknown = [w for w in which if not any(w in n for n in known)]
        if unknown:
            p.error(f"unknown figure(s) {unknown}; known: {sorted(known)}")
    results = render_all(which)
    if not results:
        p.error("no figures rendered")
    print(f"rendered {len(results)} figure(s)")


if __name__ == "__main__":
    main()
