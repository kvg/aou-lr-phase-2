#!/usr/bin/env python3
"""Parse and write FLARE ``.model`` files.

FLARE model layout (browning-lab/flare README "Model file format").
With ``A`` ancestries and ``P`` reference panels the file has ``2A + 5``
data lines (comment / blank lines allowed anywhere):

1. Ancestry names (length A); first ancestry has index 0.
2. Reference panel names (length P); first panel has index 0.
3. Generations since admixture ``T`` (scalar).
4. Admixture proportions (length A); must sum to ~1.0.
5. Next A lines: panel-weight matrix (A x P). Entry (i, j) is the
   probability a model-state haplotype is in panel j given ancestry i.
6. Next A lines: miscopy-rate matrix (A x P). Entry (i, j) is the
   probability a model-state haplotype and the admixed haplotype carry
   different alleles given panel j and ancestry i. This is FLARE's
   true ``mu`` block — **not** the proportion vector on line 4.
7. Final data line: pre-admixture IBD-rate vector (length A).

Historical bug: ``flare_summarize_models`` labeled line 4 as ``mu_*``.
Those columns are admixture proportions; use ``prop_*`` / ``Model.props``.
"""

from __future__ import annotations

import argparse
import copy
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence, TextIO

PROP_SUM_TOL = 1e-4


@dataclass
class Model:
    """In-memory FLARE model (all named blocks)."""

    ancestries: list[str]
    panels: list[str]
    t_gen: float
    props: list[float]
    panel_weights: list[list[float]]
    miscopy: list[list[float]]
    ibd_rates: list[float]
    # Original raw text for byte-identical round-trips when unchanged.
    _raw: Optional[str] = field(default=None, repr=False, compare=False)

    @property
    def n_ancestries(self) -> int:
        return len(self.ancestries)

    @property
    def n_panels(self) -> int:
        return len(self.panels)

    def validate(self) -> None:
        a, p = self.n_ancestries, self.n_panels
        if a < 1 or p < 1:
            raise ValueError("model needs >=1 ancestry and >=1 panel")
        if len(self.props) != a:
            raise ValueError(f"props length {len(self.props)} != A={a}")
        if abs(sum(self.props) - 1.0) > PROP_SUM_TOL:
            raise ValueError(
                f"props sum to {sum(self.props)!r}, expected 1.0 ± {PROP_SUM_TOL}"
            )
        if len(self.panel_weights) != a:
            raise ValueError(f"panel_weights rows {len(self.panel_weights)} != A={a}")
        if len(self.miscopy) != a:
            raise ValueError(f"miscopy rows {len(self.miscopy)} != A={a}")
        for i, row in enumerate(self.panel_weights):
            if len(row) != p:
                raise ValueError(f"panel_weights[{i}] length {len(row)} != P={p}")
        for i, row in enumerate(self.miscopy):
            if len(row) != p:
                raise ValueError(f"miscopy[{i}] length {len(row)} != P={p}")
        if len(self.ibd_rates) != a:
            raise ValueError(f"ibd_rates length {len(self.ibd_rates)} != A={a}")

    def mean_miscopy(self) -> list[float]:
        """Per-ancestry mean miscopy rate across panels."""
        return [sum(row) / len(row) for row in self.miscopy]


def _is_data_line(text: str) -> bool:
    s = text.strip()
    return bool(s) and not s.startswith("#")


def parse_model(path: Path | str) -> Model:
    """Parse a FLARE ``.model`` file into a :class:`Model`."""
    path = Path(path)
    raw = path.read_text()
    data: list[str] = []
    for line in raw.splitlines():
        if _is_data_line(line):
            data.append(line.strip())
    if len(data) < 5:
        raise ValueError(f"{path}: expected >=5 data lines, found {len(data)}")

    ancestries = data[0].split()
    panels = data[1].split()
    a, p = len(ancestries), len(panels)
    expect = 2 * a + 5
    if len(data) != expect:
        raise ValueError(
            f"{path}: expected {expect} data lines (2A+5 with A={a}), found {len(data)}"
        )
    try:
        t_gen = float(data[2].split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{path}: bad T line {data[2]!r}") from exc
    props = [float(x) for x in data[3].split()]
    panel_weights = [[float(x) for x in data[4 + i].split()] for i in range(a)]
    miscopy = [[float(x) for x in data[4 + a + i].split()] for i in range(a)]
    ibd_rates = [float(x) for x in data[4 + 2 * a].split()]

    model = Model(
        ancestries=ancestries,
        panels=panels,
        t_gen=t_gen,
        props=props,
        panel_weights=panel_weights,
        miscopy=miscopy,
        ibd_rates=ibd_rates,
        _raw=raw if raw.endswith("\n") else raw + "\n",
    )
    model.validate()
    return model


def _fmt_row(xs: Sequence[float]) -> str:
    return "\t".join(f"{x:.8g}" for x in xs)


def write_model(model: Model, path: Path | str, *, preserve_raw: bool = True) -> None:
    """Write ``model`` to ``path``.

    If ``preserve_raw`` and the model still carries unchanged original bytes,
    write those bytes (byte-identical round-trip). Otherwise emit a canonical
    whitespace-delimited file with a trailing newline.
    """
    path = Path(path)
    model.validate()
    if preserve_raw and model._raw is not None:
        path.write_text(model._raw)
        return
    lines = [
        "\t".join(model.ancestries),
        "\t".join(model.panels),
        f"{model.t_gen:g}",
        _fmt_row(model.props),
    ]
    for row in model.panel_weights:
        lines.append(_fmt_row(row))
    for row in model.miscopy:
        lines.append(_fmt_row(row))
    lines.append(_fmt_row(model.ibd_rates))
    path.write_text("\n".join(lines) + "\n")


def pop_of_model_path(path: Path | str) -> str:
    """Extract population from ``*.<POP>.model`` basename."""
    name = Path(path).name
    if name.endswith(".model"):
        name = name[: -len(".model")]
    if "." not in name:
        raise ValueError(
            f"model basename {Path(path).name!r} must look like prefix.<POP>.model"
        )
    return name.rsplit(".", 1)[-1]


# Back-compat alias used by summarize / tests.
pop_from_model_path = pop_of_model_path


def rewrite_t(model: Model, t_gen: float) -> Model:
    out = copy.deepcopy(model)
    out.t_gen = float(t_gen)
    out._raw = None
    out.validate()
    return out


def set_props(model: Model, props_by_anc: dict[str, float], *, min_sum: float = 0.5) -> Model:
    """Overwrite proportions from ``{ancestry: weight}``; normalize to sum 1.

    Ancestry keys are matched case-insensitively to model labels.
    """
    out = copy.deepcopy(model)
    lower_to_canon = {a.lower(): a for a in out.ancestries}
    normalized: dict[str, float] = {}
    unknown: list[str] = []
    for key, val in props_by_anc.items():
        canon = lower_to_canon.get(key.lower())
        if canon is None:
            unknown.append(key)
        else:
            normalized[canon] = float(val)
    if unknown:
        raise ValueError(
            f"props name ancestries not in model: {sorted(unknown)}; "
            f"model has {out.ancestries}"
        )
    raw = [float(normalized.get(a, 0.0)) for a in out.ancestries]
    s = sum(raw)
    if s < min_sum:
        raise ValueError(
            f"props sum to {s} before normalization (< {min_sum}); likely a typo"
        )
    out.props = [x / s for x in raw]
    out._raw = None
    out.validate()
    return out


def _identity_weights(a: int, p: int, ancestries: list[str], panels: list[str]) -> list[list[float]]:
    """Diagonal-ish defaults: match ancestry name to panel name when possible."""
    panel_idx = {name: j for j, name in enumerate(panels)}
    rows: list[list[float]] = []
    for anc in ancestries:
        row = [0.0] * p
        if anc in panel_idx:
            row[panel_idx[anc]] = 1.0
        else:
            # Uniform fallback if no same-named panel.
            row = [1.0 / p] * p
        rows.append(row)
    return rows


def reset_panel_weights_default(model: Model) -> Model:
    out = copy.deepcopy(model)
    out.panel_weights = _identity_weights(
        out.n_ancestries, out.n_panels, out.ancestries, out.panels
    )
    out._raw = None
    out.validate()
    return out


def reset_miscopy_default(model: Model, rate: float = 0.0) -> Model:
    out = copy.deepcopy(model)
    out.miscopy = [
        [float(rate)] * out.n_panels for _ in range(out.n_ancestries)
    ]
    out._raw = None
    out.validate()
    return out


def parse_props_by_pop(raw: str) -> dict[str, dict[str, float]]:
    """Parse ``AFR:AFR=0.80,EUR=0.20;EAS:EAS=0.97`` into ``{pop: {anc: w}}``."""
    out: dict[str, dict[str, float]] = {}
    text = (raw or "").strip()
    if not text:
        return out
    for chunk in text.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if ":" not in chunk:
            raise ValueError(
                f"props_by_pop entry {chunk!r} must look like POP:ANC=w,ANC=w"
            )
        pop, rest = chunk.split(":", 1)
        pop = pop.strip()
        weights: dict[str, float] = {}
        for part in rest.split(","):
            part = part.strip()
            if not part:
                continue
            if "=" not in part:
                raise ValueError(f"props_by_pop weight {part!r} must look like ANC=w")
            anc, val = part.split("=", 1)
            weights[anc.strip()] = float(val.strip())
        out[pop] = weights
    return out


def parse_gen_by_pop(raw: str) -> dict[str, float]:
    overrides: dict[str, float] = {}
    text = (raw or "").strip()
    if not text:
        return overrides
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise ValueError(f"gen_by_pop entry {part!r} must look like POP:T")
        key, val = part.split(":", 1)
        overrides[key.strip()] = float(val.strip())
    return overrides


def build_input_model(
    *,
    pop: str,
    gen: float,
    gen_by_pop: str,
    pop_model_path: Optional[Path],
    template_model_path: Optional[Path],
    props_by_pop: str,
    inherit_props: bool,
    inherit_panel_weights: bool,
    inherit_mu: bool,
) -> tuple[Optional[Model], dict[str, str]]:
    """Assemble the model FLARE should be given for one shard.

    Returns ``(model_or_None, sources)`` where sources maps block name to
    ``default`` / ``template`` / ``pop_model`` / ``override``. ``None`` means
    omit ``model=`` and let FLARE use built-in defaults.
    """
    sources = {
        "t_gen": "override",
        "props": "default",
        "panel_weights": "default",
        "miscopy": "default",
        "ibd_rates": "default",
    }
    overrides = parse_gen_by_pop(gen_by_pop)
    t_gen = overrides.get(pop, gen)
    props_map = parse_props_by_pop(props_by_pop)
    want_props_override = pop in props_map
    need_source = (
        want_props_override
        or not inherit_props
        or not inherit_panel_weights
        or not inherit_mu
        or pop_model_path is not None
        or template_model_path is not None
    )

    source_path: Optional[Path] = None
    source_kind = "default"
    if pop_model_path is not None:
        source_path = Path(pop_model_path)
        source_kind = "pop_model"
    elif template_model_path is not None:
        source_path = Path(template_model_path)
        source_kind = "template"

    if source_path is None:
        if need_source and (
            want_props_override
            or not inherit_props
            or not inherit_panel_weights
            or not inherit_mu
        ):
            raise ValueError(
                f"{pop}: props/inherit overrides require pop_models or template_model"
            )
        # No model file: FLARE defaults; T still passed via gen=.
        return None, {**sources, "t_gen": "gen_arg"}

    model = parse_model(source_path)
    for block in ("props", "panel_weights", "miscopy", "ibd_rates"):
        sources[block] = source_kind

    model = rewrite_t(model, t_gen)

    if want_props_override:
        model = set_props(model, props_map[pop])
        sources["props"] = "override"
    elif not inherit_props:
        # Uniform proportions.
        n = model.n_ancestries
        model = set_props(model, {a: 1.0 / n for a in model.ancestries}, min_sum=0.5)
        sources["props"] = "default"

    if not inherit_panel_weights:
        model = reset_panel_weights_default(model)
        sources["panel_weights"] = "default"

    if not inherit_mu:
        model = reset_miscopy_default(model)
        sources["miscopy"] = "default"

    return model, sources


def match_pop_model(
    pop: str,
    pop_models: Sequence[Path | str],
    *,
    require: bool,
) -> Optional[Path]:
    """Return the model path for ``pop``, or fail if ``require`` and missing."""
    paths = [Path(p) for p in pop_models if str(p).strip()]
    if not paths:
        if require:
            raise ValueError(f"{pop}: pop_models is non-empty requirement but list empty")
        return None
    parsed = []
    matches = []
    for path in paths:
        try:
            got = pop_of_model_path(path)
        except ValueError as exc:
            raise ValueError(f"unparseable pop_models entry {path}: {exc}") from exc
        parsed.append((got, path))
        if got == pop:
            matches.append(path)
    if len(matches) > 1:
        raise ValueError(f"{pop}: multiple model files: {matches}")
    if not matches:
        if require or paths:
            detail = ", ".join(f"{g}<-{p.name}" for g, p in parsed)
            raise ValueError(
                f"{pop}: no matching entry in pop_models (parsed: {detail})"
            )
        return None
    return matches[0]


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    show = sub.add_parser("show", help="Print parsed blocks")
    show.add_argument("model", type=Path)

    rt = sub.add_parser("roundtrip", help="parse+write and compare bytes")
    rt.add_argument("model", type=Path)
    rt.add_argument("--out", type=Path, required=True)

    build = sub.add_parser("build-input", help="Build rewritten input model for a shard")
    build.add_argument("--pop", required=True)
    build.add_argument("--gen", type=float, required=True)
    build.add_argument("--gen-by-pop", default="")
    build.add_argument("--pop-model", type=Path)
    build.add_argument("--template-model", type=Path)
    build.add_argument("--props-by-pop", default="")
    build.add_argument("--inherit-props", action=argparse.BooleanOptionalAction, default=True)
    build.add_argument(
        "--inherit-panel-weights", action=argparse.BooleanOptionalAction, default=True
    )
    build.add_argument("--inherit-mu", action=argparse.BooleanOptionalAction, default=True)
    build.add_argument("--out", type=Path, required=True)

    args = p.parse_args(argv)
    if args.cmd == "show":
        m = parse_model(args.model)
        print(f"ancestries={m.ancestries}")
        print(f"panels={m.panels}")
        print(f"t_gen={m.t_gen}")
        print(f"props={m.props} sum={sum(m.props)}")
        print(f"mean_miscopy={m.mean_miscopy()}")
        print(f"ibd_rates={m.ibd_rates}")
        return 0
    if args.cmd == "roundtrip":
        m = parse_model(args.model)
        write_model(m, args.out, preserve_raw=True)
        a = Path(args.model).read_bytes()
        b = Path(args.out).read_bytes()
        if a != b:
            raise SystemExit("round-trip byte mismatch")
        print("ok", args.out)
        return 0
    if args.cmd == "build-input":
        model, sources = build_input_model(
            pop=args.pop,
            gen=args.gen,
            gen_by_pop=args.gen_by_pop,
            pop_model_path=args.pop_model,
            template_model_path=args.template_model,
            props_by_pop=args.props_by_pop,
            inherit_props=args.inherit_props,
            inherit_panel_weights=args.inherit_panel_weights,
            inherit_mu=args.inherit_mu,
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        if model is None:
            args.out.write_text("")  # sentinel: no model=
            print("sources", sources, file=sys.stderr)
            print("NO_MODEL", file=sys.stderr)
            return 0
        write_model(model, args.out, preserve_raw=False)
        print("sources", sources, file=sys.stderr)
        print(args.out.read_text(), file=sys.stderr)
        return 0
    raise SystemExit(f"unknown cmd {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
