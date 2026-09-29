"""Sample peak RSS / CPU for a process tree. Stdlib only (prepended into Batch -c)."""

from __future__ import annotations

import os
import threading
import time
from pathlib import Path


def _proc_pids(root: int) -> list[int]:
    by_ppid: dict[int, list[int]] = {}
    try:
        for entry in Path("/proc").iterdir():
            if not entry.name.isdigit():
                continue
            pid = int(entry.name)
            try:
                stat = (entry / "stat").read_text()
            except OSError:
                continue
            rparen = stat.rfind(")")
            if rparen < 0:
                continue
            rest = stat[rparen + 2 :].split()
            if len(rest) < 2:
                continue
            ppid = int(rest[1])
            by_ppid.setdefault(ppid, []).append(pid)
    except OSError:
        return [root]
    seen = {root}
    stack = [root]
    while stack:
        cur = stack.pop()
        for kid in by_ppid.get(cur, []):
            if kid not in seen:
                seen.add(kid)
                stack.append(kid)
    return list(seen)


def _rss_bytes(pid: int) -> int:
    try:
        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) * 1024
    except (OSError, ValueError, IndexError):
        return 0
    return 0


def _cpu_ticks(pid: int) -> int:
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
    except OSError:
        return 0
    rparen = stat.rfind(")")
    rest = stat[rparen + 2 :].split()
    if len(rest) < 13:
        return 0
    try:
        return int(rest[11]) + int(rest[12])
    except ValueError:
        return 0


def _dir_bytes(path: Path) -> int:
    total = 0
    try:
        for p in path.rglob("*"):
            if p.is_file() and not p.is_symlink():
                try:
                    total += p.stat().st_size
                except OSError:
                    pass
    except OSError:
        return 0
    return total


class ResourceMonitor:
    def __init__(self, work: Path, interval: float = 1.0):
        self.work = work
        self.interval = interval
        self.peak_rss_bytes = 0
        self.peak_disk_bytes = 0
        self.peak_nproc = 1
        self.n_samples = 0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._t0 = 0.0
        self._cpu0 = 0
        self._root = os.getpid()

    def start(self) -> None:
        self._root = os.getpid()
        self._t0 = time.time()
        self._cpu0 = sum(_cpu_ticks(p) for p in _proc_pids(self._root))
        self._sample()
        self._thread = threading.Thread(target=self._loop, name="resource-monitor", daemon=True)
        self._thread.start()

    def _sample(self) -> None:
        pids = _proc_pids(self._root)
        rss = sum(_rss_bytes(p) for p in pids)
        self.peak_rss_bytes = max(self.peak_rss_bytes, rss)
        self.peak_nproc = max(self.peak_nproc, len(pids))
        self.peak_disk_bytes = max(self.peak_disk_bytes, _dir_bytes(self.work))
        self.n_samples += 1

    def _loop(self) -> None:
        while not self._stop.wait(self.interval):
            try:
                self._sample()
            except Exception:
                return

    def stop(self) -> dict[str, float | int]:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=3)
        try:
            self._sample()
        except Exception:
            pass
        wall = max(time.time() - self._t0, 1e-6)
        ticks = sum(_cpu_ticks(p) for p in _proc_pids(self._root)) - self._cpu0
        try:
            clk = float(os.sysconf("SC_CLK_TCK") or 100)
        except (OSError, ValueError):
            clk = 100.0
        cpu_sec = max(ticks, 0) / clk
        return {
            "peak_rss_bytes": int(self.peak_rss_bytes),
            "peak_disk_bytes": int(self.peak_disk_bytes),
            "wall_sec": float(wall),
            "cpu_sec": float(cpu_sec),
            "cpu_cores_avg": float(cpu_sec / wall),
            "peak_nproc": int(self.peak_nproc),
            "n_samples": int(self.n_samples),
            "n_cpus_os": int(os.cpu_count() or 0),
        }


RESOURCE_FIELDS = (
    "sample_id",
    "stage",
    "peak_rss_bytes",
    "peak_disk_bytes",
    "wall_sec",
    "cpu_sec",
    "cpu_cores_avg",
    "peak_nproc",
    "n_samples",
    "n_cpus_os",
)


def write_resource_tsv(path: Path, row: dict) -> None:
    path.write_text(
        "\t".join(RESOURCE_FIELDS)
        + "\n"
        + "\t".join(str(row.get(k, "")) for k in RESOURCE_FIELDS)
        + "\n",
        encoding="utf-8",
    )


def resource_monitor_cls():
    cls = globals().get("ResourceMonitor")
    if cls is not None:
        return cls
    try:
        from resource_monitor import ResourceMonitor as imported
    except ImportError:
        return None
    return imported
