#!/usr/bin/env python3
"""Patch SAIGE v1.3.3 for oneTBB: tbb::concurrent_vector needs an explicit include."""

from pathlib import Path

p = Path("/opt/SAIGE/src/SAIGE_fitGLMM_fast.cpp")
text = p.read_text()
if "concurrent_vector.h" in text:
    print("SAIGE already includes concurrent_vector header")
    raise SystemExit(0)

needle = "#include <RcppParallel.h>"
insert = needle + "\n#include <oneapi/tbb/concurrent_vector.h>"
if needle not in text:
    raise SystemExit("RcppParallel include not found for SAIGE TBB patch")
p.write_text(text.replace(needle, insert, 1))
print("Patched SAIGE_fitGLMM_fast.cpp for oneTBB concurrent_vector")
