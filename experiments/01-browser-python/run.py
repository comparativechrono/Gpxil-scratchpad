"""Runs gpxiltools in the browser's Python on real gpxil files, three ways, and measures time and memory.

    python run.py --data DIR [--files chr21.gpxil.tif,chr1.gpxil.tif,...]

DIR holds the gpxil files (served at /data/ for the web runs, and given to the page as chosen files for the
others). Needs Playwright with Chromium, and the wheels in wheels/ (see README.md). Writes results.json.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
import urllib.request

import psutil
from playwright.sync_api import sync_playwright

HERE = os.path.dirname(os.path.abspath(__file__))
MB = 1 << 20

# Python run in the browser. P is the file's path there.
OPEN = """
import gpxiltools
g = gpxiltools.GpxilFile(P)
c = g[g.contigs[0]]
{"contigs": len(g.contigs), "first": c.name, "length": c.length, "planes": list(c.pages)}
"""
REGION_10K = """
g = gpxiltools.GpxilFile(P)
c = g[g.contigs[0]]
a = c.length // 2
s = c.sequence(a, a + 10_000)
{"bases": len(s), "gc": round(sum(b in "GCgc" for b in s) / len(s), 3)}
"""
REGION_1M = """
g = gpxiltools.GpxilFile(P)
c = g[g.contigs[0]]
a = c.length // 3
v = c.read("ref", a, a + 1_000_000)
{"bases": int(v.size)}
"""
CLI_INFO = """
import io, contextlib
from gpxiltools.cli import main
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    rc = main(["info", P])
{"rc": rc, "lines": len(buf.getvalue().splitlines()), "head": buf.getvalue().splitlines()[:3]}
"""
RENDER = """
import os
from gpxiltools.cli import main
g = gpxiltools.GpxilFile(P)
c = g[g.contigs[0]]
a = c.length // 2
mode = "genotype" if "gt" in c.pages else "reference"
rc = main(["render", P, f"{c.name}:{a + 1}-{a + 200_000}", "-o", "/tmp/region.png", "--mode", mode])
{"rc": rc, "png": os.path.getsize("/tmp/region.png")}
"""
FASTA_HAP = """
import os
from gpxiltools.cli import main
g = gpxiltools.GpxilFile(P)
c = g[g.contigs[0]]
a = c.length // 2
rc = main(["fasta", P, "--which", "1", "--region", f"{c.name}:{a + 1}-{a + 100_000}", "-o", "/tmp/hap1.fa"])
{"rc": rc, "bytes": os.path.getsize("/tmp/hap1.fa")}
"""
LATE = """
g = gpxiltools.GpxilFile(P)
out = {}
for name in ("chr17", "chrY", g.contigs[-1]):
    if name in g:
        c = g[name]
        out[name] = len(c.sequence(c.length // 2, min(c.length, c.length // 2 + 1000)))
out
"""
VALIDATE = """
import io, contextlib
from gpxiltools.cli import main
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    rc = main(["validate", P, "--no-digests"])
{"rc": rc, "tail": buf.getvalue().strip().splitlines()[-1:]}
"""


def renderer_mb() -> float:
    """Resident memory of Chromium's page processes (the page and its workers live in one), in MB."""
    total = 0
    for proc in psutil.process_iter(["cmdline", "memory_info"]):
        cmd = proc.info["cmdline"] or []
        if any("--type=renderer" in c for c in cmd) and proc.info["memory_info"]:
            total += proc.info["memory_info"].rss
    return total / MB


def server(data: str, port: int) -> subprocess.Popen:
    proc = subprocess.Popen([sys.executable, os.path.join(HERE, "serve.py"), "--data", data, "--port", str(port)], stdout=subprocess.PIPE)
    for _ in range(50):
        try:
            urllib.request.urlopen(f"http://localhost:{port}/index.html", timeout=1)
            return proc
        except Exception:
            time.sleep(0.2)
    raise RuntimeError("server did not start")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--files", default="chr21.gpxil.tif,HG002_chr21.gpxil.tif,chr1.gpxil.tif,GRCh38.gpxil.tif")
    ap.add_argument("--port", type=int, default=8010)
    ap.add_argument("--skip", default="", help="comma-separated: copy,web,validate")
    a = ap.parse_args()
    skip = set(filter(None, a.skip.split(",")))
    files = [f for f in a.files.split(",") if os.path.exists(os.path.join(a.data, f))]
    wheels = sorted(f"wheels/{w}" for w in os.listdir(os.path.join(HERE, "wheels")) if w.endswith(".whl"))
    url = f"http://localhost:{a.port}/index.html"
    results = []
    srv = server(a.data, a.port)
    try:
        with sync_playwright() as p:
            # a persistent profile, as people use: an incognito-style context has a small storage quota
            profile = tempfile.mkdtemp(prefix="gpxil-profile-")
            ctx = p.chromium.launch_persistent_context(profile)

            def fresh():
                page = ctx.new_page()
                page.goto(url)
                page.set_input_files("#files", [os.path.join(a.data, f) for f in files])
                boot = page.evaluate("(w) => call('boot', {wheels: w})", wheels)
                return page, boot

            def run(page, path, code, label, how, name):
                t = time.time()
                try:
                    r = page.evaluate("([p, c]) => call('run', {code: 'P = ' + JSON.stringify(p) + '\\n' + c})", [path, code])
                except Exception as exc:  # noqa: BLE001 - record and go on
                    lines = [ln for ln in str(exc).splitlines() if ln.strip() and not ln.strip().startswith("at ")]
                    r = {"error": lines[-1][:300] if lines else str(exc)[:300]}
                mem = renderer_mb()
                row = {"file": name, "how": how, "step": label, **r, "renderer_mb": round(mem), "wall_s": round(time.time() - t, 2)}
                results.append(row)
                print(f"{name:22s} {how:5s} {label:10s} {r.get('ms', 0) / 1000:7.2f}s  wasm {r.get('wasm', 0) / MB:7.0f} MB"
                      f"  disk {r.get('diskBytes', 0) / MB:7.2f} MB in {r.get('diskReads', 0):5d}  web {r.get('webBytes', 0) / MB:6.2f} MB in {r.get('webRequests', 0):3d}"
                      f"  renderer {mem:6.0f} MB  {r.get('error') or r.get('out')}", flush=True)
                return r

            page, boot = fresh()
            print(f"boot: {boot['ms'] / 1000:.1f}s, wasm {boot['wasm'] / MB:.0f} MB, {boot['versions']}, renderer {renderer_mb():.0f} MB")
            results.append({"step": "boot", **boot})
            page.close()

            for name in files:
                size = os.path.getsize(os.path.join(a.data, name))
                # 1. the file as the browser keeps it, mounted lazily
                page, _ = fresh()
                kept = page.evaluate("(n) => keepAndReturn(chosen(n)).then((b) => b ? b.size : 0)", name)
                print(f"{name:22s} kept in the browser's storage: {kept / MB:.0f} MB")
                page.evaluate("(n) => keepAndReturn(chosen(n)).then((b) => call('mountBlob', {dir: '/mnt/blob', name: n, blob: b || chosen(n)}))", name)
                path = f"/mnt/blob/{name}"
                steps = [("open", OPEN), ("10 kb", REGION_10K), ("1 Mb", REGION_1M), ("late", LATE), ("info", CLI_INFO), ("render", RENDER), ("fasta h1", FASTA_HAP)]
                if "validate" not in skip and size < 200 * MB:
                    steps.append(("validate", VALIDATE))
                for label, code in steps:
                    if label == "fasta h1" and not name.startswith("HG"):
                        continue
                    run(page, path, code, label, "blob", name)
                page.close()
                # 2. the whole file copied into Python's memory, for comparison
                if "copy" not in skip and size < 300 * MB:
                    page, _ = fresh()
                    t = time.time()
                    r = page.evaluate("(n) => call('copyBlob', {dir: '/mem', name: n, blob: chosen(n)})", name)
                    print(f"{name:22s} copy  copied in {time.time() - t:.1f}s: wasm {r['wasm'] / MB:.0f} MB, renderer {renderer_mb():.0f} MB")
                    for label, code in (("open", OPEN), ("10 kb", REGION_10K)):
                        run(page, f"/mem/{name}", code, label, "copy", name)
                    page.close()
                # 3. the file on a web server, read with range requests
                if "web" not in skip:
                    page, _ = fresh()
                    page.evaluate("(n) => call('mountUrl', {dir: '/mnt/web', name: n, url: location.origin + '/data/' + n})", name)
                    for label, code in (("open", OPEN), ("10 kb", REGION_10K), ("1 Mb", REGION_1M), ("late", LATE)):
                        run(page, f"/mnt/web/{name}", code, label, "web", name)
                    page.close()
            ctx.close()
    finally:
        srv.terminate()
    with open(os.path.join(HERE, "results.json"), "w") as fh:
        json.dump(results, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
