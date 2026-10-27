# 01. gpxiltools in the browser's Python, without running out of memory

**Question.** Can gpxiltools run in a web page (Pyodide, Python compiled to
WebAssembly), on gpxil files of hundreds of megabytes, without loading them
into memory? This decides how GpxilStudio handles files.

**Answer: yes.** Mounted the right way, Python reads only the strips of a
file it needs, from the disk or from the web, and a whole-genome file of
850 MB is as quick to work with as an 11 MB one.

## What was tried

`worker.js` runs Pyodide 0.29.5 (Python 3.13) in a web worker and installs
numpy, Pillow and jsonschema from Pyodide, then tifffile 2026.10.15 and
gpxiltools 0.2.0 from their wheels. gpxiltools is installed without its
dependencies, because it asks for imagecodecs, which Pyodide does not have;
tifffile reads and writes deflate, the gpxil default, with Python's own
`zlib`, so nothing is lost but Zstandard.

Each file was given to Python three ways:

| | How | What it costs |
|---|---|---|
| **blob** | The file is kept in the browser's storage (IndexedDB), which returns it as a Blob on disk, and mounted with Emscripten's `WORKERFS`. | Each read is a slice of the Blob, read from disk with `FileReaderSync`. Nothing is copied into memory. |
| **web** | The file stays on a web server and is mounted with `FS.createLazyFile`, which reads it with HTTP range requests, 1 MB at a time. | Only the megabytes touched are downloaded. |
| **copy** | The whole file is copied into Python's in-memory file system, as most in-browser Pythons do. | The file's size in memory, at least. |

`run.py` drives it in headless Chromium (Playwright) with a normal browser
profile, and measures the time, the WebAssembly memory, what was read, and
the resident memory of Chromium's page processes. The files: the GRCh38
reference's chr21 (11 MB) and chr1 (64 MB) from the Gpxil browser, HG002's
chr21 with its genotypes (11 MB, from gpxiltools' worked example), and the
whole GRCh38 reference in one file (809 MB, the Gpxil browser's release).

## Results

Python started in 5 s (packages cached), with 72 MB of WebAssembly memory
and Chromium's page processes at about 350 MB.

Mounted as **blobs**, reading stayed small however big the file:

| File | Open | 10 kb of sequence | 1 Mb | `gpxil info` | `gpxil render` 200 kb | Memory |
|---|---|---|---|---|---|---|
| chr21, 11 MB | 0.01 MB read | 0.08 MB | 0.38 MB | 0.01 MB | 0.09 MB | 72 MB |
| HG002 chr21, 11 MB | 0.01 MB | 0.08 MB | 0.38 MB | 0.01 MB | 0.12 MB | 72 MB |
| chr1, 64 MB | 0.01 MB | 0.03 MB | 0.39 MB | 0.01 MB | 0.05 MB | 72 MB |
| **GRCh38, 809 MB** | **0.05 MB** | **0.06 MB** | **0.42 MB** | **0.05 MB** (194 contigs) | **0.12 MB** | **72 MB** |

Every step took well under half a second. Reaching contigs far into the
whole-genome file (chr17, chrY and the last scaffold) read 0.71 MB in 158
small reads, in 0.06 s: tifffile walks the chain of TIFF directories to find a
page.

Whole-file work reads the whole file, but a strip at a time:
`gpxil validate` read all of chr1 (64 MB) with WebAssembly memory at 124 MB,
and HG002's chr21 with its three planes at 258 MB.

The browser's storage kept the 809 MB file and gave it back as a Blob that
mounted like the others. (In a private window, or Playwright's default
context, the quota is small and the file is not kept: the page must say so
and fall back to the file chosen.)

**Copying** the file in costs its size: chr1 took Chromium's page processes
from 350 to 477 MB before any work was done. That is fine for small files,
and impossible for whole genomes.

Over **the web**, a chromosome file costs 1 to 2 MB to open and 1 to 2 MB a
region, because `createLazyFile` reads whole megabytes. The whole-genome file
is a poor fit for the web: reaching chr17, chrY and the last scaffold took
27 MB, one megabyte for each directory in the chain. The Gpxil browser's
per-chromosome files are the right ones to read over the web.

## What the Studio does with this

- The user's files live in the browser's storage, and Python mounts them as
  blobs: never copied into memory, whatever their size.
- The GRCh38 reference is read over the web, a chromosome file at a time,
  from the Gpxil browser's site (which allows it: CORS and range requests).
  A reader with smaller blocks than `createLazyFile`'s 1 MB, and a cache,
  should cut the downloads further.
- Tools work on regions where they can, and offer working copies (a few
  contigs or regions cut out of a big file) for work on whole planes.
- Python's output goes into its in-memory file system and is moved to the
  browser's storage when a command ends, so outputs should stay modest.

Two things for gpxiltools: it should not require imagecodecs (only Zstandard
needs it), so that `pip install` works in the browser; and `gpxil fasta`
writes to `sys.stdout.buffer`, which a browser console must provide.

## Run it

```sh
./fetch_wheels.sh
pip install playwright psutil && python -m playwright install chromium
python run.py --data DIR       # DIR: gpxil files, e.g. from the Gpxil browser's downloads
```

`serve.py` serves the page with range requests and the headers of
cross-origin isolation. `results.json` holds the numbers of the last run.
