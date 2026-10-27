# Tools to work with gpxil: what next

What building GpxilStudio showed is missing, in the order it would help most.
Each is small enough to try here first, in `packages/` or `experiments/`.

## Done for the studio

- **gpxil2pxil** (`packages/gpxil2pxil`): pictures of genes, transcripts and
  proteins, from the reference or either haplotype, with the annotation from
  the Gpxil browser or a GTF. It says how a person's protein differs (P72R).
- **gpxiltools 0.3.0**: files on web servers (range requests, with a HEAD
  request when a browser hides the size), a gpxil file as the reference,
  `vcf --region`, and no imagecodecs needed (Pyodide).
- **pxil**: a `pxil` command, with `-o`.
- **The studio's lazy files**: files of 8 MB or more are read a piece at a
  time, never copied into Python (experiment 01).

## 1. Working copies: `gpxil subset`

Copy some contigs, and some planes, of a file into a new one, without
encoding again: chromosome 17 out of a whole genome, or a file without its
`gq`, `dp` and `ps` planes. Each contig is self-contained (its pages and its
entry in the metadata), so this is a copy of pages, fast and lossless, and
valid under the specification as it is.

A *region* as a working copy (TP53 and 100 kb either side) needs the
specification to say what a contig that is a window of a reference sequence is:
for example `"window": {"of": "chr17", "start": 7561778}` in the contig's
entry, the digests being those of the window's sequence. Worth proposing for
version 0.2 of the specification: it would let the studio keep small files of
the regions people work on, and let the Viewer save what it shows as gpxil.

## 2. `gpxil stats`

What a file holds, counted from its planes, contig by contig and for a region:
bases called and not called, heterozygous and homozygous substitutions,
deletions and insertions, transitions over transversions, the phased fraction,
and the distributions of GQ and DP. Strip by strip, so a genome needs little
memory. JSON for programs, a table for people; the Tools tab and the Viewer (for
the region in view) would show it.

## 3. `gpxil compare`

Two gpxil files on the same reference, compared base by base: where their
genotypes agree, differ, or one has no call, with a picture of the differences
(a plane of codes, drawn like the variants mode). Because the planes are aligned
to the reference, this is an array comparison, without VCF normalisation. It is
what benchmarking a call set against HG002 needs, and a striking picture.

## 4. Genes in the Viewer

The genes of the Gpxil browser's files (or a GTF) drawn over the pixels:
exons, coding sequence, names; pointing at one names it, and clicking opens
gpxil2pxil for it. The Gpxil browser's JavaScript already reads these files. A
BED file as a track too (the callable regions, say).

## 5. Overviews

Whole chromosomes need a smaller picture than one pixel a base: counts of each
colour class in bins (`gpxiltools.render.class_counts` does it), written into
the file as the specification's planned *overview pyramids*, or beside it. The
Viewer could then open any chromosome whole, at once, and zoom in.

## 6. One JavaScript reader

The studio's `gpxil-read.js` and the Gpxil browser's reader do the same thing.
One small library (`gpxil-js`), tested against gpxiltools' files as pxil.py is,
would serve both, and anyone else's web page.

## 7. pxil in the page

The Files tab shows a pxil file as an ordinary picture. Decoding its palette
indices in the page would show it large, with its sequence and metadata beside
it, and a reference protein and a person's side by side, the differences
marked.

## Later

- **Cohorts**: several samples on one reference (the specification's planned
  cohort files), and `gpxil stack` to build one.
- **Variant effects**: every coding change of a gene as HGVS (`p.Pro72Arg`),
  from gpxil2pxil's sequences.
- **Depth from reads**: a `dp` plane from a BAM or CRAM file, for files made
  from VCFs without DP.
