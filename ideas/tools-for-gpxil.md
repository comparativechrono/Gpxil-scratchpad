# Tools to work with gpxil: what next

What building GpxilStudio showed is missing, in the order it would help most.
Each is small enough to try here first, in `packages/` or `experiments/`.

## Done for the studio

- **gpxil2pxil** ([its repository](https://github.com/comparativechrono/gpxil2pxil)): pictures of genes, transcripts and
  proteins, from the reference or either haplotype, with the annotation from
  the Gpxil browser or a GTF. It says how a person's protein differs (P72R).
- **gpxiltools 0.3.0**: files on web servers (range requests, with a HEAD
  request when a browser hides the size), a gpxil file as the reference,
  `vcf --region`, and no imagecodecs needed (Pyodide).
- **pxil**: a `pxil` command, with `-o`.
- **The studio's lazy files**: files of 8 MB or more are read a piece at a
  time, never copied into Python (experiment 01).
- **pxil in the page**: the studio's Files tab decodes a pxil file's palette
  indices itself and shows it large, with its sequence and metadata, naming
  the residue under the pointer; pxil's FASTA names a person's sequences
  apart (`ENSP00000269305_HG002-hap1`).
- **gpxiltools 0.4.0**: `gpxil subset` (working copies of some contigs and
  planes, the strips copied as stored) and `gpxil stats` (what a file holds,
  counted; in the studio's Tools tab and under the Viewer's key).
- **gpxiltools 0.5.0**: `gpxil compare`, two genomes base by base, or a child
  with both parents (experiment 03).

## 1-3. Done: `gpxil subset`, `gpxil stats`, `gpxil compare`

In gpxiltools 0.4 and 0.5, and the studio. What they showed the specification
needs, for version 0.2:

- **A region as a working copy** (TP53 and 100 kb either side): a contig that
  is a window of a reference sequence, for example
  `"window": {"of": "chr17", "start": 7561778}` in the contig's entry, the
  digests being those of the window's sequence. It would let the studio keep
  small files of the regions people work on, and let the Viewer save what it
  shows as gpxil.
- **Overlapping unphased calls** (experiment 03): two heterozygous calls that
  overlap can't be on one haplotype, but spec 0.1 applies both as written, to
  haplotype 2, and joins them. Like the rule that moves an unphased `*`
  allele to the haplotype with the deletion, the second could go to
  haplotype 1 where that is free.

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
One reader, shared by the two (not published on npm), tested against
gpxiltools' files as pxil.py is, to come with the work on version 0.2 of the
specification.

## 7. Proteins side by side

The studio shows one pxil picture at a time. A reference protein and a
person's, side by side and aligned, with the differences marked, would make
gpxil2pxil's notes (P72R) visible at a glance.

## Later

- **Cohorts**: several samples on one reference (the specification's planned
  cohort files), and `gpxil stack` to build one.
- **Variant effects**: every coding change of a gene as HGVS (`p.Pro72Arg`),
  from gpxil2pxil's sequences.
- **Depth from reads**: a `dp` plane from a BAM or CRAM file, for files made
  from VCFs without DP.
