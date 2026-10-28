# Gpxil scratchpad

A place to try out ideas for [gpxil](https://github.com/comparativechrono/Gpxil),
the genome pixel image format, before they become tools of their own.
Anything that works and is new moves to a repository of its own (and to PyPI,
if it is pure Python) and into
[GpxilStudio](https://github.com/comparativechrono/GpxilStudio), the online
workbench for gpxil files.

Nothing here is a release: code is written to answer a question, and each
experiment's README says what the question was, what was found and what was
done with it.

## Experiments

| | Question | Finding |
|---|---|---|
| [01-browser-python](experiments/01-browser-python/) | Can gpxiltools run in the browser on files of hundreds of megabytes? | Yes. Files mounted lazily are read a strip at a time: the 809 MB whole-genome file opens and gives up a megabase of sequence after reading 0.5 MB, in 72 MB of memory. |
| [02-gpxil2pxil](experiments/02-gpxil2pxil/) | Can gene, transcript and protein pictures be made from a gpxil genome, for each of a person's haplotypes? | Yes: [gpxil2pxil](https://github.com/comparativechrono/gpxil2pxil). HG002 turns out homozygous for TP53's Arg72, one pixel different from the reference. |
| [03-compare](experiments/03-compare/) | Can genomes be compared base by base, two call sets or a child with both parents, and still see an indel written at another place in a repeat? | Yes: `gpxil compare` in gpxiltools 0.5. All 806 indels moved along repeats are found, and no controls; HG002's v4.2.1 and v5.0q agree on > 99.99 % of variant bases, and HG002 with HG003 and HG004 shows 22 candidate new mutations on chromosome 17 inside all three benchmarks' regions, and one artefact of the specification. |

## Packages

New tools are written in [`packages/`](packages/), each ready to become a
repository of its own and a PyPI package. When one is ready, it moves to its
repository, with its history:

| Package | What it does | Status |
|---|---|---|
| [gpxil2pxil](https://github.com/comparativechrono/gpxil2pxil) | gene, transcript and protein pxil files from a gpxil genome (reference or haplotypes), with annotation from the Gpxil browser or a GTF | moved to [its own repository](https://github.com/comparativechrono/gpxil2pxil), released on PyPI from 0.1.0; in GpxilStudio |

## Ideas

[`ideas/tools-for-gpxil.md`](ideas/tools-for-gpxil.md): the tools that building
GpxilStudio showed are missing, in the order they would help most. Working
copies (`gpxil subset`), `gpxil stats` and `gpxil compare` are done; next are
genes in the Viewer, overviews, one JavaScript reader, and proteins side by
side, with what the specification's version 0.2 needs.

## Related

- [Gpxil](https://github.com/comparativechrono/Gpxil): the specification
- [Gpxil-tools](https://github.com/comparativechrono/Gpxil-tools): gpxiltools, the Python toolkit
- [pxil](https://github.com/comparativechrono/pxil): single sequences as pictures
- [Gpxil browser](https://comparativechrono.github.io/Gpxil-browser/): the GRCh38 reference as pixels
- [GpxilStudio](https://comparativechrono.github.io/GpxilStudio/): the workbench, with the tools in the browser

## Licence

Copyright © 2026 Chronomic Medicine Group, Department of Genomic Medicine,
University of Cambridge. Released under the [MIT licence](LICENSE).
