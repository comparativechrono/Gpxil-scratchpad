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

## Related

- [Gpxil](https://github.com/comparativechrono/Gpxil): the specification
- [Gpxil-tools](https://github.com/comparativechrono/Gpxil-tools): gpxiltools, the Python toolkit
- [pxil](https://github.com/comparativechrono/pxil): single sequences as pictures
- [Gpxil browser](https://comparativechrono.github.io/Gpxil-browser/): the GRCh38 reference as pixels

## Licence

Copyright © 2026 Chronomic Medicine Group, Department of Genomic Medicine,
University of Cambridge. Released under the [MIT licence](LICENSE).
