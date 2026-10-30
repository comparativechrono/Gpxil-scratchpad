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
| [03-compare](experiments/03-compare/) | Can genomes be compared base by base, two call sets or a child with both parents, and still see an indel written at another place in a repeat? | Yes: `gpxil compare` in gpxiltools 0.5. All 806 indels moved along repeats are found, and no controls; HG002's v4.2.1 and v5.0q agree on > 99.99 % of variant bases, and HG002 with HG003 and HG004 shows 22 candidate new mutations on chromosome 17 inside all three benchmarks' regions. One artefact of spec 0.1, overlapping unphased calls joined on one haplotype, led to version 0.2, after which every difference agrees with a check from the VCF records. |
| [04-ga-sequences](experiments/04-ga-sequences/) | Albrecht-Buehler's pure GA-sequences, with a person's genome and 76,215 others: where do they come from, and do they really have fewer point mutations? | Many began as the poly(A) tails of Alu elements (37 %, and 71 % of those that start with a long poly(A)): the older the Alu, the more often its tail is a GA-sequence and the less of it is A. New point mutations reach them at about the usual rate (gnomAD's singletons: 1.0 to 1.2 times what their contexts predict), but they have 19-28 times the indels, and their transitions copy a neighbouring unit of their own sequence, the more so the commoner they are (80 % of HG002's transitions; 63 % of gnomAD's singleton transitions to 84 % of its common ones; chance gives 59 %). |
| [05-ga-de-novo](experiments/05-ga-de-novo/) | Do new mutations in pure GA-sequences look like the old ones, and are they "many times fewer"? Found in the long reads of a four-generation family (CEPH 1463), as published short-read de novo sets hardly reach them, and set beside the family's own published calls. | New insertions and deletions hit them 80 to 110 times as often as the genome, a base, in the family's tandem-repeat calls (67 in 8 transmissions: 15 % of all, in 0.14 % of the bases; 17 times in this experiment's stricter calls), the more often the longer the sequence, and add or remove whole units: 28 of 29 of this experiment's are copies of the unit beside them. New substitutions are consistent with the genome's rate (1, where 1.5 are expected) and not 160 times fewer (p = 0.009); that one copies its neighbouring unit. |
| [06-repeat-mutability](experiments/06-repeat-mutability/) | Are pure GA-sequences more mutable than other tandem repeats as long, with units as short, and as pure? Locus by locus, in the 7.3 million loci TRGT genotyped in the family, with its 440 new TR mutations. | No: matched for length, unit and purity they mutate as other repeats do (0.92 times; 1.26 in a Poisson model, p = 0.14), where repeats of A and T mutate twice as often. They are hotspots because they are long, pure repeats, half of the genome's: 3,245 of the 6,546 loci longer than 100 bases with units of up to 6 bases and a purity of 0.85 or more. |
| [07-copied-units](experiments/07-copied-units/) | Do new mutations copy a repeat's own units, as experiment 04's population variants do? The new substitutions of the CEPH families and the Platinum Pedigree in every tandem repeat of the catalogue, with controls, and gnomAD's singletons. | Not as such. At first sight they do (158 of 737 copy, where their contexts predict 96; 2.4 times as often as the same substitution in the same 6-base window where it wouldn't copy), but 82 % of those that copy are alleles gnomAD already has (27 % of the others), and those it lacks copy no more than chance (0.85). gnomAD's singletons copy 1.02 (GA-sequences) to 1.26 (other repeats) times as often, rising to 2 to 2.4 times for common alleles: the copying belongs mostly to alleles seen more than once, and new mutations at large copy little. |
| [08-forensic-markers](experiments/08-forensic-markers/) | Which forensic STR markers are pure GA-sequences, where did they come from, and do they mutate faster? | 7 of 27 lie in GA-sequences (FGA, SE33, D18S51, Penta D, D19S433, D2S1338, D10S1248), and 5 begin at the 3' end of an Alu (D18S51, D19S433, Penta D, D3S1358, Penta E). There is no sign that they mutate faster than the other markers (Liu et al. 2021: Mann-Whitney p = 0.89). |
| [09-alphabets](experiments/09-alphabets/) | Why are the genome's long pure runs purines, when it is 59 % A and T? Pure runs of more than 50 bases in each two-letter alphabet, the Alu and L1 tails they start from, by age, and the spectrum of new substitutions. | Purine runs outnumber A/T runs 2.5 to 1, where composition predicts the reverse. Alu tails start both, but from the youngest Alus to the oldest their pure purine runs rise from 4.8 to 23.8 per 1,000 while their A/T runs fall from 13.4 to 3.5, and in Alu tails only purine runs lengthen. Transitions keep a purine run pure and break the others: 63 % of new substitutions keep an (AAAG)n run pure, 16 % an A/T run. |
| [10-repeat-edges](experiments/10-repeat-edges/) | Do tandem repeats raise substitution in the DNA around them? Every base outside the 7.3 million catalogue loci, by distance from the nearest repeat and its instability, in HG002's two haplotypes and in new mutations. | No broad shadow: from 11 bases away, HG002's variation is what its trinucleotides predict, except as its own genotype at the repeat sets it. But at the base right beside a repeat, substitutions that extend the repeat run at 2.0 to 6.9 times the expected, rising with the repeat's instability, while the others run near the usual rate. The reference makes 1.6 to 1.7 times of that, by deciding where repeats end. New mutations, free of it, extend repeats 3.1 times as often as expected (40 against 13.1), though all of those are short-read calls. |
| [11-along-arrays](experiments/11-along-arrays/) | Where along a repeat array does change happen: do arrays grow or decay from one end? The catalogue's loci at Alu and L1 tails (A-rich, at complete 3' ends), lined up at the element, with the family's and HG002's alleles. | Alu tails mutate like other loci of their length, unit and purity. In the reference they are richest in A next to the Alu at every age, even in the human-specific subfamilies, and lose A all along as they age. Their alleles differ in length next to the Alu 2.7 to 3.5 times as often as at the far end, mostly because their purest stretch lies there: length changes follow the purest stretch wherever it is. Beyond it, Alu tails 4 to 12 % diverged change more at the Alu end (1.75 to 2.54 times, in the family and HG002), an excess that fades with age; L1 tails show none. |
| [12-transcription](experiments/12-transcription/) | Do repeats whose strands differ mutate differently depending on which strand a gene transcribes? Every intronic locus of the catalogue, sense or antisense, against loci between genes, in the family's parents and HG002. | Repeats in genes vary less than repeats between genes (odds 0.84 to 0.97), and least when the polymerase reads their A- or purine-rich strand: sense over antisense, 1.03 for poly(A), in both sets of haplotypes, within Alu tails and outside them; 1.02 to 1.05 for purine runs and 1.08 to 1.11 for A-rich A/T repeats, some intervals including 1; none for G-rich or strand-symmetric repeats. The family's new TR mutations are too few to tell. |
| [13-format](experiments/13-format/) | Is a gpxil file the lossless, compact and quickly read record the software article says, for a whole genome, checked by tools that share none of gpxiltools' code? HG002's v5.0q benchmark encoded on GRCh38, read back and set against bcftools, sized and timed. | Yes. It encodes in 7 minutes into 39.6 MB (the VCF.gz is 42.1 MB) and validates. Its haplotypes match those bcftools consensus builds from the VCF at all 5.35 billion bases compared but 3,720 sites, every one where the VCF contradicts itself. 99.44 % of the VCF's alleles come back from the export as written, the rest the same haplotypes in other records, and the export encodes back to the same haplotypes. A region reads about as fast as with faidx and tabix: a megabase of calls in 1 ms. |

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
