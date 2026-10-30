# 11. Where along a repeat array does change happen?

**Question.** A tandem repeat is an array of units, and it changes by
gaining and losing them, and by substitutions that break them. Does that
happen anywhere along the array, or at one end: do arrays grow, or decay,
from one end? An array has no ends to tell apart unless something marks
one. The tail of an Alu or an L1 element does: it begins at the element's
3' end. Lined up there, as Albrecht-Buehler stacked GA-sequences, every
tail array has a near end (next to the element) and a far end. gpxil's
planes give each array's bases on the element's strand, and a person's two
haplotypes base for base beside them.

**Answers.**

1. **TE tails hold almost half of the family's new TR mutations.** Of the
   catalogue's 7,270,921 autosomal loci, 760,533 begin at an Alu's 3' end
   and 252,379 at an L1's, by the rule experiments 04 and 09 use (within 30
   bases, on the element's strand). The Alu tails are 10.5 % of loci and
   hold 154 of the 440 new TR mutations, 35 %: 1.37 times what loci of their
   length, motif size and purity have elsewhere (1.16 to 1.61). The L1
   tails hold 41, 1.14 times (0.82 to 1.54). The strata leave out
   composition, which experiment 06 found matters (AT-rich repeats mutate
   about twice as much); an Alu tail is mostly A.
2. **In the reference, Alu tails are purest next to the Alu, at every age.**
   In the youngest Alu tails (0 to 4 % divergence), 97.5 % of the bases in
   the tenth next to the Alu are A, against 81 to 84 % in the last two
   tenths. In the oldest (20 % or more), 73.9 % against 60 to 66 %. With
   age every tenth loses A alike: the near tenth falls by 24 points, the
   middle by 23, the far one by 18. The tails decay along their length; the
   gradient is there in the youngest Alus already, so it was set when they
   arrived, not by decay from one end. In experiment 09's pure purine runs
   at Alu tails, the share of G, the letter a transition brings, rises
   from 4 % next to the Alu to 27 % at the far end in the youngest Alus,
   and from 33 % to 40 % (44 % mid-run) in the oldest. A/T runs, with T,
   behave the same way. L1 tails show no gradient.
3. **Length changes happen where the array is purest, and an Alu tail's
   purest stretch lies next to the Alu.** Where an allele differs in length
   from the reference, the comparison places the difference in a stretch of
   the array; in 102,466 of the 264,812 differences of the family's parents
   (TRGT) at Alu tails, that stretch is at most half the array. Of those,
   58,265 lie in the half next to the Alu and 25,549 in the far half, 2.3
   to 1, on either strand alike (so not an artefact of how variants are
   placed). HG002's haplotypes, called independently, give 80,898 against
   27,210, 3.0 to 1. But the difference follows the array's longest stretch
   in frame (its purest; for a poly(A), its longest run of one base). Where
   that stretch lies in the near half, the family's differences fall there
   4.4 to 1; where it lies in the far half, 0.64 to 1, that is in the far
   half; and the differences outside it split almost evenly (0.86 to 1).
   They touch the purest stretch 1.48 times as often as stretches as wide
   placed at random would (1.67 in HG002; 1.44 over all 355,505 localised
   differences at any catalogue locus). At L1 tails the two halves are
   even (1.05 and 1.08 to 1).
4. **As the near end decays, the polarity fades.** The family's near-to-far
   ratio at Alu tails is 2.1 to 2.8 in Alus up to 12 % divergent, 2.09 at
   12 to 16 %, 1.75 at 16 to 20 % and 1.31 beyond (HG002: 2.9 to 3.7, then
   2.69, 2.11 and 1.48).

**What it means.** Arrays change where they are purest: interruptions are
known to slow slippage (Ananda et al. 2014), and here the changes in the
alleles of two families of haplotypes follow the purest stretch wherever it
lies. Alu tails look as if they change at one end only because their purest
stretch lies next to the Alu: their poly(A) begins there, and their far end
is less pure from the start, perhaps because it is where reverse
transcription began on the target's sequence. Neither growth nor decay goes
from one end: substitutions wear the tails down all along, and the
polarity of length change fades as the near end loses its purity. The
landmark is what makes this visible. Without the Alu to line them up on, an
array's two ends can't be told apart, and polarity can't be asked about.

![Alu tails stacked at the Alu, by age; the share of A in each tenth of a tail array, by age; and the share of length differences that lie in the half of the array next to the element](figures/along_arrays.png)

## How it was done

| Script | What it does | Output |
|---|---|---|
| `tails.py` | Finds the catalogue's loci that begin at an Alu's or an L1's 3' end (RepeatMasker; within 30 bases, on the element's strand), reads each on the element's strand from GRCh38's gpxil files (A in each tenth; the longest stretch in frame), and the second letter of experiment 09's pure tail runs by tenth; new TR mutations at TE tails against experiment 06's strata | `results/tail_loci.tsv.gz`, `results/tails.json` |
| `lengths.py` | The parents' alleles (TRGT's open-access calls for the Platinum Pedigree) and HG002's haplotypes (GIAB v5.0q, from its gpxil files: the gt plane, the insertion flags and the sidecar's inserted sequences) that differ in length from the reference at each tail array; where along the array each difference can lie, read from the element | `results/windows.tsv.gz`, `results/lengths.json` |
| `stacks.py` | gpxiltools' stack: the Alu tails lined up at their near end, youngest Alu first; each column's bases by age | `figures/stack_alu_tails.png`, `results/stack_rows.npy`, `results/stack_profile.json` |
| `figures.py` | The figure | `figures/along_arrays.png` |

The per-array tables, `results/tail_loci.tsv.gz` (31 MB) and
`results/windows.tsv.gz` (7 MB), are written by the scripts and not kept in
the repository.

**Placing a length difference.** An allele that differs from the reference
in length differs somewhere in a stretch of the reference that the two
sequences can't narrow further. For one insertion or deletion, it is the
span of the places the event could be put; otherwise, the bases between
the longest start and the longest end the two share. In a pure array the
stretch is the whole array. A difference is localised when its stretch is
at most half the array, near when the stretch lies in the half next to the
element, far when in the other half. A stretch placed at random would be
near as often as far, and would touch the longest stretch in frame with the
chance `lengths.py` computes for each one.

## Results

New TR mutations of the Platinum Pedigree (eight G3 children) at TE tails,
against loci of the same length, motif size and purity elsewhere:

| Loci | Loci (share) | Mutations (share) | Expected | Ratio (95 % CI) |
|---|---|---|---|---|
| Alu tails | 760,533 (10.5 %) | 154 (35.0 %) | 112.1 | 1.37 (1.16 to 1.61) |
| L1 tails | 252,379 (3.5 %) | 41 (9.3 %) | 36.0 | 1.14 (0.82 to 1.54) |
| either | 1,012,912 (13.9 %) | 195 (44.3 %) | 145.3 | 1.34 (1.16 to 1.54) |

Alu tail arrays in the reference, by the Alu's divergence (share of A by
tenth, from the tenth next to the Alu to the far one):

| Divergence | Arrays | Mean length | A, near tenth | A, middle (6th) | A, far tenth | Purest stretch: share of array |
|---|---|---|---|---|---|---|
| 0-4 % | 9,064 | 38.8 | 0.975 | 0.849 | 0.844 | 0.774 |
| 4-8 % | 90,042 | 37.8 | 0.943 | 0.790 | 0.802 | 0.694 |
| 8-12 % | 304,714 | 39.3 | 0.909 | 0.751 | 0.773 | 0.650 |
| 12-16 % | 219,979 | 43.7 | 0.853 | 0.703 | 0.733 | 0.615 |
| 16-20 % | 113,513 | 45.8 | 0.797 | 0.662 | 0.697 | 0.600 |
| 20 %+ | 23,221 | 52.6 | 0.739 | 0.622 | 0.660 | 0.579 |

The far tenth rises a little above the eighth and ninth (0.81 against 0.83
in the youngest): the catalogue's loci tend to end on a unit.

Length differences at tail arrays (localised: the stretch they could lie in
is at most half the array):

| | Differences | Localised | Near half | Far half | Near : far |
|---|---|---|---|---|---|
| Alu, family | 264,812 | 102,466 | 58,265 | 25,549 | 2.28 |
| Alu, HG002 | 322,083 | 127,252 | 80,898 | 27,210 | 2.97 |
| L1, family | 63,341 | 29,243 | 11,635 | 11,091 | 1.05 |
| L1, HG002 | 58,407 | 27,699 | 11,471 | 10,656 | 1.08 |

Alu tails, by where the array's purest stretch lies (family; HG002 in
brackets):

| | Localised | Near half | Far half |
|---|---|---|---|
| purest stretch in the near half | 72,670 (97,691) | 49,134 (72,864) | 11,290 (11,403) |
| purest stretch in the far half | 29,796 (29,561) | 9,131 (8,034) | 14,259 (15,807) |
| differences not touching it | 17,179 (15,702) | 6,209 (5,688) | 7,190 (6,701) |

The near-to-far ratio at Alu tails by the Alu's divergence, family (HG002):
0-4 % 2.11 (2.89); 4-8 % 2.83 (3.65); 8-12 % 2.64 (3.45); 12-16 % 2.09
(2.69); 16-20 % 1.75 (2.11); 20 %+ 1.31 (1.48).

## Caveats

- The arrays are the catalogue's loci, whose ends TRsolve chose; an Alu
  tail's far end may take in the target's sequence. The rule that makes an
  array a tail allows its start to lie up to 30 bases from the end that
  RepeatMasker gives the element, whose consensus takes in part of the
  poly(A).
- Most differences can't be placed: half of them could lie anywhere in at
  least 69 % of the array. The placed ones come from interrupted arrays,
  where slippage is slower, and where an allele carries more than one event
  the stretch is the core between them.
- The family's open-access TRGT calls cover the 537,486 loci where one of
  seven people (the G2 couple and five G3 children) has a non-reference
  allele; they hold only one of the 461 loci where the family's paper
  called a new TR mutation that passed its filters, so the new mutations
  themselves can't be placed. The parents' and HG002's alleles
  differ from the reference by history, not one generation's mutation.
- HG002's alleles come from GIAB's small-variant benchmark, whose calls in
  long homopolymers are least sure; only arrays it calls on both
  haplotypes, with the base either side, are counted (968,809 of
  1,012,912).
- The purest stretch is the longest stretch in which each base equals the
  base one shortest motif on; for a poly(A), the longest run of one base.

## Sources

- Albrecht-Buehler G. Properties and distribution of pure GA-sequences of
  mammalian genomes. PLoS ONE. 2008;3(11):e3818.
  doi:[10.1371/journal.pone.0003818](https://doi.org/10.1371/journal.pone.0003818).
- Ananda G, Hile SE, Breski A, Wang Y, Kelkar Y, Makova KD, Eckert KA.
  Microsatellite interruptions stabilize primate genomes and exist as
  population-specific single nucleotide polymorphisms within individual
  human genomes. PLoS Genet. 2014;10(7):e1004498.
  doi:[10.1371/journal.pgen.1004498](https://doi.org/10.1371/journal.pgen.1004498).
- Roy-Engel AM, et al. Active Alu element "A-tails": size does matter.
  Genome Res. 2002;12(9):1333-1344.
  doi:[10.1101/gr.384802](https://doi.org/10.1101/gr.384802).
- Porubsky D, et al. Human de novo mutation rates from a four-generation
  pedigree reference. Nature. 2025;643:427-436.
  doi:[10.1038/s41586-025-08922-2](https://doi.org/10.1038/s41586-025-08922-2)
  (TRGT calls: `ceph_1463_tandem_repeats.oa.vcf.gz`, the Platinum
  Pedigree's open-access data).
- Platinum TRs v1.0 (Dashnow, [Zenodo 13178746](https://doi.org/10.5281/zenodo.13178746), CC BY 4.0).
- GIAB HG002 v5.0q small-variant benchmark (GRCh38).
- UCSC Genome Browser, RepeatMasker annotation of hg38 (`rmsk.txt.gz`).

## Reproducing

```sh
python3 -I tails.py REF_DIR RMSK FEATURES                   # 2 min
python3 -I lengths.py REF_DIR HG002_DIR TRGT_VCF FEATURES   # 4 min
python3 -I stacks.py GENOME                                 # 3 min
python3 figures.py
```

REF_DIR is the Gpxil browser's GRCh38 gpxil files and GENOME the same as
one file. RMSK is UCSC's `rmsk.txt.gz` for hg38. FEATURES is experiment 06's
locus measures (its `--features` file). HG002_DIR holds HG002_chrN.gpxil.tif
from GIAB's v5.0q small-variant VCF, made as experiment 04's
`encode_hg002.sh` makes them. TRGT_VCF is the Platinum Pedigree's
`ceph_1463_tandem_repeats.oa.vcf.gz`. `tails.py` reads experiment 09's
`results/runs.tsv.gz`. Python 3.12 or later, with gpxiltools 0.8 or later,
NumPy, pandas and Matplotlib.
