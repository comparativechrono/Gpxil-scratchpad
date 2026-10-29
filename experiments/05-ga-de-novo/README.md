# 05. New mutations in pure GA-sequences, in a four-generation family

**Question.** [Experiment 04](../04-ga-sequences/) found that pure
GA-sequences carry 19 to 28 times the genome's insertions and deletions, and
that their substitutions are mostly transitions, which copy a neighbouring
unit of the sequence the more often the commoner the variant: 63 % of
gnomAD's singleton transitions, 84 % of its common ones, where chance gives
59 %. But a polymorphism is a mutation that survived. Do new mutations
already look like that? And what do they say of Albrecht-Buehler's "many
times fewer point mutations" (160 times, before his correction)?

Published de novo mutations hardly reach GA-sequences. Of the 28,725
autosomal de novo mutations of the CEPH three-generation families (Sasani et
al. 2019, from short reads), 5 fall inside, where 40 would at the genome's
density, while their 400-base flanks get 297, where 285 would: the reads, or
the call set's filters, leave GA-sequences out. So here they are found anew,
in the long reads of a four-generation family, CEPH 1463 (the
[Platinum Pedigree](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets),
[Porubsky et al. 2025](https://doi.org/10.1038/s41586-025-08922-2)), whose
reads span GA-sequences whole.

**Answers.**

1. **New insertions and deletions hit GA-sequences about 17 times as often as
   the genome, a base.** 24 de novo indels in 14 transmissions, 1.7 a
   transmission, where the genome's rate gives 1.4 in all: 16.8 times (95 %
   interval 10.8 to 25). The calls are strict, so this is probably an
   underestimate. It is of the order of experiment 04's polymorphisms:
   HG002's indels, 19 to 24 times the genome's, and gnomAD's indel alleles,
   28 times the background's (22 times for its singletons).
2. **They add or remove whole units of the repeat.** Of the 29 simple
   insertions and deletions (de novo, or mosaic in a parent), 28 are one or
   more copies of the unit beside them, mostly a single 4-base unit (GAAA,
   TTTC, GGAA ...): 16 lose one, 8 gain one. Losses outnumber gains, 19 to
   10; 7 more changes are complex. Of the 31 whose parent can be told, 24
   arose on the father's allele (77 %; the family's paper finds 75 to 81 %
   genome-wide).
3. **New substitutions are consistent with the genome's rate, and not 160
   times fewer.** One in 14 transmissions, where the genome's rate gives
   1.46: 0.7 times (95 % interval 0.017 to 3.8). At 160 times fewer, one
   would be seen with a probability of 0.009; ten times fewer can't be ruled
   out (0.14). That one is a transition (G>A on the purine strand) that
   makes its 6 bases a 6-mer found elsewhere in the allele: it copies its
   neighbouring unit, GGAA. Two more, less sure (a parent's mosaic, and one
   whose other allele isn't seen), are transitions too, and one of them
   copies.
4. **Longer GA-sequences mutate more**, from 0.016 mutations per 1,000
   GA-sequences a transmission at 51 to 75 bases to 0.17 at 101 to 150 (de
   novo and mosaic): ten times; longer ones are too few to say more. Alu
   tails carry 18 of the 38 events (47 %, where 37 % of GA-sequences are
   Alu tails; p = 0.13).

**What it means.** New mutations in GA-sequences are mostly a unit gained
or lost, as in microsatellites, where slippage copies or drops a unit; new
substitutions are rare there, as anywhere, and not the 160 times rarer of
the first count. Whether new substitutions copy their neighbours as often as
the ones that spread do (experiment 04) one family can't say: at the
genome's rate a transmission gives 0.1, so some hundred would take about a
thousand long-read trios.

![New mutations in one family's GA-sequences](figures/rates.png)

## How it was done

The family's public data (AWS Open Data, `s3://platinum-pedigree-data`,
GRCh38): the HiFi joint calls of its 22 open members (DeepVariant and
GLnexus), and each member's HiFi reads. 16 children have both parents among
them: 2 in G2, whose parents' DNA is from cell lines only (G1 is no longer
living), and 5 in G3 and 9 in G4, whose parents' DNA is from blood
([UW Medicine](https://newsroom.uw.edu/news-releases/4-generations-help-science-explore-genome-mutation-rate)).

| Script | What it does | Output |
|---|---|---|
| [`ped.py`](ped.py) | the family, its trios, and where its files are | |
| [`short_reads.py`](short_reads.py) | how many of the CEPH three-generation families' short-read de novo mutations fall in GA-sequences (lifted to GRCh38) | `results/short_reads.json` |
| [`call_dnms.py`](call_dnms.py) | candidates genome-wide: an allele in a child, in neither parent (no read of it, GQ ≥ 20, depth 10 to 2.5 times the median), in no one else but the child's siblings and descendants | `results/candidates.tsv.gz`, `results/candidates.json` |
| [`validate.py`](validate.py) | every candidate in a GA-sequence, and 150 elsewhere, read whole in the family's HiFi reads: each read that spans the repeat, between two anchors of 15 bases found once nearby, for the child, parents, siblings, children and grandparents | `results/read_alleles.json.gz` |
| [`alleles.py`](alleles.py) | each person's alleles from those reads (by length, a base either way joined, then by sequence), and the verdict: is the child's new allele in neither parent? | `results/validated.tsv` |
| [`analyse.py`](analyse.py) | classes, what each new allele is against the parental allele closest to it, rates, mutability by length | `results/dnms_ga.tsv`, `results/analysis.json` |
| [`figures.py`](figures.py) | the pictures | `figures/` |

A candidate is **de novo** when the child's new allele (5 reads or more) is
in neither parent (at most a read of it, of 8 or more), the child's other
allele is in a parent (or the new one was passed on), no sibling has it, and,
when the child has three children or more, one of them has it. It is
**mosaic in a parent** when a parent has a few reads of it (fewer than a
fifth of theirs). A one-base change in a run of one base is taken for a read
error, so such new alleles are missed. Fetched twice, the reads gave the same
alleles and verdicts.

## Results

### Candidates

3,185 candidate alleles in the 16 trios: 81 substitutions a trio (the
paper: 74.5 a transmission), and 118 insertions and deletions, most of them
not new. Of 150 candidate sites checked the same way outside GA-sequences
(162 carriers), 43 of the 54 substitution verdicts with enough reads are
confirmed (24 of 26 in G4), but only 14 of 101 insertion and deletion
verdicts: the rest are inherited alleles a parent was called without (73),
or have a few reads in a parent (14); 25 of them are one-base changes in
runs of one base, which the reads can't settle. 117 candidate alleles (107
sites) fall in GA-sequences.

### Verdicts in GA-sequences

| | G3 and G4 (14 trios) | G2 (2 trios) |
|---|---|---|
| de novo | 25 | 8 |
| de novo, other allele not seen | 4 | |
| mosaic in a parent | 13 | 12 |
| shared by siblings | 4 | |
| too few reads of it | 2 | |
| not passed on (to any of 5 children) | | 1 |
| not new, or unclear | 30 | 18 |

The rates use G3 and G4 only: G1's DNA is from cell lines, where repeats can
drift in culture, and 12 of G2's 20 new alleles have a few reads in a parent
(G1). In G3 and G4, no new allele that could be tested for transmission
failed: NA12879's three were each passed on to two or three of her five
children.

| Per base, against the genome's rate (the paper's, per transmission) | GA-sequences | expected at the genome's rate | ratio (95 %) |
|---|---|---|---|
| insertions and deletions (7.4 + 65.3) | 24 | 1.43 | 16.8 (10.8 to 25.0) |
| substitutions (74.5) | 1 | 1.46 | 0.7 (0.017 to 3.8) |

### What the new alleles are

Each new allele against the parental allele closest to it, read as one
change (38 events of G3 and G4, de novo and mosaic):

| | |
|---|---|
| deletions | 19 (one unit of 4 bases: 16) |
| insertions | 10 (one unit of 4 bases: 8) |
| complex (both, or noisy reads) | 7 |
| substitutions | 2 |
| insertions and deletions that are copies of the unit beside them | 28 of 29 (de novo 18 of 19, mosaic 10 of 10) |
| arisen in the father / the mother / can't tell | 24 / 7 / 7 |
| in an Alu tail / simple repeat / low complexity / other | 18 / 15 / 3 / 2 |

The units: GAAA 4, TTTC 3, TCTT 3, TTCT 3, AAGA, CTTT, GGAA, AAAG and CCTT 2
each, TCCT and AGGA 1 each (read on the + strand); two 3-base units, two
copies of an 8-base unit lost at once, and one 24-base deletion that copies
nothing.

[figures/family_unit_deletion.png](figures/family_unit_deletion.png): a
unit lost. The father's 181-base allele went to three of the five children
with reads; in one of them it is one GAAA shorter.

![A unit lost](figures/family_unit_deletion.png)

[figures/family_substitution.png](figures/family_substitution.png): the de
novo substitution. The father's 158-base allele went to all five children
with reads; in one of them, one G is an A (on the purine strand), and the
run of GGAA units next to it is one unit longer.

![A new base](figures/family_substitution.png)

### By length

| GA-sequence length | GA-sequences | mutations (G3, G4; de novo and mosaic) | per 1,000 a transmission (95 %) |
|---|---|---|---|
| 51-75 | 22,245 | 5 | 0.016 (0.005 to 0.037) |
| 76-100 | 8,210 | 8 | 0.070 (0.030 to 0.137) |
| 101-150 | 7,079 | 17 | 0.172 (0.100 to 0.275) |
| 151-250 | 3,904 | 6 | 0.110 (0.040 to 0.239) |
| 251 or more | 724 | 2 | 0.197 (0.024 to 0.713) |

![Longer GA-sequences mutate more](figures/by_length.png)

## Caveats

- One family: 14 transmissions, and 1 new substitution. The numbers are
  small, and the longest GA-sequences have fewest reads that span them.
- The calls are strict: candidates from one caller, then confirmed in the
  reads. One-base changes in runs of one base are missed, and so are alleles
  the caller didn't propose, so the counts, and the ratio to the genome, are
  probably underestimates. But a call rests on one technology's reads, and
  allows a read of the new allele in a parent; and the genome's rates are
  the family's paper's, from a more thorough pipeline.
- Mosaic events (a few reads in a parent) are counted apart from de novo
  ones, and shared ones (in two siblings or more, and in neither parent)
  left out: a parent's germline mosaic, or a parental allele missed.
- G2's new alleles may be G1's cell lines drifting; they are set apart.
- The parent a change arose in is the one whose allele is closest to the new
  one; when both parents' are as close, it is left open.

## Reproducing

```sh
# in a folder of its own: the joint calls (1.6 GB)
curl -O https://platinum-pedigree-data.s3.amazonaws.com/variants/small_variants/GRCh38/CEPH1463.GRCh38.deepvariant.glnexus.oa.vcf.gz{,.tbi}
python3 -I call_dnms.py CEPH1463.GRCh38.deepvariant.glnexus.oa.vcf.gz 2        # 4 min
CURL_CA_BUNDLE=... python3 -I validate.py REF_DIR 4 150    # 22 min, from the family's BAM files; run where their
                                                           # indexes (430 MB) can be saved
python3 -I alleles.py && python3 -I analyse.py && python3 figures.py
python3 -I short_reads.py ceph-dnm-manuscript/data GRCh37_to_GRCh38.chain.gz   # Sasani et al.'s files; Ensembl's chain
```

REF_DIR is the Gpxil browser's GRCh38 gpxil files, as in experiment 04,
whose `results/` this one reads. Python 3.12 or later, with gpxiltools 0.7,
NumPy, SciPy, Matplotlib, Pillow, pysam, cyvcf2 and pyliftover.
