# Stage 1: nested barcode strand-specific RT-qPCR

This is the core CenTraL assay: quantify centromeric transcripts by qPCR,
strand by strand, in a way that is blind to genomic DNA. Stage 2 (Nanopore
sequencing of the same tagged products, see the repo root) is optional and
only needed when you want to place transcripts onto specific HOR arrays.

A full protocol paper is in preparation; this page is the operational
summary.

## Principle in one paragraph

The inner primers are locus- and strand-specific and carry a 5' universal
handle, a short sequence absent from the genome. Reverse transcription is
primed with ONE strand-specific handled primer, so the handle becomes part of
the cDNA and each reaction reads one strand only. A SINGLE PCR cycle with the
complementary handled primer then makes the molecule double stranded with the
handle at both ends. qPCR uses the handle itself as the only primer, and a
molecule amplifies exponentially only with the handle at both ends. Genomic
DNA is never reverse transcribed, so in the single tagging cycle it can gain
the handle at one end at most and stays linear. That is why the tagging step
must remain one cycle: a second cycle would tag genomic DNA at both ends and
break the exclusion.

## Primer design rules

- Each target needs a forward and a reverse locus primer, both carrying the
  same 5' handle. The handle sequence is also the universal qPCR primer.
- The handle must be absent from the reference assembly, in particular from
  alpha-satellite. Confirm in silico before ordering.
- Targets can share a handle (chr2 and chr17 below do); they stay separable
  because each target is reverse transcribed in its own reaction.
- The housekeeping reference (GAPDH by default, or GUSB) sits on a second
  handle so the centromere universal primer can never amplify it.
- The ordered primer is handle + locus sequence, for example chr17 forward is
  CCATGCCGTCGAAACAAGTT-TCGTTCGAAACGGGTATATC.
- Handles carry a slight self-dimer and suppression-PCR risk; heating primers
  at 65 C before use mitigates this.

Handles: centromere targets use CCATGCCGTCGAAACAAGTT; reference genes (GUSB,
ACTB, GAPDH) use ACTTACCGTCGAATTTACCT.

The full primer table (worked chr17/chr2/GUSB set, proposed sets for all
chromosomes, and a Pan-HOR pair) is in [`primers.tsv`](primers.tsv). The
non-worked sets are proposed; validate each against your reference assembly
before use. Do not multiplex chr11 and chr22 in one reaction (shared
sequence). Chromosome assignment ultimately also comes from read mapping in
Stage 2.

## Procedure overview

Work RNase-free through the RNA steps. Filter tips until the end of reverse
transcription. Pre-cool the centrifuge to 4 C.

**A. RNA isolation (about 3 h + precipitation).** TRI Reagent/BCP extraction
from about 2 x 10^5 to 5 x 10^6 cells with a 55 C, 5 min shake step (helps
separate DNA-RNA hybrids), double BCP phase separation, isopropanol +
7.5 M ammonium acetate + GlycoBlue precipitation (-20 C, 2 h to overnight),
two ethanol washes, resuspend in 50 ul nuclease-free water. Do not over-dry
the pellet.

**B. DNase (about 1.5 h).** 3,000 ng RNA in 50 ul; TURBO DNase 1 ul, 37 C
30 min; ADD A SECOND 1 ul and repeat. The alpha-satellite copy number makes a
single DNase round insufficient. Inactivate, spin, keep supernatant,
re-measure.

**C. Strand-specific RT (about 3 h).** Per target set up two reactions, one
with the forward handled primer, one with the reverse, plus a no-RT control
for each. 200 ng RNA + 1 ul of 10 uM primer + 2 ul of 10 mM dNTPs in 10 ul;
65 C 5 min then chill (melts secondary structure of AT-rich centromeric RNA);
add 10 ul MultiScribe RT mix; 37 C 120 min, 85 C 5 min.

**D. Single tagging cycle + primer cleanup (about 1 h).** Add the
COMPLEMENTARY handled primer and PfuUltra II mix to 40 ul. Run ONE cycle
(95 C 1 min; 95 C 20 s, 60 C 20 s, 72 C 20 s; 72 C 3 min). Keep it to one
cycle. Then digest leftover single-stranded primers with Exonuclease I
(10 U, 37 C, 3 min; inactivate 80 C, 1 min). Do not skip the ExoI step;
leftover handled primers prime in qPCR and inflate background.

**E. qPCR (about 1.5 h).** Heat universal primers at 65 C before use.
20 ul reactions: 10 ul 2X SensiFAST SYBR No-ROX, 300 nM universal primer,
about 10 ng template. Cycle 95 C 2 min; 40 x (95 C 5 s, 62 C 10 s, 72 C
20 s); melt 50 to 72 C. Include +RT and -RT for every target plus the GUSB
reference. The +RT to -RT Ct gap should be at least 5 cycles; treat take-off
below 0.3 as no amplification.

## Quantification

Average technical replicates. dCt = Ct(target) - Ct(reference). Pick a calibrator
condition; ddCt = dCt(sample) - mean dCt(calibrator). Relative level =
2^-ddCt. Test significance on dCt values (t test or one-way ANOVA).

GAPDH is the default reference. In ASO experiments use GUSB instead: ASO
treatment had an off-target effect on the GAPDH transcript in our hands, so
GAPDH is unreliable there.

[`analyze_ddct.py`](analyze_ddct.py) does this from a Ct table and checks
the +RT/-RT gap:

```bash
python3 qpcr/analyze_ddct.py qpcr/example_ct.csv \
  --calibrator untreated --out results_qpcr
# GAPDH is the default reference; add --reference GUSB for ASO experiments
```

## Troubleshooting (qPCR stage)

| Problem | Likely cause | Fix |
|---|---|---|
| Signal in -RT | Residual gDNA or handled-primer carryover | Repeat/extend DNase; confirm ExoI ran; heat primers |
| No signal in +RT | Degraded/dilute RNA, RT failure, secondary structure | Check integrity and yield; confirm 65 C denature+chill; standardise input |
| Primer dimers | Universal primer self-priming | Heat primers at 65 C; lower concentration |
| Fwd = Rev signal | Real bidirectional transcription or barcode cross-talk | Check primer assignment; run single-primer controls |

## Continuing to Stage 2

The same single-cycle tagged product feeds sequencing: amplify 31 cycles
with the universal primers only (genomic DNA still cannot gain both
handles), SPRI clean, end-prep, native barcode per sample (each strand and
the reference control gets its own barcode), pool, adaptor ligation with Short
Fragment Buffer washes, and sequence on a MinION/PromethION with live
basecalling and live alignment in MinKNOW. Then run the pipeline in the repo
root: `./run_dcs_workflow.sh /path/to/bam_pass`.
