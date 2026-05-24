# CenTraL

_**Cen**tromeric **Tra**nscript capturing and **L**ong read sequencing._

**Status:** v1.0 alpha — in active use, API may change before the paper is published.

A small post-sequencing pipeline for **per-barcode Nanopore cDNA amplicon data** from centromeric regions. Written for our centromere RNA work in RPE1 cells — native-barcoded libraries on a MinION/PromethION with live basecalling and alignment in MinKNOW.

After the run, you have a `bam_pass/` folder full of chunked per-barcode BAMs. CenTraL turns that into:

- one merged BAM per barcode,
- DCS spike-in counts + scale factors for cross-sample normalization,
- DCS-normalized bigwigs for IGV,
- 1-bp-per-read BEDs + bedgraphs (for karyoploteR / pyGenomeTracks),
- a QC plot of before/after normalization,
- and (optionally) per-chromosome karyoplots with HOR / centromere highlighting.

**The whole thing is one command:**

```bash
./run_dcs_workflow.sh /path/to/your/bam_pass
```

---

## Before you run

CenTraL doesn't basecall, demultiplex, or align — MinKNOW does that live during sequencing. CenTraL assumes:

1. Library was prepared with the **ONT native barcoding kit** (e.g. SQK-NBD114). DCS lambda spike-in comes standard with that kit — that's what step 2 uses for normalization.
2. The run was set up with **live basecalling + live alignment** to your reference in MinKNOW.
3. The run ended with the standard `bam_pass/barcodeNN/*.bam` chunked layout. Reads are already basecalled and aligned.

If you have something different (e.g. only fast5/pod5, or unaligned BAMs), you'll need to basecall and align with dorado yourself before pointing CenTraL at the result.

---

## Install

CenTraL needs these tools on `$PATH`. Nothing else — the pipeline doesn't care HOW you install them.

| Tool | What for |
|---|---|
| `samtools` (≥1.10) | merge / sort / index / view / fastq |
| `dorado` *or* `minimap2` (≥2.20) | re-align unmapped reads to DCS |
| `deepTools` ≥3.5 (`bamCoverage`) | normalized bigwigs |
| `python3` ≥3.9, `matplotlib`, `numpy` | QC plot |
| R ≥4.0 + `karyoploteR` + `regioneR` | karyoplot (step 7, optional — auto-installs on first run) |
| `bash`, `awk`, `gzip` | standard Unix (always pre-installed) |

### Get the tools — pick one

**Option A — one-shot conda install**

```bash
conda create -y -n central \
    -c bioconda -c conda-forge \
    samtools minimap2 deeptools \
    python=3.10 matplotlib numpy \
    r-base bioconductor-karyoploter bioconductor-regioner

conda activate central
```

(`central` is just the env name I picked — pass `-n whatever_you_like` to use a different one. The pipeline doesn't read the env name; it only looks for the tools on `$PATH`.)

**Option B — from the bundled `environment.yml`**

```bash
git clone https://github.com/santoshbiowarrior333/CenTraL.git
cd CenTraL
conda env create -f environment.yml      # creates env "central"
conda activate central
```

**Option C — system modules (HPC clusters)**

```bash
module load samtools dorado deeptools python R
pip install --user matplotlib numpy      # only if your Python doesn't have them
```

**Option D — mix and match**

E.g. conda env for everything Python/R + `module load dorado` for the aligner. Works fine.

### Sanity check

Paste this once to confirm all five are reachable:

```bash
samtools --version | head -1
command -v dorado >/dev/null && dorado --version || minimap2 --version
bamCoverage --version
python3 -c "import matplotlib, numpy; print('python OK')"
Rscript -e 'suppressMessages(library(karyoploteR)); cat("R OK\n")'
```

Five lines, five OKs → you're ready. If something's missing, the pipeline will tell you exactly which tool at startup before doing any work.

---

## Quickstart

```bash
# 1. clone the repo (one-time)
git clone https://github.com/santoshbiowarrior333/CenTraL.git

# 2. activate whatever env you set up above
conda activate central          # or module load …, etc.

# 3. run the whole pipeline in one command
cd /scratch/myrun               # wherever you want results to land
./CenTraL/run_dcs_workflow.sh /path/to/your/bam_pass

# or as a SLURM job (no time/RAM caps set by default):
sbatch ./CenTraL/run_dcs_workflow.sh /path/to/your/bam_pass
```

Results land in `./dcs_analysis_<timestamp>/` in your current directory. The merged BAMs from step 1 land alongside your `bam_pass/` (so you can re-use them without re-running the merge).

---

## What it does, step by step

The orchestrator (`run_dcs_workflow.sh`) calls these in order. You can also run any single step from `scripts/` if you want to redo just one thing.

### Step 1 — merge chunked per-barcode BAMs

Your raw `bam_pass/` looks like this after a MinKNOW run:

```
bam_pass/
├── barcode01/   *.bam, *.bam.bai   (often hundreds of chunks)
├── barcode02/   *.bam, *.bam.bai
├── ...
└── unclassified/   *.bam, *.bam.bai
```

Step 1 pipes the chunks of each folder through `samtools merge → samtools sort → samtools index` (no big intermediate file on disk) and produces:

```
bam_pass/merged_bam/
├── barcode01.bam + .bam.bai
├── barcode02.bam + .bam.bai
└── ...
```

Run just this step:

```bash
./scripts/barcode_merge_nanopore.sh -i /path/to/bam_pass -t 8
```

### Step 2 — count DCS spike-in reads

The DCS (DNA Control Strand — a 3.6 kb lambda phage fragment) added during native-barcoding library prep doesn't map to your target reference, so DCS reads sit in the unmapped fraction of each merged BAM.

For each barcode, step 2 pulls out the unmapped reads and re-aligns them against `scripts/DCS_Lambda_3.6kb.fa`:

```bash
samtools view -b -f 4 barcode01.bam \
  | dorado aligner DCS_Lambda_3.6kb.fa -                 \
  | samtools view -c -F 2308 -        # primary mapped only
```

Output: `dcs_counts.tsv` with columns

```
barcode | total_reads | target_mapped | unmapped | dcs_mapped | scale_factor
```

Scale factor is `min(dcs > 0) / dcs_this_sample`, so the smallest-spike-in sample gets `1.0` and every other sample is scaled down — depths become comparable across samples.

Run just this step:

```bash
./scripts/count_dcs_spikein.sh -b /path/to/bam_pass/merged_bam -t 8
```

(`dorado` is picked automatically when it's on PATH; `minimap2 -ax map-ont` is the fallback.)

### Step 3 — filter to primary mapped reads

`samtools view -F 2308` drops unmapped + secondary + supplementary alignments. One alignment per read, no double-counting from split reads. Writes `primary_bams/`:

```bash
./scripts/filter_primary_bams.sh -b /path/to/bam_pass/merged_bam -t 8
```

### Step 4 — DCS-normalized bigwigs

`bamCoverage --scaleFactor <per-barcode factor>` per barcode. 50 bp bins by default; bump down to 10 bp if your amplicons are short, or up to 100+ for whole-genome views.

```bash
./scripts/make_normalized_bigwigs.sh \
    -b dcs_analysis_*/primary_bams \
    -c dcs_analysis_*/dcs_counts.tsv \
    -s 50 -t 8
```

Drop the `.bw` files straight into IGV — heights are comparable across samples.

### Step 5 — 1-bp-per-read BEDs + bedgraphs

For each primary-mapped read, one tiny entry at its leftmost coordinate:

```
barcode01.readstart.bed.gz       # one row per read (BED6 with strand)
barcode01.startcount.bedgraph    # chr/start/end/count — collapsed for karyoploteR
```

Useful when you want a clean "where does each amplicon fire?" view without overlapping intervals cluttering the plot.

```bash
./scripts/make_readstart_bed.sh -b dcs_analysis_*/primary_bams -t 8
```

### Step 6 — QC plot

Two-panel matplotlib bar chart — DCS counts up top, raw vs normalized target-mapped reads below. Quick visual confirmation that the normalization actually equalized your samples.

```bash
./scripts/plot_normalization_qc.py -c dcs_analysis_*/dcs_counts.tsv
```

### Step 7 — chromosome karyoplots (manual)

Not in the orchestrator because it needs **your** chrom.sizes + centromere/HOR BEDs (varies per genome build). Run it for whatever subset of barcodes and whatever chromosome you want:

```bash
Rscript scripts/karyoplot_bedgraph.R \
    /path/to/hg38.chrom.sizes \
    /path/to/centromere_horAll.bed \            # sharp HOR (or NA)
    /path/to/centromere_broad.bed \             # faint backdrop (or NA)
    plots/chr1_six_samples \                    # output prefix → PDF + PNG
    chr1 \                                       # or "all" for genome-wide
    auto \                                       # zoom: auto / full / chr:start-end
    dcs_analysis_*/dcs_counts.tsv \             # DCS TSV (or NA = raw counts)
    dcs_analysis_*/readstart_beds/barcode01.startcount.bedgraph \
    dcs_analysis_*/readstart_beds/barcode02.startcount.bedgraph \
    dcs_analysis_*/readstart_beds/barcode05.startcount.bedgraph \
    dcs_analysis_*/readstart_beds/barcode09.startcount.bedgraph \
    dcs_analysis_*/readstart_beds/barcode11.startcount.bedgraph \
    dcs_analysis_*/readstart_beds/barcode15.startcount.bedgraph
```

The script installs `karyoploteR` + `regioneR` on first run if they're not already there. Output is a stacked multi-track plot — one track per barcode, color-coded, sharing a y-axis (so heights are directly comparable), with the centromere/HOR highlighted on the chromosome ideogram.

---

## All the flags for `run_dcs_workflow.sh`

```bash
./run_dcs_workflow.sh -h
```

Most common patterns:

```bash
# default
./run_dcs_workflow.sh /path/to/bam_pass

# more threads, smaller bigwig bins (for short amplicons)
./run_dcs_workflow.sh /path/to/bam_pass -t 16 -s 10

# specify where the analysis folder lands
./run_dcs_workflow.sh /path/to/bam_pass -o /scratch/myrun

# include the "unclassified" pool (default skips it — it's not a real sample)
./run_dcs_workflow.sh /path/to/bam_pass -u

# force re-run, overwriting any existing outputs
./run_dcs_workflow.sh /path/to/bam_pass -f
```

---

## Output layout

```
bam_pass/merged_bam/                      ← created in place by step 1; kept for re-use
    barcode01.bam + .bam.bai
    ...

./dcs_analysis_<timestamp>/                ← created in your CWD by the orchestrator
├── dcs_counts.tsv                         per-barcode counts + scale factors
├── dcs_logs/                              per-barcode aligner logs (dorado/minimap2)
├── primary_bams/                          -F 2308 filtered + indexed BAMs
│   ├── barcode01.bam + .bam.bai
│   └── ...
├── bw/                                    DCS-normalized bigwigs (for IGV)
│   ├── barcode01.bw + .bamCoverage.log
│   └── ...
├── readstart_beds/                        per-read BED + collapsed bedgraph
│   ├── barcode01.readstart.bed.gz
│   └── barcode01.startcount.bedgraph
├── dcs_normalization_qc.png + .pdf        QC plot (before/after normalization)
└── run.log                                full transcript of this run
```

---

## Repo layout

```
CenTraL/
├── README.md                              you are here
├── LICENSE                                MIT
├── .gitignore                             keeps run outputs out of git
├── environment.yml                        conda env spec for reproducibility
├── run_dcs_workflow.sh                    ← ENTRY POINT — one command, does steps 1–6
├── scripts/
│   ├── barcode_merge_nanopore.sh          step 1 — merge per-barcode chunks
│   ├── count_dcs_spikein.sh               step 2 — DCS spike-in count + scale factor
│   ├── filter_primary_bams.sh             step 3 — primary-only (-F 2308) BAMs
│   ├── make_normalized_bigwigs.sh         step 4 — normalized bigwigs
│   ├── make_readstart_bed.sh              step 5 — 1 bp BED + bedgraph per barcode
│   ├── plot_normalization_qc.py           step 6 — DCS QC bar plot
│   ├── karyoplot_bedgraph.R               step 7 — chromosome karyoplot (manual)
│   └── DCS_Lambda_3.6kb.fa                ONT DCS reference (auto-resolved by step 2)
└── examples/
    ├── demo_dcs_counts.tsv                example step 2 output
    └── dcs_normalization_qc_DEMO.png      example step 6 output
```

---

## Design choices worth knowing

**No MAPQ filtering.** Reads in pericentromeric regions multi-map across HOR copies and segmental duplications — that's real signal for centromere work, not noise. Filtering on MAPQ would throw away the data you actually care about.

**No deduplication.** Every read in PCR amplicon data is a PCR product by definition. The duplicates ARE the amplification signal. Removing them would destroy the readout.

**DCS-based normalization corrects for what it can.** DCS proxies for library prep + sequencing variation (everything that happens after DCS is added). It does **not** correct for PCR amplification efficiency differences (which happen before).

**Idempotent.** Re-running any step skips work that's already done. Pass `-f` to force.

**SLURM with no resource caps.** `#SBATCH --mem=0 --time=0` in every script header — the job uses the whole node and doesn't get killed early.

---

## Citation — required

If you use CenTraL in any published work, you **must** cite it.

The companion manuscript is **in preparation**. The full citation (authors,
journal, DOI) will be added to this section as soon as the paper is out —
please re-check this README before submitting your manuscript.

For now, link to the repository in your Methods:

> https://github.com/santoshbiowarrior333/CenTraL

shashisantosh2007@gmail.com
path1327@ox.ac.uk


## License

MIT — see `LICENSE`. You're free to use, modify, and redistribute the code;
the only firm condition is the citation requirement above.

---

In active use for centromere RNA work. Issues and PRs welcome.
