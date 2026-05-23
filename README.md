# CenTraL

_**Cen**tromeric **Tra**nscript capturing and **L**ong read sequencing — a one-command analysis pipeline for per-barcode Nanopore cDNA amplicon data from centromere regions._

End-to-end post-sequencing pipeline for **per-barcode Nanopore cDNA amplicon data**
targeting centromeric transcripts: merge chunked BAMs → DCS spike-in normalization
→ primary-only filtering → normalized bigwigs → 1 bp read-start BEDs + bedgraphs
→ QC plot. A companion R script (`scripts/karyoplot_bedgraph.R`) makes
per-chromosome karyoplots of the bedgraphs with centromere / HOR highlighting.

Tailored for **centromere / repeat work** — α-satellite arrays, HOR regions, and
the pericentromeric segmental duplications around them. No MAPQ filtering by
default, so reads with ambiguous repeat-array placement are kept (they're real
signal, not noise, when the whole point of your experiment is to read those
regions).

---

## What this pipeline starts from

CenTraL is a **post-sequencing** pipeline. It assumes you've already:

1. **Prepared the library with the ONT native barcoding kit** (e.g. SQK-NBD114).
   The DCS lambda spike-in that comes standard with that kit is what makes the
   per-sample DCS normalization in step 2 possible.
2. **Sequenced on a Nanopore device with basecalling + alignment running live
   in MinKNOW** (target reference loaded at run start), so each barcode's reads
   are already basecalled and aligned by the time the run ends.
3. **Ended the run with the standard `bam_pass/barcodeNN/*.bam` layout** —
   chunked, already-aligned, per-barcode BAMs.

CenTraL picks up from there: it merges the chunked BAMs into one BAM per
barcode and runs every downstream step. It does **not** do basecalling,
demultiplexing, or alignment — those are MinKNOW's job during sequencing.

---

## Quickstart

After `git clone`-ing this repo onto your cluster:

```bash
# 1. Make sure the required tools are on $PATH — see the Dependencies section.
#    Use whatever you already use: conda env, system modules, container,
#    system install. The pipeline doesn't care. The orchestrator does a
#    preflight check at startup and aborts cleanly if something's missing.

# 2. Run the whole post-sequencing pipeline in ONE command
cd /scratch/myproject     # wherever you want results to land
sbatch /path/to/CenTraL/run_dcs_workflow.sh /path/to/your/bam_pass
# or interactive:
/path/to/CenTraL/run_dcs_workflow.sh /path/to/your/bam_pass
```

Where `/path/to/your/bam_pass` is the standard MinKNOW/Dorado output dir with
`barcodeNN/` subfolders of chunked `*.bam` files. The script does everything
from there.

---

## Repo layout

```
CenTraL/
├── README.md                          you are here
├── LICENSE                            MIT
├── .gitignore                         keeps run outputs out of git
├── environment.yml                    conda env spec for reproducibility
├── run_dcs_workflow.sh                ← ENTRY POINT — one command, does steps 1-6
├── scripts/                           worker scripts (don't run these yourself unless you want a single step)
│   ├── barcode_merge_nanopore.sh        step 1 — merge per-barcode chunks
│   ├── count_dcs_spikein.sh             step 2 — DCS spike-in count + scale factor
│   ├── filter_primary_bams.sh           step 3 — primary-only (-F 2308) BAMs
│   ├── make_normalized_bigwigs.sh       step 4 — normalized bigwigs
│   ├── make_readstart_bed.sh            step 5 — 1 bp BED + bedgraph per barcode
│   ├── plot_normalization_qc.py         step 6 — DCS QC bar plot
│   ├── karyoplot_bedgraph.R             step 7 (manual) — chromosome karyoplot of bedgraphs
│   └── DCS_Lambda_3.6kb.fa              ONT DCS reference (auto-resolved by step 2)
└── examples/
    ├── demo_dcs_counts.tsv             example step 2 output
    └── dcs_normalization_qc_DEMO.png   example step 6 output
```

---

## What you get out

`run_dcs_workflow.sh /path/to/bam_pass` produces:

```
/path/to/bam_pass/merged_bam/                 (created in place, kept for re-use)
   barcode01.bam + .bai
   barcode02.bam + .bai
   ...

./dcs_analysis_<timestamp>/                   (created in your CWD)
├── dcs_counts.tsv                            per-barcode counts + scale factors
├── dcs_logs/                                 per-barcode aligner logs (dorado/minimap2)
├── primary_bams/                             -F 2308 filtered + indexed BAMs
│   └── *.bam + *.bam.bai                      (clean input for IGV / downstream tools)
├── bw/                                       DCS-normalized bigwigs
│   └── *.bw + *.bamCoverage.log
├── readstart_beds/                           1 bp BED + collapsed bedgraph per barcode
│   ├── *.readstart.bed.gz                    one row per read (BED6: chr/start/end/name/mapq/strand)
│   └── *.startcount.bedgraph                 chr/start/end/count — drop into karyoploteR
├── dcs_normalization_qc.png + .pdf           before/after normalization QC
└── run.log                                   full transcript of this run
```

---

## Flags for `run_dcs_workflow.sh`

| Flag | Meaning | Default |
|------|---------|---------|
| (positional) | `bam_pass/` directory — required | — |
| `-i PATH` | Same as positional — `bam_pass/` directory | top-of-script value |
| `-o PATH` | Output dir for analysis results | `./dcs_analysis_<timestamp>/` |
| `-t N` | Threads | `$SLURM_CPUS_PER_TASK` if set, else `4` |
| `-s N` | bigwig bin size (bp) | `50` |
| `-u` | Include `unclassified` barcode in DCS/bw/BED | off (it's not a real sample) |
| `-f` | Force overwrite of any existing outputs | off (idempotent by default) |
| `-h` | Show help | |

**Examples:**

```bash
# Default — let the script create a timestamped folder in CWD
./run_dcs_workflow.sh /path/to/bam_pass

# Specify output dir
./run_dcs_workflow.sh -i /path/to/bam_pass -o /scratch/myrun

# Bigger bigwig bins for whole-chromosome views, more threads
./run_dcs_workflow.sh -i /path/to/bam_pass -s 100 -t 16

# Re-run and overwrite everything (e.g. after fixing something)
./run_dcs_workflow.sh -f /path/to/bam_pass

# Include the unclassified pool (default skips it since it's not a real sample)
./run_dcs_workflow.sh -u /path/to/bam_pass
```

---

## Step 7 — per-chromosome karyoplot (separate, project-specific)

The orchestrator stops at step 6 because step 7 needs **your project's**
chrom-sizes and centromere/HOR BEDs, which vary per genome build. Run it
manually for any subset of barcodes you want to compare:

```bash
Rscript /path/to/CenTraL/scripts/karyoplot_bedgraph.R \
    /path/to/chrom.sizes \
    /path/to/centromere_horAll.bed     \   # sharp HOR highlight (or NA)
    /path/to/centromere_recent.bed     \   # broader centromere backdrop (or NA)
    plots/chr1_compare \                   # output prefix (PDF + PNG)
    chr1 \                                 # target_chr — or "all" for genome-wide
    auto \                                 # zoom — "auto" uses backdrop > regions > data extent
    dcs_analysis_*/dcs_counts.tsv \        # DCS TSV for normalization, or NA for raw counts
    dcs_analysis_*/readstart_beds/barcode01.startcount.bedgraph \
    dcs_analysis_*/readstart_beds/barcode05.startcount.bedgraph \
    dcs_analysis_*/readstart_beds/barcode09.startcount.bedgraph
```

**Karyoplot args (7 fixed, then 1+ bedgraphs):**

| # | Arg | Notes |
|---|-----|-------|
| 1 | `chrom.sizes` | UCSC `chr<TAB>length` |
| 2 | `regions`     | sharp HOR / α-satellite arrays BED — drawn dark rose on ideogram. `NA` to skip. |
| 3 | `backdrop`    | broader centromere context BED — drawn faint pink behind the HOR. `NA` to skip. |
| 4 | `out_prefix`  | writes `<prefix>.pdf` and `<prefix>.png`. Parent dir auto-created. |
| 5 | `target_chr`  | `all` for genome-wide, or e.g. `chr5`. |
| 6 | `zoom`        | `auto` (prefers backdrop > regions > data extent) / `full` / `chr:start-end`. Ignored for `all`. |
| 7 | `dcs_tsv`     | `dcs_counts.tsv` for DCS-normalized tracks (heights comparable across samples). `NA` for raw counts. |
| 8+ | bedgraphs    | one or more `*.startcount.bedgraph` from step 5 — order determines top-to-bottom track order. |

Use `NA` for any optional slot you don't have.

---

## Before you start — dependencies

Everything below must be on `$PATH` **before** you run the pipeline. The
scripts don't activate an env for you (you do that once in your login shell).

| Tool | Min version | Used by | Why |
|------|-------------|---------|-----|
| **bash** | 4+ | all | basic shell |
| **samtools** | 1.10+ | all bash scripts | merge / sort / index / view / fastq |
| **dorado** *or* **minimap2** | dorado 0.5+ / minimap2 2.20+ | step 2 | align unmapped reads to DCS. Auto-detected — dorado preferred. |
| **deepTools** | 3.5+ | step 4 (`bamCoverage`) | bigwig generation |
| **awk** | any | several | TSV / SAM column math |
| **gzip** | any | step 5 | compress BED output |
| **Python** | 3.9+ | step 6 | QC plot |
| **matplotlib** | 3.7+ | step 6 | plotting |
| **numpy** | 1.24+ | step 6 | numeric ops |
| **R** | 4.0+ | step 7 (optional) | karyoplot — installs `regioneR` + `karyoploteR` on first run |

### Get the dependencies on PATH

Use whatever you already use for managing scientific software — the pipeline
doesn't care HOW you do it, only that the tools above are on `$PATH` when you
invoke `run_dcs_workflow.sh`. A few common options:

- **Conda env** — there's a bundled `environment.yml` you can build:
  ```bash
  conda env create -f environment.yml
  conda activate central
  ```
- **System modules** (Lmod / Environment Modules) — load the modules your
  cluster provides, e.g.:
  ```bash
  module load samtools dorado deeptools python R
  ```
- **Mix-and-match** — e.g. a conda env for samtools/deepTools/Python +
  `module load dorado` for the aligner.
- **System install** — if your sysadmin already installed everything globally,
  you don't need to do anything.

Either way, the orchestrator runs a preflight check at startup and **aborts
cleanly before SLURM allocates anything** if any tool is missing — so a wrong
env costs nothing.

### Sanity check

Paste this once to confirm everything is in place:

```bash
samtools --version | head -1
command -v dorado >/dev/null && dorado --version || minimap2 --version
bamCoverage --version
python3 -c "import matplotlib, numpy; print('OK')"
Rscript -e 'suppressMessages(library(karyoploteR)); cat("OK\n")'
```

If all five lines print without errors, you're ready to go.

The orchestrator runs the same dependency checks at startup and **aborts
cleanly before SLURM allocates anything** if any tool is missing.

---

## Reference data

- **Your target genome** — already used by MinKNOW during basecalling /
  alignment, so the merged BAMs from step 1 already carry the alignment. You
  don't need a separate FASTA for steps 1-6.
- **DCS reference** — shipped at `scripts/DCS_Lambda_3.6kb.fa` (3,587 bp of
  ONT's DNA Control Strand, lambda phage). Step 2 finds it automatically. If
  the file is missing for any reason, step 2 will recreate it from a sequence
  embedded in `scripts/count_dcs_spikein.sh` — repo is self-contained.

For **karyoplots (step 7)** you additionally need a chrom-sizes file and
ideally centromere / HOR BEDs for your build. UCSC `fetchChromSizes` works
for chrom-sizes; HOR / centromere coordinates come from your assembly
annotation (T2T, your custom build, etc.).

---

## Pipeline design notes (in case you want to know why)

- **PCR amplicon data** — every read is a PCR product; PCR duplicates are
  the readout of amplification yield, not an artifact. The pipeline never
  deduplicates.
- **DCS spike-in normalization** — `scale = min(dcs_across_samples) /
  dcs_this_sample`. Smallest spike-in sample gets `scale = 1.0`; others get
  scaled down so per-amplicon depths are comparable across samples.
- **`-F 2308` filter** at step 3 keeps only primary mapped reads (drops
  unmapped + secondary + supplementary). One alignment per read, no
  double-counting. Read-start BEDs and bigwigs are derived from this.
- **No MAPQ filtering**. For centromere / repeat regions, reads with
  ambiguous repeat-copy placement (MAPQ 0-1) are real signal. Filtering them
  out for "cleanliness" throws away the data you care about.
- **Idempotent** — re-running any step skips work that's already done unless
  you pass `-f`. So you can incrementally re-run only what changed.
- **No time/RAM caps** in the SBATCH headers (`--time=0 --mem=0`) — the job
  takes whatever the partition allows. If your cluster doesn't honour
  `--time=0`, set it to your partition's `MaxTime` explicitly.

---

## SLURM defaults

`#SBATCH` directives baked in at the top of every `.sh` script:

| Directive         | Value | Why                                       |
|-------------------|-------|-------------------------------------------|
| `--cpus-per-task` | `8`   | Threads passed to all tools               |
| `--mem`           | `0`   | All RAM on the allocated node — no cap    |
| `--time`          | `0`   | No wall-time limit — run to completion    |

About `--time=0`: most clusters interpret it as "no limit". On a few it
falls back to the partition's `DefaultTime`. If a job gets killed early,
run `scontrol show partition <name>` and replace `--time=0` with the
`MaxTime` value (e.g. `--time=14-00:00:00`).

---

## Examples directory

`examples/` has:

- `demo_dcs_counts.tsv` — what a typical step 2 output looks like
- `dcs_normalization_qc_DEMO.png` — what a typical step 6 plot looks like

---

## License

MIT — see `LICENSE`.
