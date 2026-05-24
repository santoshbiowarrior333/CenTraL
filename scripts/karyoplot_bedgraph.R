#!/usr/bin/env Rscript
# Karyoplot of one OR many bedgraph tracks. Drop in one bedgraph to get a
# single-track plot; drop in several to get them stacked above each chromosome
# bar (each labelled, sharing the same y-axis for easy comparison).
#
# Same input pattern as the old fwd/rev primer plot, just with bedgraphs:
#   - chrom sizes file (UCSC chr<TAB>length)
#   - regions BED for centromere / HOR highlight (or NA)
#   - one or more bedgraphs (chr/start/end/count)
#
# Usage:
#   Rscript karyoplot_bedgraph.R <chrom.sizes> <regions> <backdrop> <out_prefix> \
#                                <target_chr> <zoom> <dcs_tsv> \
#                                <bedgraph1> [bedgraph2 ...]
#
# Args (seven fixed, then 1+ bedgraphs at the end):
#   chrom.sizes : chr<TAB>length file (UCSC fetchChromSizes hg38 > hg38.chrom.sizes)
#   regions     : BED of SHARP regions to highlight (HOR / specific centromere arrays).
#                 Drawn on top of the ideogram. Pass "NA" to skip.
#   backdrop    : BED of BROADER regions to draw faintly behind `regions` on the
#                 ideogram (e.g. whole centromere area). Same layered look as your
#                 old fwd/rev script. Pass "NA" to skip.
#   out_prefix  : output filename prefix (writes <prefix>.pdf and <prefix>.png).
#                 Parent dir is auto-created.
#   target_chr  : "all" for genome-wide (every chr stacked), or a chromosome
#                 name (e.g. "chr17") for a single-chr plot.
#   zoom        : "auto" / "full" / "chr:start-end". Ignored when target_chr=all.
#                 Pass "NA" if not relevant.
#   dcs_tsv     : path to dcs_counts.tsv to multiply each track by its
#                 scale_factor — heights become directly comparable across
#                 samples. Pass "NA" for raw counts (no normalization).
#   bedgraph*   : one or more bedgraph paths. Track name comes from the
#                 filename ("barcode05.startcount.bedgraph" -> "barcode05").
#
# Examples:
#   # Same shape as your old fwd/rev plot — one barcode, chr2, with HOR
#   # (sharp) layered over a broader centromere backdrop (faint)
#   Rscript karyoplot_bedgraph.R \
#       ../chrom.sizes \
#       ../../../../centromere_horAll.bed \
#       ../../../../../recent_rpe1/centromere_recent_rpe1.bed \
#       chr2_plot_primers chr2 auto NA \
#       chr2_starts.startcount.bedgraph
#
#   # Pick 3 barcodes, genome-wide, DCS-normalized
#   Rscript karyoplot_bedgraph.R \
#       hg38.chrom.sizes centromere_horAll.bed centromere_broad.bed \
#       plots/three_genome all NA merged_bam/dcs_counts.tsv \
#       readstart_beds/barcode01.startcount.bedgraph \
#       readstart_beds/barcode05.startcount.bedgraph \
#       readstart_beds/barcode09.startcount.bedgraph
#
#   # All 18 barcodes, zoomed to chr5p, DCS-normalized
#   Rscript karyoplot_bedgraph.R \
#       hg38.chrom.sizes chr5_HOR.bed NA \
#       plots/all18_chr5p chr5 chr5:1-50000000 merged_bam/dcs_counts.tsv \
#       readstart_beds/barcode*.startcount.bedgraph

# ---- user library + install missing packages ----
user_lib <- Sys.getenv("R_LIBS_USER", unset = file.path(Sys.getenv("HOME"), "R", "library"))
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
cran_repo <- "https://cloud.r-project.org"
install_if_missing <- function(pkg, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible())
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", lib = user_lib, repos = cran_repo)
    BiocManager::install(pkg, lib = user_lib, ask = FALSE, update = FALSE)
  } else {
    install.packages(pkg, lib = user_lib, repos = cran_repo)
  }
}
suppressPackageStartupMessages({
  install_if_missing("BiocManager")
  install_if_missing("regioneR",    bioc = TRUE)
  install_if_missing("karyoploteR", bioc = TRUE)
  library(karyoploteR); library(regioneR); library(GenomicRanges)
})

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8)
  stop("Usage: Rscript karyoplot_bedgraph.R <chrom.sizes> <regions> <backdrop> <out_prefix> <target_chr> <zoom> <dcs_tsv> <bedgraph1> [bedgraph2 ...]")

chrom_sizes_file <- args[1]
regions_file     <- args[2]
backdrop_file    <- args[3]
out_prefix       <- args[4]
target_chr       <- args[5]
zoom_arg         <- args[6]
dcs_tsv_file     <- args[7]
bedgraph_files   <- args[8:length(args)]

all_chr       <- target_chr %in% c("all", "ALL", "*", "genome")
skip_regions  <- regions_file  %in% c("NA", "-", "none", "NONE")
skip_backdrop <- backdrop_file %in% c("NA", "-", "none", "NONE")
skip_dcs      <- dcs_tsv_file  %in% c("NA", "-", "none", "NONE")

cat("Tracks to plot: ", length(bedgraph_files), "\n", sep = "")

# ---- chromosome sizes ----
sizes <- read.table(chrom_sizes_file, header = FALSE,
                    col.names = c("chr", "length"), stringsAsFactors = FALSE)
chr_order <- function(x) {
  num <- suppressWarnings(as.numeric(sub("^chr", "", x)))
  ifelse(is.na(num),
         99 + match(toupper(sub("^chr", "", x)), c("X","Y","M","MT")),
         num)
}
sizes <- sizes[order(chr_order(sizes$chr)), ]
if (!all_chr) {
  sizes <- sizes[sizes$chr == target_chr, ]
  if (nrow(sizes) == 0)
    stop("target_chr '", target_chr, "' not found in ", chrom_sizes_file)
}
custom_genome <- toGRanges(data.frame(chr = sizes$chr, start = 1, end = sizes$length))

# ---- optional DCS scale factor lookup ----
scale_lookup <- list()
if (!skip_dcs) {
  dcs <- read.table(dcs_tsv_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  if (!all(c("barcode", "scale_factor") %in% colnames(dcs)))
    stop("DCS TSV must have columns 'barcode' and 'scale_factor': ", dcs_tsv_file)
  for (i in seq_len(nrow(dcs))) {
    sf <- suppressWarnings(as.numeric(dcs$scale_factor[i]))
    if (!is.na(sf)) scale_lookup[[ dcs$barcode[i] ]] <- sf
  }
  cat("DCS scale factors loaded for ", length(scale_lookup), " barcode(s)\n", sep = "")
}

# ---- read each bedgraph; apply scale factor if available ----
track_names <- character(length(bedgraph_files))
track_data  <- vector("list", length(bedgraph_files))

for (i in seq_along(bedgraph_files)) {
  f <- bedgraph_files[i]
  if (!file.exists(f)) stop("Bedgraph not found: ", f)
  bc <- sub("\\.bedgraph$", "",
       sub("\\.startcount\\.bedgraph$", "", basename(f)))
  track_names[i] <- bc
  bg <- read.table(f, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("chr", "start", "end", "count"))
  if (!all_chr) bg <- bg[bg$chr == target_chr, ]
  bg <- bg[bg$chr %in% sizes$chr, ]
  if (nrow(bg) == 0) {
    cat("  ", bc, ": NO entries on target — track will be blank\n", sep = "")
    track_data[[i]] <- list(gr = GRanges(), count = numeric(0), name = bc, sf = 1)
    next
  }
  sf <- if (bc %in% names(scale_lookup)) scale_lookup[[bc]] else 1
  scaled <- bg$count * sf
  cat("  ", bc, ": ", nrow(bg), " positions, count ",
      round(min(scaled), 2), "-", round(max(scaled), 2),
      "  (scale=", round(sf, 4), ")\n", sep = "")
  track_data[[i]] <- list(
    gr    = toGRanges(bg[, 1:3]),
    count = scaled,
    name  = bc,
    sf    = sf
  )
}

nonempty <- Filter(function(t) length(t$count) > 0, track_data)
if (length(nonempty) == 0) stop("All tracks empty — nothing to plot.")
# ymax is computed AFTER we know the zoom — see below.

# ---- optional region highlight ----
read_regions <- function(path) {
  if (!file.exists(path)) return(GRanges())
  first <- trimws(readLines(path, n = 1, warn = FALSE))
  if (grepl("^[^:]+:[0-9]+-[0-9]+$", first)) {
    regs <- read.table(path, header = FALSE, stringsAsFactors = FALSE)$V1
    m    <- regmatches(regs, regexec("^([^:]+):(\\d+)-(\\d+)$", regs))
    df   <- do.call(rbind, lapply(m, function(x) if (length(x) == 4)
              data.frame(chr = x[2], start = as.numeric(x[3]),
                         end = as.numeric(x[4]), stringsAsFactors = FALSE)))
  } else {
    df <- read.table(path, header = FALSE, stringsAsFactors = FALSE)[, 1:3]
    colnames(df) <- c("chr", "start", "end")
  }
  if (!all_chr) df <- df[df$chr == target_chr, ]
  df <- df[df$chr %in% sizes$chr, ]
  if (nrow(df)) toGRanges(df) else GRanges()
}
regions  <- if (!skip_regions)  read_regions(regions_file)  else GRanges()
backdrop <- if (!skip_backdrop) read_regions(backdrop_file) else GRanges()
cat("region intervals: ", length(regions),
    "   backdrop intervals: ", length(backdrop), "\n", sep = "")

# ---- zoom (single-chr only) ----
zoom_gr <- NULL
if (!all_chr) {
  pad <- 2e5
  if (zoom_arg %in% c("full", "NA", "-", "none", "NONE")) {
    zoom_gr <- custom_genome
  } else if (zoom_arg == "auto") {
    # Prefer the user-supplied region files (they signal where they care):
    #   1) backdrop BED (broader centromere) - if present, zoom to its extent on target_chr
    #   2) regions BED (sharp HOR) - else zoom to its extent
    #   3) fall back to the extent of the data itself
    if (length(backdrop) > 0) {
      ref_gr <- backdrop
      src    <- "backdrop"
    } else if (length(regions) > 0) {
      ref_gr <- regions
      src    <- "regions"
    } else {
      ref_gr <- unlist(GRangesList(lapply(nonempty, function(t) t$gr)))
      src    <- "data"
    }
    if (length(ref_gr) == 0) {
      zoom_gr <- custom_genome
    } else {
      z_start <- max(1,            min(start(ref_gr)) - pad)
      z_end   <- min(sizes$length, max(end(ref_gr))   + pad)
      zoom_gr <- toGRanges(data.frame(chr = target_chr,
                                      start = z_start, end = z_end))
      cat("auto-zoom from ", src, ": ", target_chr, ":",
          z_start, "-", z_end, "\n", sep = "")
    }
  } else {
    mm <- regmatches(zoom_arg, regexec("^([^:]+):(\\d+)-(\\d+)$", zoom_arg))[[1]]
    if (length(mm) != 4) stop("Bad zoom string: ", zoom_arg, " (use chr:start-end)")
    zoom_gr <- toGRanges(data.frame(chr = mm[2], start = as.numeric(mm[3]),
                                    end = as.numeric(mm[4])))
  }
}

# widen narrow bars for visibility at this zoom
zoomspan <- if (!is.null(zoom_gr)) {
              end(zoom_gr) - start(zoom_gr) + 1
            } else {
              sum(as.numeric(width(custom_genome)))
            }
draw_w   <- max(round(zoomspan / 1500), 50)

# --- Clip each track to the visible window so the y-axis reflects ONLY
# what is actually shown in the plot (otherwise a spike outside the zoom
# would dominate the scale and flatten the in-view signal).
clip_gr <- if (!is.null(zoom_gr)) zoom_gr else custom_genome
for (i in seq_along(track_data)) {
    t <- track_data[[i]]
    if (length(t$gr) == 0) next
    keep <- which(overlapsAny(t$gr, clip_gr))
    track_data[[i]]$gr    <- t$gr[keep]
    track_data[[i]]$count <- t$count[keep]
}
nonempty <- Filter(function(t) length(t$count) > 0, track_data)
ymax <- if (length(nonempty) > 0) {
  max(unlist(lapply(nonempty, function(t) t$count)))
} else {
  1
}
cat("Shared ymax across tracks (within plotted region): ",
    round(ymax, 2), "\n", sep = "")

N <- length(track_data)
palette_cols <- if (N == 1) "#1f77b4" else hcl.colors(N, palette = "Dynamic")

# ---- plot ----
plot_kp <- function() {
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$ideogramheight  <- if (all_chr) 18 else 10
  pp$data1height     <- if (all_chr) max(80, 18 * N) else max(120, 22 * N)
  pp$data1inmargin   <- 14
  pp$data1outmargin  <- 24
  pp$topmargin       <- 30
  pp$bottommargin    <- 30
  pp$leftmargin      <- 0.10

  norm_tag <- if (!skip_dcs) "DCS-normalized" else "raw counts"
  title_txt <- if (all_chr) {
                 paste0("genome-wide   ", N, " track(s)   ", norm_tag,
                        "   (ymax=", round(ymax, 2), ")")
               } else {
                 paste0(target_chr, ":", start(zoom_gr), "-", end(zoom_gr),
                        "   ", N, " track(s)   ", norm_tag,
                        "   (ymax=", round(ymax, 2), ")")
               }

  kp <- if (all_chr) {
          plotKaryotype(genome = custom_genome, plot.type = 1,
                        plot.params = pp, main = title_txt)
        } else {
          plotKaryotype(genome = custom_genome, plot.type = 1, zoom = zoom_gr,
                        plot.params = pp, main = title_txt)
        }

  tick <- if (all_chr) 5e7 else 5e5
  kpAddBaseNumbers(kp, tick.dist = tick, add.units = TRUE, cex = 0.5,
                   tick.len = 1.5, minor.tick.len = 0.8)

  # chromosome bar
  kpRect(kp, data = custom_genome, y0 = 0, y1 = 1,
         col = "#eaebe4", border = "#888888", lwd = 0.5, data.panel = "ideogram")
  # faint broad centromere backdrop (drawn first, under the sharp regions)
  if (length(backdrop) > 0)
    kpRect(kp, data = backdrop, y0 = 0, y1 = 1,
           col = "#EAD7D7", border = NA, data.panel = "ideogram")
  # sharp region highlight (HOR / specific arrays) drawn on top
  if (length(regions) > 0)
    kpRect(kp, data = regions, y0 = 0, y1 = 1,
           col = "#C18787", border = "#C18787", data.panel = "ideogram")

  # One sub-track per bedgraph, top-to-bottom in arg order.
  # Leave a small gap between tracks so axis labels don't collide.
  slot_h   <- 1 / N
  gap_frac <- if (N > 1) 0.18 else 0   # 18% of each slot used as inter-track gap
  for (i in seq_len(N)) {
    t  <- track_data[[i]]
    r0 <- 1 - i * slot_h + (gap_frac * slot_h / 2)
    r1 <- 1 - (i - 1) * slot_h - (gap_frac * slot_h / 2)
    col <- palette_cols[i]

    # Only label the max value (single tick) — avoids the 0-tick of this track
    # overlapping with the max-tick of the track below it.
    kpAxis(kp, ymin = 0, ymax = ymax, data.panel = 1,
           r0 = r0, r1 = r1,
           tick.pos = ymax, labels = formatC(ymax, format = "g", digits = 3),
           cex = 0.45, col = "#666666")

    if (N > 1)
      kpAddLabels(kp, labels = t$name, r0 = r0, r1 = r1,
                  data.panel = 1, side = "right",
                  cex = 0.55, label.margin = 0.005)

    if (length(t$gr) > 0) {
      drawn <- resize(t$gr, width = draw_w, fix = "center")
      kpBars(kp, data = drawn, y0 = 0, y1 = t$count,
             ymin = 0, ymax = ymax,
             r0 = r0, r1 = r1,
             col = col, border = col, lwd = 0.2,
             data.panel = 1)
    }
  }
}

# ---- ensure output dir exists, then write ----
out_dir <- dirname(out_prefix)
if (nzchar(out_dir) && !dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

W <- if (all_chr) 14 else 12
H <- if (all_chr) max(6, 0.35 * N + 4) else max(4, 0.4 * N + 2)

pdf(paste0(out_prefix, ".pdf"), width = W, height = H); plot_kp(); dev.off()
tryCatch({
  png(paste0(out_prefix, ".png"), width = W * 150, height = H * 150,
      res = 150, type = "cairo")
  plot_kp(); dev.off()
}, error = function(e) message("PNG cairo failed: ", conditionMessage(e)))

cat("Wrote ", out_prefix, ".pdf / .png   tracks: ",
    paste(track_names, collapse = ", "), "\n", sep = "")
