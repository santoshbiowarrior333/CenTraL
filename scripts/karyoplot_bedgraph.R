#!/usr/bin/env Rscript
# Karyoplot of one OR many bedgraph tracks. Drop in one bedgraph to get a
# single-track plot; drop in several to get them stacked above each chromosome
# bar (each labelled, sharing the same y-axis for easy comparison).
#
# Same input pattern as before:
#   - chrom sizes file (UCSC chr<TAB>length)
#   - regions BED for centromere / HOR highlight (or NA)
#   - one or more bedgraphs (chr/start/end/count)
#
# Optional overlay BEDs (anywhere in args, drawn as full-height vertical
# tick rows at the top of the stack — independent of the y-axis):
#   --primer-fwd FILE   navy ticks at every primer site
#   --primer-rev FILE   dark-red ticks
#   --aso-plus FILE     dark-green ticks
#   --aso-minus FILE    dark-purple ticks
#
# Usage:
#   Rscript karyoplot_bedgraph.R <chrom.sizes> <regions> <backdrop> <out_prefix> \
#                                <target_chr> <zoom> <dcs_tsv> \
#                                <bedgraph1> [bedgraph2 ...] \
#                                [--primer-fwd F] [--primer-rev F] [--aso-plus F] [--aso-minus F]

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

# ---- args (parse named flags first, then positional) ----
raw_args   <- commandArgs(trailingOnly = TRUE)
primer_fwd <- NA_character_
primer_rev <- NA_character_
aso_plus   <- NA_character_
aso_minus  <- NA_character_
primer_fwd_name <- "primer-fwd"
primer_rev_name <- "primer-rev"
aso_plus_name   <- "aso-plus"
aso_minus_name  <- "aso-minus"
custom_names    <- NA_character_   # pipe-separated, positional for bedgraph tracks
custom_title    <- NA_character_   # full override of the plot heading
args <- character(0)
i <- 1
while (i <= length(raw_args)) {
  a <- raw_args[i]
  if      (a == "--primer-fwd")      { primer_fwd      <- raw_args[i + 1]; i <- i + 2 }
  else if (a == "--primer-rev")      { primer_rev      <- raw_args[i + 1]; i <- i + 2 }
  else if (a == "--aso-plus")        { aso_plus        <- raw_args[i + 1]; i <- i + 2 }
  else if (a == "--aso-minus")       { aso_minus       <- raw_args[i + 1]; i <- i + 2 }
  else if (a == "--primer-fwd-name") { primer_fwd_name <- raw_args[i + 1]; i <- i + 2 }
  else if (a == "--primer-rev-name") { primer_rev_name <- raw_args[i + 1]; i <- i + 2 }
  else if (a == "--aso-plus-name")   { aso_plus_name   <- raw_args[i + 1]; i <- i + 2 }
  else if (a == "--aso-minus-name")  { aso_minus_name  <- raw_args[i + 1]; i <- i + 2 }
  else if (a == "--names")           { custom_names    <- raw_args[i + 1]; i <- i + 2 }
  else if (a == "--title")           { custom_title    <- raw_args[i + 1]; i <- i + 2 }
  else                               { args <- c(args, a);                 i <- i + 1 }
}

if (length(args) < 8)
  stop("Usage: Rscript karyoplot_bedgraph.R <chrom.sizes> <regions> <backdrop> <out_prefix> <target_chr> <zoom> <dcs_tsv> <bedgraph1> [bedgraph2 ...]  [--primer-fwd F] [--primer-rev F] [--aso-plus F] [--aso-minus F]")

chrom_sizes_file <- args[1]
regions_file     <- args[2]
backdrop_file    <- args[3]
out_prefix       <- args[4]
target_chr       <- args[5]
zoom_arg         <- args[6]
dcs_tsv_file     <- args[7]
bedgraph_files   <- args[8:length(args)]

is_na_arg <- function(x) is.null(x) || length(x) == 0 || is.na(x) || tolower(x) %in% c("na", "", "-", "none")

# Two color versions per overlay:
#   *_RGB   full opacity — used for labels so the text stays readable
#   *_COL   alpha 0.5     — used for the ticks so dense regions don't form
#                            solid dark blocks (you keep density information)
PRIMER_FWD_RGB <- "#1e3a8a"
PRIMER_REV_RGB <- "#991b1b"
ASO_PLUS_RGB   <- "#166534"
ASO_MINUS_RGB  <- "#6b21a8"
PRIMER_FWD_COL <- adjustcolor(PRIMER_FWD_RGB, alpha.f = 0.5)
PRIMER_REV_COL <- adjustcolor(PRIMER_REV_RGB, alpha.f = 0.5)
ASO_PLUS_COL   <- adjustcolor(ASO_PLUS_RGB,   alpha.f = 0.5)
ASO_MINUS_COL  <- adjustcolor(ASO_MINUS_RGB,  alpha.f = 0.5)

read_bed3_loose <- function(path) {
  d <- read.table(path, header = FALSE, sep = "\t", quote = "",
                  stringsAsFactors = FALSE, comment.char = "#")
  if (ncol(d) < 3) stop("Overlay BED needs >= 3 columns: ", path)
  GRanges(seqnames = d[[1]], ranges = IRanges(start = d[[2]] + 1, end = d[[3]]))
}

# Split a single arg into (path, inline_name).
# Syntax: "path/to/file.bed||Display name with spaces"
# Returns a list with $path and $name (NA if no inline name was given).
split_pathname <- function(s) {
  if (is_na_arg(s)) return(list(path = NA_character_, name = NA_character_))
  parts <- strsplit(s, "||", fixed = TRUE)[[1]]
  list(
    path = parts[1],
    name = if (length(parts) >= 2 && nzchar(parts[2])) parts[2] else NA_character_
  )
}

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
  if (nrow(sizes) == 0) stop("target_chr '", target_chr, "' not in ", chrom_sizes_file)
}
custom_genome <- toGRanges(data.frame(chr = sizes$chr, start = 1, end = sizes$length))

# ---- optional DCS scale factor lookup ----
scale_lookup <- list()
if (!skip_dcs) {
  dcs <- read.table(dcs_tsv_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  if (!all(c("barcode", "scale_factor") %in% colnames(dcs)))
    stop("DCS TSV needs columns 'barcode' and 'scale_factor': ", dcs_tsv_file)
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
  # Each bedgraph arg can carry an inline display name: "path||Sample A".
  pn <- split_pathname(bedgraph_files[i])
  f  <- pn$path
  if (!file.exists(f)) stop("Bedgraph not found: ", f)
  if (!is.na(pn$name)) {
    bc <- pn$name
  } else {
    bc <- sub("\\.bedgraph$", "", sub("\\.startcount\\.bedgraph$", "", basename(f)))
    bc <- sub("\\.sites\\.bedgraph$", "", bc)
  }
  track_names[i] <- bc
  bg <- read.table(f, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("chr", "start", "end", "count"))
  if (!all_chr) bg <- bg[bg$chr == target_chr, ]
  bg <- bg[bg$chr %in% sizes$chr, ]
  if (nrow(bg) == 0) {
    cat("  ", bc, ": no entries on target — track will be blank\n", sep = "")
    track_data[[i]] <- list(gr = GRanges(), count = numeric(0), name = bc, sf = 1)
    next
  }
  sf <- if (bc %in% names(scale_lookup)) scale_lookup[[bc]] else 1
  scaled <- bg$count * sf
  cat("  ", bc, ": ", nrow(bg), " positions, count ",
      round(min(scaled), 2), "-", round(max(scaled), 2),
      "  (scale=", round(sf, 4), ")\n", sep = "")
  track_data[[i]] <- list(gr = toGRanges(bg[, 1:3]), count = scaled, name = bc, sf = sf)
}
nonempty <- Filter(function(t) length(t$count) > 0, track_data)
if (length(nonempty) == 0) stop("All tracks empty — nothing to plot.")

# Apply user-supplied bedgraph names if provided (pipe-separated, positional).
# Blank or missing slots fall back to whatever was parsed from the filename.
if (!is_na_arg(custom_names)) {
  parsed_names <- strsplit(custom_names, "|", fixed = TRUE)[[1]]
  for (i in seq_along(track_data)) {
    if (i <= length(parsed_names) && nzchar(parsed_names[i])) {
      track_data[[i]]$name <- parsed_names[i]
      track_names[i]       <- parsed_names[i]
    }
  }
  cat("Renamed tracks (from --names): ",
      paste(track_names, collapse = ", "), "\n", sep = "")
}

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
    if (length(backdrop) > 0)      { ref_gr <- backdrop; src <- "backdrop" }
    else if (length(regions) > 0)  { ref_gr <- regions;  src <- "regions"  }
    else                           { ref_gr <- unlist(GRangesList(lapply(nonempty, function(t) t$gr))); src <- "data" }
    if (length(ref_gr) == 0) {
      zoom_gr <- custom_genome
    } else {
      # Pick the largest cluster of regions instead of the full extent.
      # Otherwise an outlier HOR (e.g. the tiny 132 Mb island on chr2) drags
      # the zoom out across the whole chromosome and crushes the real signal.
      # Cluster = contiguous group of regions separated by gaps <= max_gap.
      # Largest = the cluster with the most total covered bases.
      dominant_cluster <- function(gr, max_gap = 5e6) {
        gr <- sort(gr)
        if (length(gr) <= 1) return(c(min(start(gr)), max(end(gr))))
        st <- start(gr); en <- end(gr)
        gaps <- st[-1] - en[-length(en)]
        ids <- cumsum(c(0, gaps > max_gap)) + 1
        widths <- tapply(en - st + 1, ids, sum)
        best <- as.integer(names(which.max(widths)))
        keep <- which(ids == best)
        c(min(st[keep]), max(en[keep]))
      }
      dc <- dominant_cluster(ref_gr)
      n_drop <- length(ref_gr) - sum(start(ref_gr) >= dc[1] & end(ref_gr) <= dc[2])
      if (n_drop > 0)
        cat("auto-zoom: dropped ", n_drop, " outlier region(s) from zoom calc\n", sep = "")
      z_start <- max(1,            dc[1] - pad)
      z_end   <- min(sizes$length, dc[2] + pad)
      zoom_gr <- toGRanges(data.frame(chr = target_chr, start = z_start, end = z_end))
      cat("auto-zoom from ", src, ": ", target_chr, ":", z_start, "-", z_end, "\n", sep = "")
    }
  } else {
    mm <- regmatches(zoom_arg, regexec("^([^:]+):(\\d+)-(\\d+)$", zoom_arg))[[1]]
    if (length(mm) != 4) stop("Bad zoom: ", zoom_arg)
    zoom_gr <- toGRanges(data.frame(chr = mm[2], start = as.numeric(mm[3]), end = as.numeric(mm[4])))
  }
}

# widen narrow bars for visibility at this zoom
zoomspan <- if (!is.null(zoom_gr)) end(zoom_gr) - start(zoom_gr) + 1 else sum(as.numeric(width(custom_genome)))
draw_w   <- max(round(zoomspan / 1500), 50)

# clip each track to the visible window so y-axis reflects what's shown
clip_gr <- if (!is.null(zoom_gr)) zoom_gr else custom_genome
for (i in seq_along(track_data)) {
  t <- track_data[[i]]
  if (length(t$gr) == 0) next
  keep <- which(overlapsAny(t$gr, clip_gr))
  track_data[[i]]$gr    <- t$gr[keep]
  track_data[[i]]$count <- t$count[keep]
}
nonempty <- Filter(function(t) length(t$count) > 0, track_data)
ymax <- if (length(nonempty) > 0) max(unlist(lapply(nonempty, function(t) t$count))) else 1
cat("Shared ymax across bedgraph tracks (within plotted region): ", round(ymax, 2), "\n", sep = "")

N <- length(track_data)
palette_cols <- if (N == 1) "#1f77b4" else hcl.colors(N, palette = "Dynamic")

# ---- overlay BEDs (loaded after zoom so we can clip them too) ----
load_overlay <- function(path) {
  if (is_na_arg(path)) return(NULL)
  if (!file.exists(path)) { message("Overlay not found: ", path); return(NULL) }
  gr <- read_bed3_loose(path)
  if (!all_chr) gr <- gr[seqnames(gr) == target_chr]
  if (!is.null(zoom_gr)) gr <- gr[overlapsAny(gr, zoom_gr)]
  gr
}
overlays <- list()
addov <- function(name, gr, col, label_col) {
  if (is.null(gr)) return(invisible())
  overlays[[length(overlays) + 1]] <<- list(
    name      = name,
    gr        = gr,
    col       = col,         # alpha 0.5 — for tick segments
    label_col = label_col    # full opacity — for the label text
  )
}
# Allow inline "path||Display name" on each --primer-/--aso- arg. The
# explicit --primer-fwd-name (etc.) flag wins; otherwise the inline name
# wins; otherwise the default ("primer-fwd" etc.) is used.
pf <- split_pathname(primer_fwd); pr <- split_pathname(primer_rev)
ap <- split_pathname(aso_plus);   am <- split_pathname(aso_minus)
pick_name <- function(explicit, inline, default) {
  if (!identical(explicit, default)) return(explicit)   # user set --..-name
  if (!is.na(inline)) return(inline)
  default
}
primer_fwd_name <- pick_name(primer_fwd_name, pf$name, "primer-fwd")
primer_rev_name <- pick_name(primer_rev_name, pr$name, "primer-rev")
aso_plus_name   <- pick_name(aso_plus_name,   ap$name, "aso-plus")
aso_minus_name  <- pick_name(aso_minus_name,  am$name, "aso-minus")

addov(primer_fwd_name, load_overlay(pf$path), PRIMER_FWD_COL, PRIMER_FWD_RGB)
addov(primer_rev_name, load_overlay(pr$path), PRIMER_REV_COL, PRIMER_REV_RGB)
addov(aso_plus_name,   load_overlay(ap$path), ASO_PLUS_COL,   ASO_PLUS_RGB)
addov(aso_minus_name,  load_overlay(am$path), ASO_MINUS_COL,  ASO_MINUS_RGB)
n_overlay <- length(overlays)
if (n_overlay > 0) cat("overlay BED rows: ", n_overlay, "\n", sep = "")

# ---- plot ----
plot_kp <- function() {
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$ideogramheight  <- if (all_chr) 18 else 10
  rows_total         <- N + n_overlay
  pp$data1height     <- if (all_chr) max(80, 18 * rows_total) else max(120, 22 * rows_total)
  pp$data1inmargin   <- 14
  pp$data1outmargin  <- 24
  pp$topmargin       <- 30
  pp$bottommargin    <- 30
  pp$leftmargin      <- 0.10
  # Give right-side labels (track / overlay names with custom text) enough
  # canvas so long names like "Chr17_forward(+RT)+ASO treated" don't clip.
  # Scaled by the longest label among the tracks + overlay rows we're drawing.
  all_label_names <- c(track_names,
                       vapply(overlays, function(o) o$name, character(1)))
  longest_label <- if (length(all_label_names) > 0) max(nchar(all_label_names)) else 0
  # ~0.007 of canvas width per character at cex 0.55, capped between 0.18 and 0.40
  pp$rightmargin <- max(0.18, min(0.40, 0.06 + 0.007 * longest_label))

  # Simple two-state label: if any track was actually scaled by DCS we call it
  # DCS-normalized; otherwise raw counts. Same logic whether the TSV was NA
  # or just had all-NA / all-1.0 rows.
  any_scaled <- any(vapply(track_data,
                           function(t) !is.null(t$sf) && t$sf != 1,
                           logical(1)))
  norm_tag <- if (any_scaled) "DCS-normalized" else "raw counts"

  default_title <- if (all_chr) {
                     paste0("genome-wide   ", N, " track(s)   ", norm_tag,
                            "   (ymax=", round(ymax, 2), ")")
                   } else {
                     paste0(target_chr, ":", start(zoom_gr), "-", end(zoom_gr),
                            "   ", N, " track(s)   ", norm_tag,
                            "   (ymax=", round(ymax, 2), ")")
                   }
  title_txt <- if (is_na_arg(custom_title)) default_title else custom_title

  kp <- if (all_chr) {
          plotKaryotype(genome = custom_genome, plot.type = 1, plot.params = pp, main = title_txt)
        } else {
          plotKaryotype(genome = custom_genome, plot.type = 1, zoom = zoom_gr,
                        plot.params = pp, main = title_txt)
        }

  tick <- if (all_chr) 5e7 else 5e5
  kpAddBaseNumbers(kp, tick.dist = tick, add.units = TRUE, cex = 0.5,
                   tick.len = 1.5, minor.tick.len = 0.8)

  # chrom bar + highlights on the ideogram
  kpRect(kp, data = custom_genome, y0 = 0, y1 = 1,
         col = "#eaebe4", border = "#888888", lwd = 0.5, data.panel = "ideogram")
  if (length(backdrop) > 0)
    kpRect(kp, data = backdrop, y0 = 0, y1 = 1,
           col = "#EAD7D7", border = NA, data.panel = "ideogram")
  if (length(regions) > 0)
    kpRect(kp, data = regions, y0 = 0, y1 = 1,
           col = "#C18787", border = "#C18787", data.panel = "ideogram")

  total_slots <- n_overlay + N
  slot_h      <- 1 / total_slots
  gap_frac    <- if (total_slots > 1) 0.18 else 0
  band <- function(k) list(
    r0 = 1 - k * slot_h + (gap_frac * slot_h / 2),
    r1 = 1 - (k - 1) * slot_h - (gap_frac * slot_h / 2)
  )

  # Overlay rows on top — short vertical ticks rising from a baseline at
  # the bottom of the band, same visual idiom as the bedgraph bar tracks
  # below. Baseline drawn first (under the ticks), so ticks read as peaks
  # standing on the line, not bars bisected by it.
  TICK_HEIGHT <- 0.35   # fraction of the row band the tick rises to
  for (j in seq_len(n_overlay)) {
    ov <- overlays[[j]]; b <- band(j)
    # Baseline at the bottom of the band — ticks will rise from here.
    kpAbline(kp, h = 0, r0 = b$r0, r1 = b$r1, col = "#dddddd", lwd = 0.4)
    if (length(ov$gr) > 0) {
      kpSegments(kp,
                 chr = as.character(seqnames(ov$gr)),
                 x0 = start(ov$gr), x1 = start(ov$gr),
                 y0 = 0, y1 = TICK_HEIGHT,
                 r0 = b$r0, r1 = b$r1,
                 col = ov$col, lwd = 0.6, data.panel = 1)
    }
    kpAddLabels(kp,
                labels = sprintf("%s (n=%d)", ov$name, length(ov$gr)),
                r0 = b$r0, r1 = b$r1,
                data.panel = 1, side = "right",
                cex = 0.55, col = ov$label_col, label.margin = 0.005)
  }

  # Bedgraph rows below — bars with shared ymax
  for (i in seq_len(N)) {
    t   <- track_data[[i]]
    b   <- band(n_overlay + i)
    r0  <- b$r0; r1 <- b$r1
    col <- palette_cols[i]

    # Label both 0 and ymax. Slightly smaller cex so the 0 of this track and
    # the ymax of the track below stay legible even when bands are tight.
    kpAxis(kp, ymin = 0, ymax = ymax, data.panel = 1,
           r0 = r0, r1 = r1,
           tick.pos = c(0, ymax),
           labels = c("0", formatC(ymax, format = "g", digits = 3)),
           cex = 0.4, col = "#666666")

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

W <- if (all_chr) 16 else 14
H <- if (all_chr) max(6, 0.35 * (N + n_overlay) + 4) else max(4, 0.4 * (N + n_overlay) + 2)

pdf(paste0(out_prefix, ".pdf"), width = W, height = H); plot_kp(); dev.off()
tryCatch({
  png(paste0(out_prefix, ".png"), width = W * 150, height = H * 150,
      res = 150, type = "cairo")
  plot_kp(); dev.off()
}, error = function(e) message("PNG cairo failed: ", conditionMessage(e)))

cat("Wrote ", out_prefix, ".pdf / .png   tracks: ",
    paste(track_names, collapse = ", "),
    if (n_overlay > 0) paste0("   overlays: ", n_overlay) else "",
    "\n", sep = "")
