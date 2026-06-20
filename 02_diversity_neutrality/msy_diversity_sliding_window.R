#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# MSY nucleotide diversity, Tajima's D, FST, and sliding-window analysis
#
# Paper: "Sex-biased introgression during South American camelid domestication revealed
#         by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold"
#
# Inputs:
#   MSY_FASTA  — complete MSY scaffold alignment without outgroup (44 individuals)
#   GROUPS_TSV — sample-to-population mapping (see data/groups_example.tsv)
#   MT_CSV     — π per site per population for mtDNA (from a prior mtDNA analysis)
#
# Outputs (saved to OUT_DIR):
#   msy_popgenome_summary.csv          — π, θW, Tajima's D per population
#   msy_fst_pairwise.csv               — 6×6 pairwise nucleotide FST matrix
#   msy_ratio_pi_mt_over_msy.csv       — π_mtDNA / π_MSY ratio per population
#   msy_tajimaD_windows_kb.csv         — Tajima's D per 50-kb window
#   msy_pi_windows_kb.csv              — π per 50-kb window
#   msy_sliding_window_scan.pdf/.tiff  — combined sliding-window figure
#
# Dependencies: PopGenome, ggplot2, patchwork
# Runtime: 20–60 min (readData) + 10–20 min (sliding windows) on a desktop
# Memory : 4–8 GB RAM
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(PopGenome)
  library(ggplot2)
})

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
MSY_FASTA  <- "path/to/csa_male_msy_wo_outgroup.fasta"   # 44-seq MSY alignment
GROUPS_TSV <- "data/groups.tsv"
MT_CSV     <- "path/to/mtdna_pi_per_pop.csv"   # columns: group, pi_per_site
OUT_DIR    <- "results/diversity"
FIG_DIR    <- "results/figures"
# ─────────────────────────────────────────────────────────────────────────────

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(MSY_FASTA), file.exists(GROUPS_TSV))

# Population order and display labels
POP_ORDER <- c("L_g_cacsilensis", "L_g_guanicoe", "L_glama",
               "V_v_mensalis",    "V_v_vicugna",  "V_pacos")
POP_LABELS <- c(
  L_g_cacsilensis = "L. g. cacsilensis",
  L_g_guanicoe    = "L. g. guanicoe",
  L_glama         = "L. glama",
  V_v_mensalis    = "V. v. mensalis",
  V_v_vicugna     = "V. v. vicugna",
  V_pacos         = "V. pacos"
)

# Wong (2011) colorblind-safe palette
PAL <- c(
  L_g_cacsilensis = "#000000",
  L_g_guanicoe    = "#E69F00",
  L_glama         = "#009E73",
  V_v_mensalis    = "#0072B2",
  V_v_vicugna     = "#CC79A7",
  V_pacos         = "#56B4E9"
)
LTY_MAP <- c(
  L_g_cacsilensis = "solid",  L_g_guanicoe = "dashed",
  L_glama         = "solid",  V_v_mensalis  = "solid",
  V_v_vicugna     = "dashed", V_pacos       = "solid"
)

# ── 1) Load population groups ─────────────────────────────────────────────────
gt <- read.table(GROUPS_TSV, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, quote = "", comment.char = "")
names(gt) <- tolower(names(gt))
gt$sample <- trimws(gt$sample)
gt$group  <- trimws(gt$group)
gt <- gt[gt$sample != "" & gt$group != "", ]
groups    <- split(gt$sample, gt$group)
groups    <- groups[POP_ORDER[POP_ORDER %in% names(groups)]]
pop_names <- names(groups)
cat("Populations:", paste(pop_names, collapse = ", "), "\n")

# ── 2) Read alignment ─────────────────────────────────────────────────────────
# PopGenome requires the FASTA to be alone in a directory; create a temp symlink
tmp_dir <- file.path(tempdir(), paste0("pg_msy_", as.integer(Sys.time())))
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
link_path <- file.path(tmp_dir, basename(MSY_FASTA))
if (!file.exists(link_path)) file.symlink(MSY_FASTA, link_path)
if (!file.exists(link_path)) file.copy(MSY_FASTA, tmp_dir, overwrite = TRUE)

cat(sprintf("Reading MSY alignment (%s) — may take 20–60 min...\n",
            format(file.size(MSY_FASTA), big.mark = ",")))
t0 <- proc.time()
GENOME <- readData(tmp_dir, format = "FASTA")
cat(sprintf("  Done in %.1f min\n", (proc.time() - t0)[3] / 60))

# Validate sample IDs
all_ids <- names(GENOME@region.data@n.nucleotides[[1]])
missing <- setdiff(unlist(groups), all_ids)
if (length(missing) > 0) stop("IDs in groups.tsv not found in FASTA: ",
                               paste(missing, collapse = ", "))
GENOME <- set.populations(GENOME, groups, diploid = FALSE)

# ── 3) Global diversity and neutrality statistics ─────────────────────────────
cat("Computing diversity.stats...\n")
GENOME <- diversity.stats(GENOME)
cat("Computing neutrality.stats...\n")
GENOME <- neutrality.stats(GENOME)

L <- min(unique(as.numeric(GENOME@region.data@n.nucleotides[[1]])))
cat(sprintf("Alignment length used for normalisation: %s bp\n",
            format(L, big.mark = ",")))

res <- data.frame(
  group           = pop_names,
  n               = sapply(groups, length),
  pi              = as.numeric(GENOME@nuc.diversity.within),
  thetaW          = as.numeric(GENOME@theta_Watterson),
  TajimaD         = as.numeric(GENOME@Tajima.D),
  stringsAsFactors = FALSE
)
res$pi_per_site     <- signif(res$pi / L, 4)
res$thetaW_per_site <- signif(res$thetaW / L, 4)
res$TajimaD         <- round(res$TajimaD, 4)
res$group <- factor(res$group, levels = POP_ORDER)
res <- res[order(res$group), ]; res$group <- as.character(res$group)

cat("\n=== Nucleotide diversity — MSY scaffold ===\n")
print(res[, c("group", "n", "pi_per_site", "thetaW_per_site", "TajimaD")])
write.csv(res, file.path(OUT_DIR, "msy_popgenome_summary.csv"), row.names = FALSE)

# ── 4) Pairwise FST ───────────────────────────────────────────────────────────
cat("\nComputing pairwise FST...\n")
GENOME <- F_ST.stats(GENOME, mode = "nucleotide")
fst_mat <- GENOME@nucleotide.F_ST
rownames(fst_mat) <- colnames(fst_mat) <- pop_names
cat("\n=== Pairwise FST ===\n"); print(round(fst_mat, 4))
write.csv(as.data.frame(fst_mat),
          file.path(OUT_DIR, "msy_fst_pairwise.csv"))

# ── 5) π_mtDNA / π_MSY ratio ─────────────────────────────────────────────────
if (file.exists(MT_CSV)) {
  mt <- read.csv(MT_CSV, stringsAsFactors = FALSE)
  ratio_df <- merge(mt[, c("group", "pi_per_site")],
                    res[, c("group", "pi_per_site")],
                    by = "group", suffixes = c("_mt", "_msy"))
  ratio_df$ratio_pi_mt_msy <- ratio_df$pi_per_site_mt / ratio_df$pi_per_site_msy
  ratio_df$group <- factor(ratio_df$group, levels = POP_ORDER)
  ratio_df <- ratio_df[order(ratio_df$group), ]
  cat("\n=== π_mtDNA / π_MSY ratio ===\n")
  print(ratio_df[, c("group", "pi_per_site_mt", "pi_per_site_msy", "ratio_pi_mt_msy")])
  write.csv(ratio_df, file.path(OUT_DIR, "msy_ratio_pi_mt_over_msy.csv"),
            row.names = FALSE)
} else {
  cat("MT_CSV not found — skipping ratio computation.\n")
}

# ── 6) Sliding-window analysis (50 kb / 25 kb) ───────────────────────────────
cat("\nComputing sliding windows (50 kb / 25 kb) — may take 10–20 min...\n")
t1 <- proc.time()
GENOMEw <- sliding.window.transform(GENOME, width = 50000, jump = 25000, type = 2)
GENOMEw <- diversity.stats(GENOMEw)
GENOMEw <- neutrality.stats(GENOMEw)
cat(sprintf("  Done in %.1f min\n", (proc.time() - t1)[3] / 60))

Dmat  <- GENOMEw@Tajima.D;           Dmat[!is.finite(Dmat)] <- NA
Pimat <- GENOMEw@nuc.diversity.within; Pimat[!is.finite(Pimat)] <- NA
colnames(Dmat) <- colnames(Pimat) <- pop_names
Pimat_site <- Pimat / 50000

parse_mid_kb <- function(lbl) {
  nums <- as.numeric(regmatches(lbl, gregexpr("[0-9]+", lbl))[[1]])
  (nums[1] + nums[2]) / 2 / 1000
}
midpoints_kb <- sapply(rownames(Dmat), parse_mid_kb)

write.csv(data.frame(pos_kb = midpoints_kb, as.data.frame(Dmat)),
          file.path(OUT_DIR, "msy_tajimaD_windows_kb.csv"), row.names = FALSE)
write.csv(data.frame(pos_kb = midpoints_kb, as.data.frame(Pimat_site)),
          file.path(OUT_DIR, "msy_pi_windows_kb.csv"), row.names = FALSE)
cat("Saved: msy_tajimaD_windows_kb.csv, msy_pi_windows_kb.csv\n")

# ── 7) Sliding-window figure ──────────────────────────────────────────────────
n_win <- length(midpoints_kb)
long_df <- function(mat, stat) {
  do.call(rbind, lapply(pop_names, function(p) {
    data.frame(pos_kb = midpoints_kb, value = mat[, p],
               population = p, stat = stat)
  }))
}
df_pi <- long_df(Pimat_site, "pi"); df_D <- long_df(Dmat, "D")
df_pi$population <- factor(df_pi$population, levels = POP_ORDER)
df_D$population  <- factor(df_D$population,  levels = POP_ORDER)

make_panel <- function(df, ylab_expr, reflines = NULL) {
  p <- ggplot(df, aes(x = pos_kb, y = value,
                      color = population, linetype = population)) +
    geom_line(linewidth = 0.65, alpha = 0.85, na.rm = TRUE)
  if (!is.null(reflines))
    for (rl in reflines)
      p <- p + geom_hline(yintercept = rl$y, linetype = rl$lty,
                           color = rl$col, linewidth = 0.4)
  p +
    scale_color_manual(values = PAL, labels = POP_LABELS, name = NULL) +
    scale_linetype_manual(values = LTY_MAP, labels = POP_LABELS, name = NULL) +
    scale_x_continuous(
      name   = "Position along MSY scaffold (Mb)",
      breaks = seq(0, max(df$pos_kb, na.rm = TRUE) + 500, by = 2000),
      labels = function(x) sprintf("%.0f", x / 1000),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(y = ylab_expr) +
    theme_classic(base_size = 11) +
    theme(axis.text = element_text(size = 9),
          panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
          legend.position = "top", legend.direction = "horizontal",
          legend.text = element_text(size = 9, face = "italic"),
          legend.key.width = unit(1.4, "cm"))
}

p_pi <- make_panel(df_pi, expression(pi[MSY]~(per~site)),
                   list(list(y = 0, lty = "longdash", col = "gray70")))
p_D  <- make_panel(df_D, "Tajima's D",
                   list(list(y =  0, lty = "longdash", col = "gray70"),
                        list(y =  2, lty = "dotted",   col = "gray60"),
                        list(y = -2, lty = "dotted",   col = "gray60")))

combined <- tryCatch({
  library(patchwork)
  (p_pi / p_D) + plot_layout(guides = "collect") &
    theme(legend.position = "top")
}, error = function(e) {
  library(gridExtra)
  grid.arrange(p_pi, p_D + theme(legend.position = "none"), ncol = 1)
})

ggsave(file.path(FIG_DIR, "msy_sliding_window_scan.pdf"),
       combined, width = 9, height = 6.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "msy_sliding_window_scan.tiff"),
       combined, width = 9, height = 6.5, device = "tiff",
       dpi = 300, compression = "lzw")
cat(sprintf("Saved: msy_sliding_window_scan.pdf / .tiff\n"))
cat(sprintf("\nTotal runtime: %.1f min\n", (proc.time() - t0)[3] / 60))
