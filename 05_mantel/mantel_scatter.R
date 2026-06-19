#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Mantel tests (mtDNA vs MSY pairwise distances) + scatter plots
#
# Paper: "Sex-biased demographic history of South American camelids revealed
#         by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold"
#
# Tests performed:
#   (i)  Domestic species: L. glama + V. pacos
#   (ii) Wild species: L. g. cacsilensis + L. g. guanicoe + V. v. mensalis + V. v. vicugna
#
# Outputs:
#   mantel_results.txt              — printed Mantel test summaries
#   mantel_summary.csv              — r, p-value per block
#   scatter_domestic.tiff/.pdf      — panel A
#   scatter_wild.tiff/.pdf          — panel B
#
# Dependencies: ape, vegan, ggplot2
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ape); library(vegan); library(ggplot2)
})

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
MT_FASTA  <- "path/to/mtdna_wo_outgroup_aln.fasta"
MSY_FASTA <- "path/to/csa_male_msy_wo_outgroup.fasta"
META_TSV  <- "data/metadata.tsv"   # columns: sample, group, status
OUT_DIR   <- "results/mantel"
FIG_DIR   <- "results/figures"
N_PERM    <- 9999
DOMESTIC_GROUPS <- c("L_glama", "V_pacos")
WILD_GROUPS     <- c("L_g_cacsilensis", "L_g_guanicoe", "V_v_mensalis", "V_v_vicugna")
# ─────────────────────────────────────────────────────────────────────────────

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1) Load data ──────────────────────────────────────────────────────────────
cat("Reading alignments...\n")
mt  <- read.dna(MT_FASTA,  format = "fasta")
msy <- read.dna(MSY_FASTA, format = "fasta")

meta <- read.table(META_TSV, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
names(meta) <- tolower(names(meta))
meta$sample <- trimws(meta$sample)
meta$group  <- trimws(meta$group)

ids <- intersect(rownames(mt), rownames(msy))
ids <- intersect(ids, meta$sample)
mt  <- mt[ids, ]; msy <- msy[ids, ]
meta <- meta[meta$sample %in% ids, ]
meta <- meta[match(ids, meta$sample), ]
stopifnot(all(meta$sample == ids))
cat(sprintf("Common samples: %d\n", length(ids)))

# ── 2) Compute pairwise distances ─────────────────────────────────────────────
cat("Computing K80 distances...\n")
d_mt  <- dist.dna(mt,  model = "raw", pairwise.deletion = TRUE)
d_msy <- dist.dna(msy, model = "raw", pairwise.deletion = TRUE)

# ── 3) Mantel test helper ─────────────────────────────────────────────────────
run_block <- function(group_vector, label) {
  cat(sprintf("\n=== %s ===\n", label))
  sub_ids <- meta$sample[meta$group %in% group_vector]
  if (length(sub_ids) < 4) stop("Too few samples in block: ", label)

  d_mt_sub  <- as.dist(as.matrix(d_mt)[sub_ids, sub_ids])
  d_msy_sub <- as.dist(as.matrix(d_msy)[sub_ids, sub_ids])

  set.seed(1)
  m <- mantel(d_mt_sub, d_msy_sub, method = "spearman", permutations = N_PERM)
  print(m)

  # Upper-triangle pairs for scatter
  M_mt  <- as.matrix(d_mt_sub); M_msy <- as.matrix(d_msy_sub)
  mask  <- upper.tri(M_mt) & is.finite(M_mt) & is.finite(M_msy)
  list(
    mantel = m,
    df     = data.frame(mtDNA = M_mt[mask], MSY = M_msy[mask])
  )
}

blk_dom  <- run_block(DOMESTIC_GROUPS, "Domestic")
blk_wild <- run_block(WILD_GROUPS,     "Wild")

# Summary table
mantel_summary <- data.frame(
  block   = c("Domestic", "Wild"),
  r       = c(blk_dom$mantel$statistic, blk_wild$mantel$statistic),
  p_value = c(blk_dom$mantel$signif,    blk_wild$mantel$signif)
)
cat("\n=== Mantel summary ===\n"); print(mantel_summary)
write.csv(mantel_summary, file.path(OUT_DIR, "mantel_summary.csv"), row.names = FALSE)

# Full output
sink(file.path(OUT_DIR, "mantel_results.txt"))
cat("=== Domestic ===\n"); print(blk_dom$mantel)
cat("\n=== Wild ===\n");    print(blk_wild$mantel)
sink()
cat("Saved: mantel_summary.csv, mantel_results.txt\n")

# ── 4) Scatter plots ──────────────────────────────────────────────────────────
df_dom  <- blk_dom$df;  m_dom  <- blk_dom$mantel
df_wild <- blk_wild$df; m_wild <- blk_wild$mantel

xmax <- max(c(df_dom$mtDNA, df_wild$mtDNA), na.rm = TRUE)
xlim_all <- c(0, xmax * 1.03)

pretty_breaks_max <- function(lims, max_ticks = 6) {
  rng <- diff(lims)
  if (!is.finite(rng) || rng <= 0) return(list(lims = lims, breaks = lims[1]))
  raw_step <- rng / (max_ticks - 1)
  pow10    <- 10^floor(log10(raw_step))
  step     <- (c(1, 2, 5, 10) * pow10)[which.min(abs(c(1, 2, 5, 10) * pow10 - raw_step))]
  lo <- floor(lims[1] / step) * step; hi <- ceiling(lims[2] / step) * step
  brks <- seq(lo, hi, by = step)
  while (length(brks) > max_ticks) {
    step <- step * 2
    lo <- floor(lims[1] / step) * step; hi <- ceiling(lims[2] / step) * step
    brks <- seq(lo, hi, by = step)
  }
  list(lims = c(lo, hi), breaks = brks)
}

yb       <- pretty_breaks_max(range(c(df_dom$MSY, df_wild$MSY), na.rm = TRUE))
ylim_all <- yb$lims; ybreaks <- yb$breaks

make_scatter <- function(df, m, title) {
  ann <- paste0("Mantel Spearman r = ", round(m$statistic, 3),
                "\np = ", signif(m$signif, 3))
  ggplot(df, aes(x = mtDNA, y = MSY)) +
    geom_point(alpha = 0.35, size = 1.4) +
    geom_smooth(method = "lm", color = "black", fill = "grey80",
                se = TRUE, linewidth = 0.6) +
    coord_cartesian(xlim = xlim_all, ylim = ylim_all, expand = TRUE) +
    scale_y_continuous(breaks = ybreaks) +
    annotate("text",
             x = xlim_all[1] + diff(xlim_all) * 0.05, y = ylim_all[2],
             label = ann, hjust = 0, vjust = 1.2, size = 3.5) +
    labs(title = title, x = "mtDNA genetic distance", y = "MSY genetic distance") +
    theme_classic(base_size = 12)
}

p_dom  <- make_scatter(df_dom,  m_dom,  "A — Domestic")
p_wild <- make_scatter(df_wild, m_wild, "B — Wild")

for (nm in c("domestic", "wild")) {
  p <- if (nm == "domestic") p_dom else p_wild
  ggsave(file.path(FIG_DIR, paste0("scatter_", nm, ".pdf")),
         p, width = 4, height = 4, device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0("scatter_", nm, ".tiff")),
         p, width = 4, height = 4, device = "tiff",
         dpi = 600, compression = "lzw")
}
cat("Saved: scatter_domestic.pdf/.tiff, scatter_wild.pdf/.tiff\n")
