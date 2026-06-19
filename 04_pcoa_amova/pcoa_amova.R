#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# PCoA (mtDNA and MSY) + AMOVA
#
# Paper: "Sex-biased demographic history of South American camelids revealed
#         by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold"
#
# Inputs:
#   MT_FASTA   — mtDNA alignment without outgroup (44 individuals, 16,643 bp)
#   MSY_FASTA  — MSY alignment without outgroup  (44 individuals, 20,577,031 bp)
#   META_TSV   — metadata with columns: sample, group, status
#
# Outputs:
#   PCoA_mtDNA_MSY.pdf/.tiff     — two-panel PCoA figure
#   AMOVA_summary.csv            — PhiST and p-values (4 AMOVA tests)
#   AMOVA_results.txt            — full AMOVA output
#
# Dependencies: ape, pegas, ggplot2, dplyr, patchwork
# Note: reading the 20.6 Mb MSY alignment is memory-intensive (~4–8 GB RAM)
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ape); library(pegas); library(ggplot2); library(dplyr); library(patchwork)
})

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
MT_FASTA  <- "path/to/mtdna_wo_outgroup_aln.fasta"
MSY_FASTA <- "path/to/csa_male_msy_wo_outgroup.fasta"
META_TSV  <- "data/metadata.tsv"   # columns: sample, group, status (wild/domestic)
OUT_DIR   <- "results/pcoa_amova"
FIG_DIR   <- "results/figures"
N_PERM    <- 999
DIST_MODEL <- "K80"
# ─────────────────────────────────────────────────────────────────────────────

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

PAL <- c(
  L_g_cacsilensis = "#000000", L_g_guanicoe = "#E69F00",
  L_glama         = "#009E73", V_pacos       = "#56B4E9",
  V_v_mensalis    = "#0072B2", V_v_vicugna   = "#CC79A7"
)
SPECIES_LABELS <- list(
  L_g_cacsilensis = expression(italic(L.~g.~cacsilensis)),
  L_g_guanicoe    = expression(italic(L.~g.~guanicoe)),
  L_glama         = expression(italic(L.~glama)),
  V_v_mensalis    = expression(italic(V.~v.~mensalis)),
  V_v_vicugna     = expression(italic(V.~v.~vicugna)),
  V_pacos         = expression(italic(V.~pacos))
)

# ── 1) Load data ──────────────────────────────────────────────────────────────
cat("Reading alignments...\n")
mt  <- read.dna(MT_FASTA,  format = "fasta")
msy <- read.dna(MSY_FASTA, format = "fasta")   # large — may take several minutes

meta <- read.table(META_TSV, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
meta$sample <- trimws(meta$sample)
meta$group  <- trimws(meta$group)
meta$status <- trimws(meta$status)
stopifnot(all(c("sample", "group", "status") %in% names(meta)))

meta_mt  <- meta[match(rownames(mt),  meta$sample), ]
meta_msy <- meta[match(rownames(msy), meta$sample), ]
stopifnot(all(meta_mt$sample  == rownames(mt)),
          all(meta_msy$sample == rownames(msy)))

# ── 2) Distance matrices and PCoA ─────────────────────────────────────────────
cat("Computing K80 distances...\n")
mt_dist  <- dist.dna(mt,  model = DIST_MODEL)
msy_dist <- dist.dna(msy, model = DIST_MODEL)

mt_pcoa  <- pcoa(mt_dist)
msy_pcoa <- pcoa(msy_dist)

mt_pc1  <- mt_pcoa$values$Relative_eig[1]  * 100
mt_pc2  <- mt_pcoa$values$Relative_eig[2]  * 100
msy_pc1 <- msy_pcoa$values$Relative_eig[1] * 100
msy_pc2 <- msy_pcoa$values$Relative_eig[2] * 100

cat(sprintf("mtDNA  PC1 = %.1f%%  PC2 = %.1f%%\n", mt_pc1,  mt_pc2))
cat(sprintf("MSY    PC1 = %.1f%%  PC2 = %.1f%%\n", msy_pc1, msy_pc2))

make_df <- function(meta, pcoa_obj) {
  data.frame(
    sample = meta$sample,
    group  = factor(meta$group,  levels = names(PAL)),
    status = factor(meta$status, levels = c("wild", "domestic")),
    PC1    = pcoa_obj$vectors[, 1],
    PC2    = pcoa_obj$vectors[, 2],
    stringsAsFactors = FALSE
  )
}
df_mt  <- make_df(meta_mt,  mt_pcoa)
df_msy <- make_df(meta_msy, msy_pcoa)

# ── 3) PCoA figure ────────────────────────────────────────────────────────────
theme_paper <- theme_classic(base_size = 12) +
  theme(legend.position = "bottom", axis.title = element_text(size = 12),
        axis.text = element_text(size = 11))

plot_pcoa <- function(df, pc1, pc2, title) {
  ggplot(df, aes(PC1, PC2, color = group, shape = status)) +
    geom_point(size = 4.0, alpha = 0.85) +
    scale_color_manual(values = PAL, breaks = names(PAL),
                       labels = SPECIES_LABELS, name = "Population") +
    scale_shape_manual(values = c(wild = 16, domestic = 17), name = "Status") +
    labs(title = title,
         x = sprintf("PC1 (%.1f%%)", pc1),
         y = sprintf("PC2 (%.1f%%)", pc2)) +
    theme_paper
}

p_mt  <- plot_pcoa(df_mt,  mt_pc1,  mt_pc2,  "PCoA — mtDNA")
p_msy <- plot_pcoa(df_msy, msy_pc1, msy_pc2, "PCoA — MSY scaffold")

p_combined <- (p_mt | plot_spacer() | p_msy) +
  plot_layout(guides = "collect", widths = c(1, 0.08, 1)) &
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "PCoA_mtDNA_MSY.pdf"),
       p_combined, width = 10, height = 4.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "PCoA_mtDNA_MSY.tiff"),
       p_combined, width = 10, height = 4.5, device = "tiff",
       dpi = 600, compression = "lzw")
cat("Saved: PCoA_mtDNA_MSY.pdf / .tiff\n")

# ── 4) AMOVA ──────────────────────────────────────────────────────────────────
cat(sprintf("\nRunning AMOVA (nperm = %d)...\n", N_PERM))
grp_mt  <- factor(meta_mt$group)
grp_msy <- factor(meta_msy$group)
stat_mt  <- factor(meta_mt$status)
stat_msy <- factor(meta_msy$status)

mt_amova_grp  <- pegas::amova(mt_dist  ~ grp_mt,  nperm = N_PERM)
msy_amova_grp <- pegas::amova(msy_dist ~ grp_msy, nperm = N_PERM)
mt_amova_dom  <- pegas::amova(mt_dist  ~ stat_mt,  nperm = N_PERM)
msy_amova_dom <- pegas::amova(msy_dist ~ stat_msy, nperm = N_PERM)

amova_row <- function(obj, marker, model) {
  data.frame(
    marker  = marker, grouping = model,
    PhiST   = as.numeric(obj$PhiST),
    p_value = as.numeric(obj$P.value)
  )
}
amova_tbl <- bind_rows(
  amova_row(mt_amova_grp,  "mtDNA", "6 populations"),
  amova_row(msy_amova_grp, "MSY",   "6 populations"),
  amova_row(mt_amova_dom,  "mtDNA", "wild vs domestic"),
  amova_row(msy_amova_dom, "MSY",   "wild vs domestic")
)
cat("\n=== AMOVA summary ===\n"); print(amova_tbl)
write.csv(amova_tbl, file.path(OUT_DIR, "AMOVA_summary.csv"), row.names = FALSE)

sink(file.path(OUT_DIR, "AMOVA_results.txt"))
cat("AMOVA mtDNA (6 groups)\n");    print(mt_amova_grp)
cat("\nAMOVA MSY (6 groups)\n");    print(msy_amova_grp)
cat("\nAMOVA mtDNA (wild vs domestic)\n"); print(mt_amova_dom)
cat("\nAMOVA MSY (wild vs domestic)\n");   print(msy_amova_dom)
sink()
cat("Saved: AMOVA_summary.csv, AMOVA_results.txt\n")
