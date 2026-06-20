#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Haplogroup assignment from tree topology (MSY and mtDNA)
#
# Paper: "Sex-biased introgression during South American camelid domestication revealed
#         by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold"
#
# Method:
#   1. Read rooted consensus trees (IQ-TREE .contree)
#   2. Find the largest MRCA clade containing ≥80% Vicugna+alpaca tips → H-Vicugna
#   3. All remaining ingroup tips → H-Lama
#   4. Individuals topologically in the "wrong" clade → Discordant
#
# Inputs:
#   MSY_TREE   — rooted MSY consensus tree (.contree)
#   MT_TREE    — rooted mtDNA consensus tree (.contree)
#   META_TSV   — metadata: columns sample, group, status
#
# Outputs:
#   haplogroup_tree_assignment.csv   — per-individual haplogroup assignment
#   haplogroup_discordance_table.csv — summary by species and marker
#   haplogroup_lineage_barplot.pdf/.tiff
#
# Dependencies: ape, ggplot2, dplyr, tidyr
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ape); library(ggplot2); library(dplyr); library(tidyr)
})

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
MSY_TREE <- "path/to/msy_20mb_full.contree"    # rooted MSY consensus tree
MT_TREE  <- "path/to/csa_male_mtdna26.contree" # rooted mtDNA consensus tree
META_TSV <- "data/metadata.tsv"   # columns: sample, group, status
OUT_DIR  <- "results/haplogroup"
FIG_DIR  <- "results/figures"
# ─────────────────────────────────────────────────────────────────────────────

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Tip-label convention: last character encodes genus
#   g = L. g. guanicoe / L. g. cacsilensis
#   l = L. glama
#   v = V. v. vicugna / V. v. mensalis
#   a = V. pacos (alpaca)
LAMA_SUFFIX   <- c("g", "l")
VICUGNA_SUFFIX <- c("v", "a")

GROUP_LABEL <- c(
  L_g_cacsilensis = "L. g. cacsilensis", L_g_guanicoe = "L. g. guanicoe",
  L_glama         = "L. glama",          V_v_mensalis  = "V. v. mensalis",
  V_v_vicugna     = "V. v. vicugna",     V_pacos       = "V. pacos"
)
SP_ORDER <- c("L. g. cacsilensis","L. g. guanicoe","L. glama",
              "V. v. mensalis","V. v. vicugna","V. pacos")
EXPECTED_LINEAGE <- c(
  L_g_cacsilensis = "H-Lama",   L_g_guanicoe = "H-Lama", L_glama = "H-Lama",
  V_v_mensalis    = "H-Vicugna",V_v_vicugna  = "H-Vicugna", V_pacos = "H-Vicugna"
)
HG_COLORS <- c("H-Vicugna" = "#0072B2", "H-Lama" = "#E69F00",
                "Discordant" = "#D55E00", "Unresolved" = "#999999")

# ── Load metadata ─────────────────────────────────────────────────────────────
meta <- read.table(META_TSV, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
names(meta) <- tolower(names(meta))
meta$sample <- trimws(meta$sample); meta$group <- trimws(meta$group)

# ── Haplogroup assignment function ───────────────────────────────────────────
assign_haplogroups <- function(tree_file, marker) {
  tr   <- read.tree(tree_file)
  tips <- tr$tip.label
  suf  <- setNames(substr(tips, nchar(tips), nchar(tips)), tips)

  ingroup_tips <- tips[suf %in% c(LAMA_SUFFIX, VICUGNA_SUFFIX)]
  suf_in       <- suf[ingroup_tips]
  vic_alp      <- ingroup_tips[suf_in %in% VICUGNA_SUFFIX]

  mrca_va   <- getMRCA(tr, vic_alp)
  clade_va  <- extract.clade(tr, mrca_va)$tip.label
  suf_clade <- substr(clade_va, nchar(clade_va), nchar(clade_va))
  pct_va    <- sum(suf_clade %in% VICUGNA_SUFFIX) / length(suf_clade)

  cat(sprintf("\n[%s] Vicugna MRCA clade: %d tips, purity %.1f%%\n",
              marker, length(clade_va), pct_va * 100))
  print(table(suf_clade))

  if (pct_va >= 0.80) {
    hv_tips <- clade_va[suf_clade %in% VICUGNA_SUFFIX]
    hl_tips <- ingroup_tips[!ingroup_tips %in% clade_va]
    introgressed <- clade_va[suf_clade %in% LAMA_SUFFIX]
    method <- "MRCA_strict"
  } else {
    n_tips <- length(tr$tip.label)
    internal_nodes <- (n_tips + 1):(n_tips + tr$Nnode)
    best_node <- NA; best_n <- 0; best_pur <- 0
    for (nd in internal_nodes) {
      cl  <- extract.clade(tr, nd)$tip.label
      sf  <- substr(cl, nchar(cl), nchar(cl))
      pur <- sum(sf %in% VICUGNA_SUFFIX) / length(cl)
      n_va <- sum(sf %in% VICUGNA_SUFFIX)
      if (pur >= 0.80 && n_va > best_n) {
        best_node <- nd; best_n <- n_va; best_pur <- pur
      }
    }
    if (!is.na(best_node)) {
      clade_va  <- extract.clade(tr, best_node)$tip.label
      suf_clade <- substr(clade_va, nchar(clade_va), nchar(clade_va))
      hv_tips   <- clade_va[suf_clade %in% VICUGNA_SUFFIX]
      hl_tips   <- ingroup_tips[!ingroup_tips %in% clade_va]
      introgressed <- clade_va[suf_clade %in% LAMA_SUFFIX]
      method <- "best_pure_subclade"
      cat(sprintf("  Best subclade: %d V+A tips, purity %.1f%%\n", best_n, best_pur*100))
    } else {
      cat("  No pure clade found; using taxonomic fallback.\n")
      hv_tips <- ingroup_tips[suf_in %in% VICUGNA_SUFFIX]
      hl_tips <- ingroup_tips[suf_in %in% LAMA_SUFFIX]
      introgressed <- character(0)
      method <- "taxonomic_fallback"
    }
  }

  hg <- setNames(rep("Unresolved", length(ingroup_tips)), ingroup_tips)
  hg[hv_tips]    <- "H-Vicugna"
  hg[hl_tips]    <- "H-Lama"
  hg[introgressed] <- "Discordant"
  va_outside <- hl_tips[substr(hl_tips, nchar(hl_tips), nchar(hl_tips)) %in% VICUGNA_SUFFIX]
  hg[va_outside] <- "Discordant"

  df_hg <- data.frame(sample = ingroup_tips, haplogroup = hg[ingroup_tips],
                      marker = marker, stringsAsFactors = FALSE)
  df_hg <- merge(df_hg, meta[, c("sample","group")], by = "sample", all.x = TRUE)
  df_hg$sp_label    <- GROUP_LABEL[df_hg$group]
  df_hg$expected_hg <- EXPECTED_LINEAGE[df_hg$group]
  df_hg$concordant  <- df_hg$haplogroup == df_hg$expected_hg

  cat(sprintf("[%s] Assignment:\n", marker))
  print(table(df_hg$group, df_hg$haplogroup))
  df_hg
}

# ── Run for both markers ──────────────────────────────────────────────────────
df_msy <- assign_haplogroups(MSY_TREE, "MSY")
df_mt  <- assign_haplogroups(MT_TREE,  "mtDNA")
df_all <- bind_rows(df_msy, df_mt)

write.csv(df_all, file.path(OUT_DIR, "haplogroup_tree_assignment.csv"), row.names = FALSE)
cat("\nSaved: haplogroup_tree_assignment.csv\n")

# ── Discordance summary ───────────────────────────────────────────────────────
disc <- df_all %>%
  group_by(marker, sp_label) %>%
  summarise(
    n_total       = n(),
    n_H_Lama     = sum(haplogroup == "H-Lama"),
    n_H_Vicugna  = sum(haplogroup == "H-Vicugna"),
    n_Discordant = sum(haplogroup == "Discordant"),
    n_Unresolved = sum(haplogroup == "Unresolved"),
    pct_discordant = round(100 * n_Discordant / n_total, 1),
    .groups = "drop"
  )
cat("\n=== Discordance summary ===\n"); print(as.data.frame(disc))
write.csv(disc, file.path(OUT_DIR, "haplogroup_discordance_table.csv"), row.names = FALSE)

# ── Proportional barplot ──────────────────────────────────────────────────────
df_plot <- df_all %>%
  mutate(
    sp_label   = factor(sp_label,   levels = SP_ORDER),
    haplogroup = factor(haplogroup, levels = c("H-Lama","H-Vicugna","Discordant","Unresolved")),
    marker     = factor(marker,     levels = c("mtDNA","MSY"))
  ) %>%
  group_by(marker, sp_label, haplogroup) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(marker, sp_label) %>%
  mutate(total = sum(n), prop = n / total) %>%
  ungroup()

df_n <- df_all %>%
  mutate(sp_label = factor(sp_label, levels = SP_ORDER),
         marker   = factor(marker,   levels = c("mtDNA","MSY"))) %>%
  group_by(marker, sp_label) %>%
  summarise(total = n(), .groups = "drop")

p <- ggplot(df_plot, aes(x = sp_label, y = prop, fill = haplogroup)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.3) +
  geom_text(data = df_n,
            aes(x = sp_label, y = 1.05, label = paste0("n=", total)),
            inherit.aes = FALSE, size = 3.0, color = "gray20") +
  scale_fill_manual(values = HG_COLORS, name = "Lineage",
                    guide = guide_legend(reverse = FALSE)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.12), expand = expansion(mult = c(0, 0))) +
  facet_wrap(~ marker, ncol = 1, strip.position = "left") +
  labs(x = NULL, y = "Proportion of individuals") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(face = "italic", size = 10, angle = 35,
                                     hjust = 1, color = "gray15"),
    axis.text.y      = element_text(size = 9.5),
    strip.background = element_rect(fill = "#444444", color = NA),
    strip.text       = element_text(color = "white", face = "bold", size = 11,
                                     margin = margin(4, 6, 4, 6)),
    strip.placement  = "outside",
    legend.position  = "top", legend.direction = "horizontal",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 10),
    panel.spacing    = unit(0.8, "lines"),
    plot.margin      = margin(10, 15, 8, 8)
  )

ggsave(file.path(FIG_DIR, "haplogroup_lineage_barplot.pdf"),
       p, width = 6.5, height = 6.0, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "haplogroup_lineage_barplot.tiff"),
       p, width = 6.5, height = 6.0, device = "tiff",
       dpi = 300, compression = "lzw")
cat("Saved: haplogroup_lineage_barplot.pdf / .tiff\n")
