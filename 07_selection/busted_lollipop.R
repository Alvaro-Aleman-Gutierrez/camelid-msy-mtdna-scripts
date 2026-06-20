#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# busted_lollipop.R — Lollipop plot of BUSTED results (22 MSY genes)
#
# Paper: "Sex-biased introgression during South American camelid domestication revealed by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold"
#
# Input : busted_summary_tableS3.csv  (output of parse_busted_results.py)
# Output: busted_lollipop.pdf / .tiff  (Supplementary Figure)
#
# Dependencies: ggplot2, ggrepel
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
INPUT_CSV <- "path/to/busted_summary_tableS3.csv"
FIG_DIR   <- "results/figures"
# ─────────────────────────────────────────────────────────────────────────────

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Wong (2011) colours
COL_SIG   <- "#D55E00"   # vermillion — significant genes
COL_NS    <- "gray65"    # non-significant genes
COL_LINE  <- "#0072B2"   # blue — threshold line

df <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)
df <- df[df$p_value != "n.d.", ]
df$p_value <- as.numeric(df$p_value)
df$neglog10p <- -log10(df$p_value)
df$significant <- df$p_value < 0.05
df$colour <- ifelse(df$significant, COL_SIG, COL_NS)

# Order genes by -log10(p) descending
df$Gene <- factor(df$Gene, levels = df$Gene[order(df$neglog10p)])

threshold <- -log10(0.05)   # 1.301

p <- ggplot(df, aes(x = Gene, y = neglog10p)) +
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = neglog10p),
               colour = df$colour[order(df$neglog10p)],
               linewidth = 0.7) +
  geom_point(aes(colour = significant), size = 3.8) +
  geom_hline(yintercept = threshold, linetype = "dashed",
             colour = COL_LINE, linewidth = 0.55) +
  geom_label_repel(
    data = subset(df, significant),
    aes(label = Gene),
    colour      = COL_SIG,
    fill        = "white",
    size        = 3.2,
    fontface    = "italic",
    label.size  = 0.25,
    box.padding = 0.4,
    nudge_y     = 0.25,
    direction   = "y",
    segment.colour = COL_SIG,
    segment.size   = 0.35
  ) +
  annotate("text", x = 1, y = threshold + 0.12,
           label = "p = 0.05", hjust = 0, vjust = 0,
           size = 3.0, colour = COL_LINE) +
  scale_colour_manual(values = c("TRUE" = COL_SIG, "FALSE" = COL_NS),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                      name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     breaks = 0:5) +
  labs(
    x = NULL,
    y = expression(-log[10](italic(p))),
    title = "Episodic diversifying selection — BUSTED (22 MSY genes)"
  ) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y        = element_text(face = "italic", size = 9.5, colour = "gray15"),
    axis.text.x        = element_text(size = 10, colour = "gray15"),
    axis.title.x       = element_text(size = 11, margin = margin(t = 8)),
    panel.grid.major.x = element_line(colour = "gray90", linewidth = 0.35),
    legend.position    = "top",
    legend.direction   = "horizontal",
    legend.background  = element_rect(fill = NA, colour = NA),
    legend.text        = element_text(size = 10.5),
    plot.title         = element_text(size = 11, hjust = 0.5,
                                      colour = "gray15", margin = margin(b = 10)),
    plot.margin        = margin(12, 18, 10, 10)
  )

ggsave(file.path(FIG_DIR, "busted_lollipop.pdf"),
       p, width = 6.5, height = 7.0, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "busted_lollipop.tiff"),
       p, width = 6.5, height = 7.0, device = "tiff",
       dpi = 300, compression = "lzw")
cat("Saved: busted_lollipop.pdf / .tiff\n")
