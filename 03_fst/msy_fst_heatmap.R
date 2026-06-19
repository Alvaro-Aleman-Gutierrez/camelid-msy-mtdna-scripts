#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Pairwise FST heat map — MSY scaffold
#
# Input : msy_fst_pairwise.csv  (output of 02_diversity_neutrality/msy_diversity_sliding_window.R)
# Output: msy_fst_heatmap.pdf / .tiff
#
# Dependencies: ggplot2, tidyr, dplyr, tibble
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2); library(tidyr); library(dplyr); library(tibble)
})

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
FST_CSV <- "results/diversity/msy_fst_pairwise.csv"
FIG_DIR <- "results/figures"
# ─────────────────────────────────────────────────────────────────────────────

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

LAB_ORDER <- c("L. g. cacsilensis", "L. g. guanicoe", "L. glama",
               "V. v. mensalis",    "V. v. vicugna",  "V. pacos")

fst_raw <- read.csv(FST_CSV, row.names = 1, check.names = FALSE)
fst_mat <- as.matrix(fst_raw)
fst_sym <- (fst_mat + t(fst_mat)) / 2
diag(fst_sym) <- NA
fst_sym <- fst_sym[LAB_ORDER, LAB_ORDER]

cat("=== Pairwise FST — MSY scaffold ===\n"); print(round(fst_sym, 4))

df <- as.data.frame(fst_sym) |>
  rownames_to_column("Pop1") |>
  pivot_longer(-Pop1, names_to = "Pop2", values_to = "FST")

pos1 <- match(df$Pop1, LAB_ORDER)
pos2 <- match(df$Pop2, LAB_ORDER)
df$FST[pos1 > pos2] <- NA
df$lbl <- ifelse(is.na(df$FST), "", sprintf("%.3f", df$FST))

df$Pop1 <- factor(df$Pop1, levels = rev(LAB_ORDER))
df$Pop2 <- factor(df$Pop2, levels = LAB_ORDER)
max_fst <- max(df$FST, na.rm = TRUE)

p <- ggplot(df, aes(x = Pop2, y = Pop1, fill = FST)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = lbl), size = 3.4, color = "white", fontface = "bold",
            na.rm = TRUE) +
  geom_text(data = subset(df, !is.na(FST) & FST < max_fst * 0.55),
            aes(label = lbl), size = 3.4, color = "gray20", fontface = "bold",
            na.rm = TRUE) +
  scale_fill_gradientn(
    colours  = c("#F7FBFF","#DEEBF7","#9ECAE1","#4292C6","#08519C","#08306B"),
    na.value = "#F2F2F2",
    name     = expression(F[ST]),
    limits   = c(0, max_fst),
    guide    = guide_colorbar(barwidth = 0.7, barheight = 5.5,
                               title.position = "top", title.hjust = 0.5)
  ) +
  coord_equal(clip = "off") +
  labs(x = NULL, y = NULL,
       title = expression(italic(F)[ST]~"pairwise — MSY scaffold (~20.6 Mb)")) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(face = "italic", size = 9.5, angle = 35, hjust = 1),
    axis.text.y  = element_text(face = "italic", size = 9.5, hjust = 1),
    panel.grid   = element_blank(),
    legend.position = "right",
    plot.title   = element_text(size = 11, hjust = 0.5, margin = margin(b = 8))
  ) +
  # Dashed line separating Lama and Vicugna genera
  annotate("segment", x = 3.5, xend = 3.5, y = 0.4, yend = 6.6,
           color = "gray50", linewidth = 0.5, linetype = "dashed") +
  annotate("segment", x = 0.4, xend = 6.6, y = 3.5, yend = 3.5,
           color = "gray50", linewidth = 0.5, linetype = "dashed")

ggsave(file.path(FIG_DIR, "msy_fst_heatmap.pdf"),
       p, width = 6.2, height = 5.8, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "msy_fst_heatmap.tiff"),
       p, width = 6.2, height = 5.8, device = "tiff",
       dpi = 300, compression = "lzw")
cat("Saved: msy_fst_heatmap.pdf / .tiff\n")
