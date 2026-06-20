#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Tajima's D dumbbell chart — mtDNA vs MSY per population
#
# Paper: "Sex-biased introgression during South American camelid domestication revealed
#         by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold"
#
# Inputs:
#   MT_CSV  — PopGenome summary CSV for mtDNA (columns: group, TajimaD)
#   MSY_CSV — PopGenome summary CSV for MSY   (columns: group, TajimaD)
#
# Output:
#   tajimaD_mtDNA_vs_MSY.pdf / .tiff
#
# Dependencies: ggplot2
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(ggplot2))

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
MT_CSV  <- "results/diversity/mtdna_popgenome_summary.csv"  # mtDNA summary
MSY_CSV <- "results/diversity/msy_popgenome_summary.csv"    # MSY summary
FIG_DIR <- "results/figures"
# ─────────────────────────────────────────────────────────────────────────────

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

POP_ORDER <- c("L_g_cacsilensis","L_g_guanicoe","L_glama",
               "V_v_mensalis","V_v_vicugna","V_pacos")
POP_LABELS <- c(
  L_g_cacsilensis = "L. g. cacsilensis", L_g_guanicoe = "L. g. guanicoe",
  L_glama         = "L. glama",          V_v_mensalis  = "V. v. mensalis",
  V_v_vicugna     = "V. v. vicugna",     V_pacos       = "V. pacos"
)

# Wong (2011): Vermillion for mtDNA (maternal), Blue for MSY (paternal)
COL_MT  <- "#D55E00"
COL_MSY <- "#0072B2"

mt  <- read.csv(MT_CSV,  stringsAsFactors = FALSE)
msy <- read.csv(MSY_CSV, stringsAsFactors = FALSE)

mt_d  <- mt[,  c("group","TajimaD")]; names(mt_d)[2]  <- "D_mt"
msy_d <- msy[, c("group","TajimaD")]; names(msy_d)[2] <- "D_msy"

df <- merge(mt_d, msy_d, by = "group")
df$group <- factor(df$group, levels = rev(POP_ORDER))
df$label <- factor(POP_LABELS[as.character(df$group)],
                    levels = POP_LABELS[rev(POP_ORDER)])

df_long <- rbind(
  data.frame(label = df$label, D = df$D_mt,  marker = "mtDNA"),
  data.frame(label = df$label, D = df$D_msy, marker = "MSY")
)
df_long$marker <- factor(df_long$marker, levels = c("mtDNA","MSY"))

x_min <- floor(min(df_long$D) * 10) / 10 - 0.15
x_max <- ceiling(max(df_long$D) * 10) / 10 + 0.18

p <- ggplot() +
  annotate("rect", xmin = -2, xmax = 2, ymin = -Inf, ymax = Inf,
           fill = "#F5F5F5", alpha = 1) +
  geom_vline(xintercept = 0, linetype = "longdash",
             color = "gray75", linewidth = 0.45) +
  geom_segment(data = df,
               aes(x = D_mt, xend = D_msy, y = label, yend = label),
               color = "gray70", linewidth = 0.7, alpha = 0.8) +
  geom_point(data = subset(df_long, marker == "mtDNA"),
             aes(x = D, y = label, color = marker),
             size = 4.5, alpha = 0.92) +
  geom_point(data = subset(df_long, marker == "MSY"),
             aes(x = D, y = label, color = marker),
             size = 4.5, shape = 17, alpha = 0.92) +
  geom_text(data = subset(df_long, marker == "mtDNA"),
            aes(x = D, y = label, label = sprintf("%.3f", D)),
            vjust = -0.85, size = 2.9, color = COL_MT) +
  geom_text(data = subset(df_long, marker == "MSY"),
            aes(x = D, y = label, label = sprintf("%.3f", D)),
            vjust = 1.85, size = 2.9, color = COL_MSY) +
  scale_color_manual(values = c("mtDNA" = COL_MT, "MSY" = COL_MSY),
                     name = NULL,
                     guide = guide_legend(override.aes = list(shape = c(16, 17), size = 4))) +
  scale_x_continuous(limits = c(x_min, x_max), breaks = seq(-0.5, 1.5, 0.5)) +
  labs(x = "Tajima's D", y = NULL,
       title = "Tajima's D per population: mtDNA vs. MSY") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y          = element_text(face = "italic", size = 10.5, color = "gray15"),
    axis.text.x          = element_text(size = 10, color = "gray15"),
    axis.title.x         = element_text(size = 12, margin = margin(t = 8)),
    panel.grid.major.x   = element_line(color = "gray88", linewidth = 0.35),
    legend.position      = "top", legend.direction = "horizontal",
    legend.background    = element_rect(fill = NA, color = NA),
    legend.margin        = margin(0, 12, 4, 12),
    legend.text          = element_text(size = 10.5),
    plot.title           = element_text(size = 11.5, hjust = 0.5,
                                        color = "gray15", margin = margin(b = 10)),
    plot.margin          = margin(12, 18, 10, 10)
  )

ggsave(file.path(FIG_DIR, "tajimaD_mtDNA_vs_MSY.pdf"),
       p, width = 6.5, height = 4.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "tajimaD_mtDNA_vs_MSY.tiff"),
       p, width = 6.5, height = 4.5, device = "tiff",
       dpi = 300, compression = "lzw")
cat("Saved: tajimaD_mtDNA_vs_MSY.pdf / .tiff\n")
