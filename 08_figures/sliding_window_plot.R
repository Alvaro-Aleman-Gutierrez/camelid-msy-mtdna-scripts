#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Sliding-window figure (plot-only from pre-computed CSVs)
#
# Paper: "Sex-biased introgression during South American camelid domestication revealed
#         by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold"
#
# Reads the pre-computed CSVs written by 02_diversity_neutrality/msy_diversity_sliding_window.R
# and regenerates the figure in seconds (no PopGenome re-computation needed).
#
# Inputs:
#   PI_CSV — msy_pi_windows_kb.csv      (columns: pos_kb, pop1, pop2, …)
#   TD_CSV — msy_tajimaD_windows_kb.csv (columns: pos_kb, pop1, pop2, …)
#
# Output:
#   msy20mb_sliding_window_scan.pdf / .tiff
#
# Dependencies: ggplot2, patchwork
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
PI_CSV  <- "results/diversity/msy_pi_windows_kb.csv"
TD_CSV  <- "results/diversity/msy_tajimaD_windows_kb.csv"
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
PAL <- c(
  L_g_cacsilensis = "#000000", L_g_guanicoe = "#E69F00",
  L_glama         = "#009E73", V_v_mensalis  = "#0072B2",
  V_v_vicugna     = "#CC79A7", V_pacos       = "#56B4E9"
)
LTY_MAP <- c(
  L_g_cacsilensis = "solid",  L_g_guanicoe = "dashed",
  L_glama         = "solid",  V_v_mensalis  = "solid",
  V_v_vicugna     = "dashed", V_pacos       = "solid"
)

pi_win <- read.csv(PI_CSV, check.names = FALSE)
td_win <- read.csv(TD_CSV, check.names = FALSE)

pop_names    <- intersect(POP_ORDER, names(pi_win))
midpoints_kb <- pi_win$pos_kb

df_pi <- do.call(rbind, lapply(pop_names, function(p)
  data.frame(pos_kb = midpoints_kb, value = pi_win[[p]], population = p)))
df_D <- do.call(rbind, lapply(pop_names, function(p)
  data.frame(pos_kb = midpoints_kb, value = td_win[[p]], population = p)))
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
    theme(
      axis.text          = element_text(size = 9, color = "gray15"),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.25),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.25),
      legend.position    = "top", legend.direction = "horizontal",
      legend.text        = element_text(size = 9, face = "italic"),
      legend.key.width   = unit(1.4, "cm"),
      legend.background  = element_rect(fill = NA, color = NA),
      plot.margin        = margin(8, 12, 4, 8)
    )
}

p_pi <- make_panel(df_pi, expression(pi[MSY] ~ (per ~ site)),
                   list(list(y = 0, lty = "longdash", col = "gray70")))
p_D  <- make_panel(df_D,  "Tajima's D",
                   list(list(y =  0, lty = "longdash", col = "gray70"),
                        list(y =  2, lty = "dotted",   col = "gray60"),
                        list(y = -2, lty = "dotted",   col = "gray60")))

combined <- (p_pi / p_D) + plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave(file.path(FIG_DIR, "msy20mb_sliding_window_scan.pdf"),
       combined, width = 9, height = 6.5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "msy20mb_sliding_window_scan.tiff"),
       combined, width = 9, height = 6.5, device = "tiff",
       dpi = 300, compression = "lzw")
cat("Saved: msy20mb_sliding_window_scan.pdf / .tiff\n")
