#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Ne_f / Ne_m ratio plot  (π_mtDNA / π_MSY)
#
# Paper: "Sex-biased introgression during South American camelid domestication revealed
#         by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold"
#
# Under equal per-site mutation rates, π_mtDNA / π_MSY ≈ Ne_f / Ne_m.
# The expected ratio for panmictic species with equal sex ratio is ~4.
#
# Input : RATIO_CSV — columns: group, pi_per_site_mt, pi_per_site_msy, ratio
#         (output of 02_diversity_neutrality/msy_diversity_sliding_window.R)
# Output: ne_ratio_plot.pdf/.tiff, ne_ratio_table.csv
#
# Dependencies: ggplot2
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(ggplot2))

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
RATIO_CSV <- "results/diversity/msy_ratio_pi_mt_over_msy.csv"
OUT_DIR   <- "results/ne_ratio"
FIG_DIR   <- "results/figures"
# ─────────────────────────────────────────────────────────────────────────────

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

PAL <- c(
  L_g_cacsilensis = "#000000", L_g_guanicoe = "#E69F00",
  L_glama         = "#009E73", V_v_mensalis  = "#0072B2",
  V_v_vicugna     = "#CC79A7", V_pacos       = "#56B4E9"
)
POP_ORDER <- c("L_g_cacsilensis","L_g_guanicoe","L_glama",
               "V_v_mensalis","V_v_vicugna","V_pacos")
POP_LABELS <- c(
  L_g_cacsilensis = "L. g. cacsilensis", L_g_guanicoe = "L. g. guanicoe",
  L_glama         = "L. glama",          V_v_mensalis  = "V. v. mensalis",
  V_v_vicugna     = "V. v. vicugna",     V_pacos       = "V. pacos"
)
STATUS <- c(L_g_cacsilensis = "Wild", L_g_guanicoe = "Wild",
            L_glama = "Domestic", V_v_mensalis = "Wild",
            V_v_vicugna = "Wild", V_pacos = "Domestic")

ratio_df <- read.csv(RATIO_CSV, stringsAsFactors = FALSE)

# Accept either column name for the ratio
if ("ratio_pi_mt_msy" %in% names(ratio_df)) {
  names(ratio_df)[names(ratio_df) == "ratio_pi_mt_msy"] <- "Ne_ratio"
} else if (!("Ne_ratio" %in% names(ratio_df))) {
  mt_col  <- grep("pi_per_site_mt",  names(ratio_df), value = TRUE)[1]
  msy_col <- grep("pi_per_site_msy", names(ratio_df), value = TRUE)[1]
  ratio_df$Ne_ratio <- ratio_df[[mt_col]] / ratio_df[[msy_col]]
}

ratio_df$label  <- POP_LABELS[ratio_df$group]
ratio_df$status <- STATUS[ratio_df$group]
ratio_df$label  <- factor(ratio_df$label, levels = POP_LABELS[POP_ORDER])

cat("=== Ne_f / Ne_m ratio ===\n")
print(ratio_df[order(match(ratio_df$group, POP_ORDER)),
               c("group","Ne_ratio")])

p <- ggplot(ratio_df, aes(x = label, y = Ne_ratio, fill = group)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = PAL, guide = "none") +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  labs(
    x     = NULL,
    y     = expression(N[e*","*f] ~ "/" ~ N[e*","*m] ~
               "  (" * pi[mtDNA] / pi[MSY] * ")"),
    title = expression("Effective population size ratio  " *
               N[e*","*f] ~ "/" ~ N[e*","*m])
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, face = "italic", size = 9),
    plot.title   = element_text(size = 11, face = "bold")
  )

ggsave(file.path(FIG_DIR, "ne_ratio_plot.pdf"),
       p, width = 6, height = 4, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "ne_ratio_plot.tiff"),
       p, width = 6, height = 4, device = "tiff",
       dpi = 300, compression = "lzw")
cat("Saved: ne_ratio_plot.pdf / .tiff\n")

out_df <- ratio_df[order(match(ratio_df$group, POP_ORDER)),
                   c("group","label","status","Ne_ratio")]
write.csv(out_df, file.path(OUT_DIR, "ne_ratio_table.csv"), row.names = FALSE)
cat("Saved: ne_ratio_table.csv\n")
