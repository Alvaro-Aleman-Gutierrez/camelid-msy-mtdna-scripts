#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Fisher's exact test for mtDNA haplogroup discordance in V. pacos
#
# Tests whether the proportion of individuals with H-Lama haplogroup differs
# significantly between V. pacos (alpaca) and wild Vicugna subspecies.
#
# Input:
#   HAPLOGROUP_CSV — output of 06_haplogroup/haplogroup_assignment.R
#                    (haplogroup_discordance_table.csv)
# Output:
#   fisher_discordance.csv — OR, 95% CI, p-value
#
# Dependencies: base R only
# ─────────────────────────────────────────────────────────────────────────────

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
HAPLOGROUP_CSV <- "results/haplogroup/haplogroup_discordance_table.csv"
OUT_DIR        <- "results/haplogroup"
# ─────────────────────────────────────────────────────────────────────────────

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

hg <- read.csv(HAPLOGROUP_CSV, stringsAsFactors = FALSE)

# Extract mtDNA rows
hg_mt <- hg[hg$marker == "mtDNA", ]

vpacos  <- hg_mt[hg_mt$sp_label == "V. pacos", ]
vvic    <- hg_mt[hg_mt$sp_label %in% c("V. v. vicugna", "V. v. mensalis"), ]

# V. pacos: individuals in H-Lama are "Discordant"; H-Vicugna are concordant
n_lama_pacos   <- vpacos$n_Discordant
n_vic_pacos    <- vpacos$n_H_Vicugna
n_lama_vicugna <- sum(vvic$n_H_Lama)
n_vic_vicugna  <- sum(vvic$n_H_Vicugna)

mat <- matrix(
  c(n_lama_pacos, n_vic_pacos, n_lama_vicugna, n_vic_vicugna),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("V. pacos", "V. vicugna (pooled)"),
                  c("H_Lama", "H_Vicugna"))
)
cat("Contingency table (mtDNA):\n"); print(mat)

ft <- fisher.test(mat)
cat(sprintf(
  "\nFisher exact test:\n  OR    = %.3f\n  95%% CI = [%.3f, %.3f]\n  p     = %.2e\n",
  ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value
))

res <- data.frame(
  test      = "Fisher_mtDNA_Vpacos_vs_Vicugna",
  OR        = round(ft$estimate, 4),
  CI_lo     = round(ft$conf.int[1], 4),
  CI_hi     = round(ft$conf.int[2], 4),
  p_value   = ft$p.value
)
write.csv(res, file.path(OUT_DIR, "fisher_discordance.csv"), row.names = FALSE)
cat("Saved: fisher_discordance.csv\n")
