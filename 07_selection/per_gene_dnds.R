#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Per-gene dN/dS (ω) — 22 MSY protein-coding genes
#
# Paper: "Sex-biased demographic history of South American camelids revealed
#         by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold"
#
# Method: All pairwise ω via seqinr::kaks() (Nei-Gojobori 1986)
#   - One _exons.fasta per gene (aligned, trimmed to CDS length)
#   - Internal stop codons masked as "---" (gap) before analysis
#   - Pairs with dS ≤ 0 or ω > 10 excluded (artifacts of near-zero dS)
#   - Summary row: weighted ω = grand_mean_dN / grand_mean_dS
#
# Input : GENE_DIR — directory with *_exons.fasta files
# Output: supp_table_dnds_per_gene.csv
#
# Dependencies: seqinr, ape
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(seqinr); library(ape) })

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
GENE_DIR <- "path/to/gene_cds_fastas/"   # one *_exons.fasta per gene
OUT_DIR  <- "results/selection"
OUTGROUP <- "NC_087473.1_mCamDro1_y-chr"  # removed before analysis
# ─────────────────────────────────────────────────────────────────────────────

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

STOP_CODONS <- c("taa","tag","tga")

# Known MSY gene functions
gene_function <- c(
  AMELY   = "Enamel matrix protein; male sex determination marker",
  CLDN34Y = "Claudin family; testis-expressed tight junction protein",
  CUL4BY  = "Cullin RING ubiquitin ligase; spermatogenesis",
  DDX3Y   = "RNA helicase; required for spermatogenesis (azoospermia factor b)",
  EIF1AY  = "Translation initiation factor; broadly expressed",
  EIF2S3Y = "Translation initiation factor; male-specific copy",
  HSFY1   = "Heat shock transcription factor; expressed in testis",
  HSFY2   = "Heat shock transcription factor; expressed in testis",
  KDM5D   = "Histone H3K4 demethylase; minor histocompatibility antigen HY",
  NLGN4Y  = "Neuroligin 4 Y-linked; synaptic adhesion molecule",
  OFD1Y   = "Oral-facial-digital syndrome 1 Y-linked; ciliogenesis",
  PRKY    = "Protein kinase Y-linked; spermatogenic function unclear",
  SRY     = "Sex-determining region Y; master regulator of male sex determination",
  STSY    = "Suppressor of tumorigenicity Y-linked; MSY-specific",
  TBL1Y   = "Transducin beta-like 1 Y-linked; broadly expressed",
  TMSB4Y  = "Thymosin beta-4 Y-linked; actin dynamics, spermatogenesis",
  TSPY1   = "Testis-specific protein Y-linked; gonocyte proliferation",
  TXNLGY  = "Thioredoxin-like glycoprotein Y-linked; sperm function",
  USP9Y   = "Ubiquitin-specific protease 9 Y; spermatogenesis (azoospermia factor a)",
  UTY     = "KDM6C; histone H3K27 demethylase; minor histocompatibility antigen",
  WWC3Y   = "WW and C2 domain-containing Y-linked; function unclear",
  ZFY     = "Zinc finger protein Y-linked; transcriptional regulator, spermatogenesis"
)

mask_stops <- function(sv) {
  for (i in seq(1, length(sv) - 2, by = 3))
    if (paste(sv[i:(i+2)], collapse = "") %in% STOP_CODONS)
      sv[i:(i+2)] <- c("-","-","-")
  sv
}

gene_files <- sort(list.files(GENE_DIR, pattern = "_exons\\.fasta$", full.names = TRUE))
cat(sprintf("Genes found: %d\n", length(gene_files)))
if (length(gene_files) == 0) stop("No *_exons.fasta files found in GENE_DIR: ", GENE_DIR)

results <- list()

for (fasta in gene_files) {
  gene <- sub("_exons\\.fasta$", "", basename(fasta))
  cat(sprintf("  Processing %s ...", gene))

  seqs_raw <- tryCatch(
    read.fasta(fasta, seqtype = "DNA", forceDNAtolower = TRUE),
    error = function(e) NULL)
  if (is.null(seqs_raw)) { cat(" ERROR reading\n"); next }

  seqs_raw <- seqs_raw[!names(seqs_raw) %in% OUTGROUP]
  n_seq    <- length(seqs_raw)
  cds_len  <- length(seqs_raw[[1]])
  seqs_masked <- lapply(seqs_raw, mask_stops)

  aln <- list(nb = n_seq, nam = names(seqs_masked),
               seq = sapply(seqs_masked, function(s) paste(s, collapse="")),
               com = rep(NA, n_seq))
  class(aln) <- "alignment"

  kk <- tryCatch(kaks(aln, verbose = FALSE), error = function(e) NULL)
  if (is.null(kk)) { cat(" ERROR kaks\n"); next }

  ka_vals <- kk$ka[upper.tri(kk$ka)]
  ks_vals <- kk$ks[upper.tri(kk$ks)]

  valid <- is.finite(ka_vals) & is.finite(ks_vals) & ka_vals >= 0 & ks_vals > 0
  ka_v  <- ka_vals[valid]; ks_v <- ks_vals[valid]; w_v <- ka_v / ks_v
  w_ok  <- is.finite(w_v) & w_v >= 0 & w_v <= 10
  w_v <- w_v[w_ok]; ka_v <- ka_v[w_ok]; ks_v <- ks_v[w_ok]

  n_pairs <- length(w_v)
  if (n_pairs == 0) {
    mean_dN <- mean_dS <- mean_w <- median_w <- NA_real_
  } else {
    mean_dN  <- round(mean(ka_v), 5)
    mean_dS  <- round(mean(ks_v), 5)
    mean_w   <- round(mean(w_v),  4)
    median_w <- round(median(w_v), 4)
  }
  cat(sprintf(" n=%d pairs, dN=%.4f, dS=%.4f, ω=%.4f\n", n_pairs, mean_dN, mean_dS, mean_w))

  results[[gene]] <- data.frame(
    Gene = gene, CDS_length_bp = cds_len, n_codons = cds_len %/% 3,
    n_pairs = n_pairs, mean_dN = mean_dN, mean_dS = mean_dS,
    mean_omega = mean_w, median_omega = median_w, Function = NA_character_,
    stringsAsFactors = FALSE
  )
}

df <- do.call(rbind, lapply(names(results), function(gene) {
  r <- results[[gene]]
  fn <- gene_function[gene]
  r$Function <- if (is.na(fn)) "Function unknown" else fn
  r
}))
df <- df[order(df$Gene), ]
rownames(df) <- NULL

g_dN <- mean(df$mean_dN, na.rm = TRUE)
g_dS <- mean(df$mean_dS, na.rm = TRUE)
overall_row <- data.frame(
  Gene = "ALL GENES (concatenated)", CDS_length_bp = sum(df$CDS_length_bp),
  n_codons = sum(df$n_codons), n_pairs = NA,
  mean_dN = round(g_dN, 5), mean_dS = round(g_dS, 5),
  mean_omega = round(g_dN / g_dS, 4),
  median_omega = round(median(df$median_omega, na.rm = TRUE), 4),
  Function = "Concatenated 22-gene alignment (reported in main text)",
  stringsAsFactors = FALSE
)
df_out <- rbind(df, overall_row)

out_csv <- file.path(OUT_DIR, "supp_table_dnds_per_gene.csv")
write.csv(df_out, out_csv, row.names = FALSE, quote = TRUE)
cat(sprintf("\nSaved: %s\n", out_csv))
cat(sprintf("Genes processed: %d\n", nrow(df)))
w_vals <- df$mean_omega[!is.na(df$mean_omega)]
cat(sprintf("omega range: %.4f – %.4f (%d genes)\n",
            min(w_vals), max(w_vals), length(w_vals)))
print(df[, c("Gene","CDS_length_bp","mean_omega","Function")])
