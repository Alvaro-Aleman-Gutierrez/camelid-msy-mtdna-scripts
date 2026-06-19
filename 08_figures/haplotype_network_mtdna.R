#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# mtDNA haplotype network (TCS-like parsimony via pegas)
#
# Paper: "Sex-biased demographic history of South American camelids revealed
#         by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold"
#
# Uses only variable (segregating) sites to avoid crowded branch labels.
# Node size is proportional to sqrt(haplotype frequency).
# Pie-chart composition per node shows population membership.
#
# Inputs:
#   MT_FASTA — mtDNA alignment without outgroup
#   META_TSV — metadata: columns sample, group
#
# Output:
#   haplotype_network_mtdna.tiff (3200×2400 px, 300 dpi)
#
# Dependencies: ape, pegas
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(ape); library(pegas) })

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
MT_FASTA <- "path/to/mtdna_wo_outgroup_aln.fasta"
META_TSV <- "data/metadata.tsv"   # columns: sample, group
FIG_DIR  <- "results/figures"
# ─────────────────────────────────────────────────────────────────────────────

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

PAL <- c(
  L_g_cacsilensis = "#000000", L_g_guanicoe = "#E69F00",
  L_glama         = "#009E73", V_pacos       = "#56B4E9",
  V_v_mensalis    = "#0072B2", V_v_vicugna   = "#CC79A7"
)

# ── Load data ─────────────────────────────────────────────────────────────────
mt <- read.dna(MT_FASTA, format = "fasta")
meta <- read.table(META_TSV, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
names(meta) <- tolower(names(meta))
meta$sample <- trimws(meta$sample); meta$group <- trimws(meta$group)

ids  <- intersect(rownames(mt), meta$sample)
mt   <- mt[ids, ]
meta <- meta[match(ids, meta$sample), ]
stopifnot(all(meta$sample == ids))

# ── Reduce to variable sites ──────────────────────────────────────────────────
var_pos <- seg.sites(mt)
cat(sprintf("Total sites: %d  |  Variable sites: %d\n", ncol(mt), length(var_pos)))
mt_var <- mt[, var_pos]

# ── Build haplotype network ───────────────────────────────────────────────────
h   <- haplotype(mt_var)
net <- haploNet(h)
cat(sprintf("Unique haplotypes: %d\n", nrow(h)))

ind.hap <- with(stack(setNames(attr(h, "index"), rownames(h))),
                table(hap = ind, pop = meta$group[values]))
freq    <- attr(net, "freq")
node_sz <- pmax(sqrt(freq) * 2.5, 1)

# ── Plot ──────────────────────────────────────────────────────────────────────
out_tiff <- file.path(FIG_DIR, "haplotype_network_mtdna.tiff")
tiff(out_tiff, width = 3200, height = 2400, res = 300, compression = "lzw")
par(mar = c(2, 2, 3, 2))

plot(net,
     size          = node_sz,
     pie           = ind.hap,
     bg            = PAL,
     legend        = FALSE,
     show.mutation = 2,
     labels        = FALSE,
     cex           = 0.8)

legend("bottomright",
       legend    = c("L. g. cacsilensis","L. g. guanicoe","L. glama",
                     "V. pacos","V. v. mensalis","V. v. vicugna"),
       fill      = PAL,
       border    = NA, bty = "n", cex = 0.8, text.font = 3)

title("Haplotype network — mtDNA (variable sites)", cex.main = 1.1)
dev.off()
cat("Saved:", out_tiff, "\n")
