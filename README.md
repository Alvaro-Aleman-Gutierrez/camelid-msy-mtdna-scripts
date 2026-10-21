# Analysis scripts — *Sex-biased introgression during South American camelid domestication revealed by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold*

**Authors:** Álvaro Alemán-Gutiérrez, Andrés Peña-Monroy, *et al.*  
**Journal:** (in review)  
**Preprint / DOI:** *pending*  
**Repository:** https://github.com/Alvaro-Aleman-Gutierrez/camelid-msy-mtdna-scripts

---

## Overview

This repository contains all custom R and shell scripts used to reproduce the phylogenomic, population-genetic, and molecular-evolution analyses presented in the paper. Analyses compare a complete mitochondrial genome alignment (16,643 bp, 45 sequences) with a chromosome-level Y-chromosome scaffold alignment (20,577,031 bp, 45 sequences) across six South American camelid populations.

| Population | Code | Status |
|---|---|---|
| *Lama guanicoe cacsilensis* | `L_g_cacsilensis` | Wild |
| *Lama guanicoe guanicoe* | `L_g_guanicoe` | Wild |
| *Lama glama* | `L_glama` | Domestic |
| *Vicugna vicugna mensalis* | `V_v_mensalis` | Wild |
| *Vicugna vicugna vicugna* | `V_v_vicugna` | Wild |
| *Vicugna pacos* | `V_pacos` | Domestic |

Outgroup: *Camelus dromedarius* (mtDNA: NC_009849.1; MSY: NC_087473.1).

---

## Repository structure

```
├── 00_read_processing/
│   ├── 01_mapping.sh               # BWA-MEM v0.7.17 — read mapping (MSY or mtDNA)
│   ├── 02_add_readgroups.sh        # SAMtools — add @RG tags required by GATK
│   ├── 03_variant_calling.sh       # GATK v4.1.5.0 HaplotypeCaller (--ploidy 1)
│   ├── 04_filter_consensus.sh      # BCFtools — SNP filtering + consensus FASTA
│   ├── 05_merge_vcf.sh             # BCFtools merge — population-level VCF
│   └── 06_align_mafft.sh           # MAFFT v7.427 — multiple sequence alignment
│
├── 01_phylogenetics/
│   ├── iqtree_msy_20mb.sh          # IQ-TREE 3 — MSY scaffold (~27 h, 30 threads)
│   └── iqtree_mtdna.sh             # IQ-TREE 3 — complete mtDNA
│
├── 02_diversity_neutrality/
│   └── msy_diversity_sliding_window.R   # π, θW, Tajima's D, FST; 50-kb sliding windows
│
├── 03_fst/
│   └── msy_fst_heatmap.R          # Pairwise FST heatmap (Fig. 3)
│
├── 04_pcoa_amova/
│   └── pcoa_amova.R               # PCoA (K80 distances) + AMOVA (Fig. 4)
│
├── 05_mantel/
│   └── mantel_scatter.R           # Mantel tests + scatter plots (Fig. 5)
│
├── 06_haplogroup/
│   ├── haplogroup_assignment.R    # Tree-topology haplogroup assignment (Fig. 6)
│   └── fisher_discordance.R       # Fisher's exact test for V. pacos discordance
│
├── 07_selection/
│   ├── per_gene_dnds.R            # Per-gene dN/dS via seqinr::kaks() (Table S2)
│   ├── run_busted_all.sh          # BUSTED episodic selection on 22 MSY genes (HyPhy v2.5.2)
│   ├── parse_busted_results.py    # Parse BUSTED JSON output → Table S3 CSV
│   └── busted_lollipop.R          # Lollipop plot of BUSTED p-values (Supp. Fig. S5)
│
├── 08_figures/
│   ├── ne_ratio_plot.R            # Ne_f / Ne_m ratio bar chart (Fig. 7)
│   ├── tajima_panel.R             # Tajima's D dumbbell chart (Fig. 2)
│   ├── sliding_window_plot.R      # Sliding-window figure from pre-computed CSVs
│   └── haplotype_network_mtdna.R  # mtDNA haplotype network via pegas (Fig. S3)
│
└── data/
    └── metadata.tsv               # Template: sample → population mapping
```

---

## Dependencies

### Shell scripts
| Tool | Version | Purpose |
|---|---|---|
| BWA | 0.7.17 | Read mapping (`00_read_processing/01`) |
| SAMtools | 1.10 | BAM processing and read groups (`00_read_processing/01–02`) |
| GATK | 4.1.5.0 | Variant calling (`00_read_processing/03`) |
| BCFtools | 1.10.2 | SNP filtering, consensus FASTA, VCF merge (`00_read_processing/04–05`) |
| MAFFT | 7.427 | Multiple sequence alignment (`00_read_processing/06`) |
| GNU Parallel | ≥ 20200101 | Parallelisation across samples |
| IQ-TREE 3 | 3.1.2 | Phylogenetic inference (`01_phylogenetics/`) |

Install via conda:
```bash
conda install -c bioconda bwa samtools gatk4 bcftools mafft parallel iqtree
```

### R scripts
| Package | Version | Purpose |
|---|---|---|
| `PopGenome` | 2.7.5 | π, θW, Tajima's D, FST, sliding windows |
| `ggrepel` | 0.9.x | Non-overlapping labels (optional; not required for BUSTED lollipop) |
| `ape` | 5.8.1 | Sequence I/O, distance matrices, PCoA, tree parsing |
| `pegas` | 1.4 | AMOVA, haplotype networks |
| `vegan` | 2.7.5 | Mantel tests |
| `seqinr` | 4.2.44 | Per-gene dN/dS via `kaks()` |
| `ggtree` | 4.0.4 | Phylogenetic tree visualisation |
| `ggplot2` | 4.0.3 | All figures |
| `patchwork` | 1.3.2 | Multi-panel figure assembly |
| `dplyr` | 1.2.1 | Data wrangling |
| `tidyr` | 1.3.2 | Data wrangling |
| `tibble` | 3.3.1 | Data wrangling |

Install in R:
```r
install.packages(c("ape", "pegas", "vegan", "seqinr",
                   "ggplot2", "patchwork", "dplyr", "tidyr", "tibble"))

# PopGenome (from CRAN or GitHub)
install.packages("PopGenome")
```

> **Note on memory:** Loading the 20.6 Mb MSY alignment into R (via `ape::read.dna` or `PopGenome::readData`) requires **4–8 GB RAM**. The `readData` step in `02_diversity_neutrality/msy_diversity_sliding_window.R` may take 20–60 min; the subsequent sliding-window computation adds 10–20 min.

---

## Quick start

### 1. Prepare your data

Edit `data/metadata.tsv` to map each sample ID (matching the FASTA header exactly) to one of the six population codes and a `wild`/`domestic` status label. See the file for the expected format.

### 2. Run phylogenetic inference

```bash
# Edit USER CONFIGURATION block in each script first
bash 01_phylogenetics/iqtree_msy_20mb.sh   # ~27 h on 30 threads
bash 01_phylogenetics/iqtree_mtdna.sh
```

The best-fit substitution models (selected by ModelFinder under BIC) are:

| Marker | Model |
|---|---|
| MSY scaffold | TVM+F+I+R5 |
| Complete mtDNA | TPM2u+F+R3 |

Both runs use UFBoot (1000 replicates), SH-aLRT (1000 replicates), 10 independent tree searches (`--runs 10`), and nearest-neighbour interchange optimisation (`-bnni`).

### 3. Population genetics (MSY)

```r
Rscript 02_diversity_neutrality/msy_diversity_sliding_window.R
```

Outputs (in `results/diversity/`):

| File | Content |
|---|---|
| `msy_popgenome_summary.csv` | π, θW, Tajima's D per population |
| `msy_fst_pairwise.csv` | 6×6 nucleotide FST matrix |
| `msy_ratio_pi_mt_over_msy.csv` | π_mtDNA / π_MSY ratio |
| `msy_tajimaD_windows_kb.csv` | Tajima's D per 50-kb window |
| `msy_pi_windows_kb.csv` | π per 50-kb window |

### 4. Downstream analyses and figures

Run scripts in the order listed under [Repository structure](#repository-structure). Each script has a `USER CONFIGURATION` block at the top where you set input/output paths.

```r
Rscript 03_fst/msy_fst_heatmap.R
Rscript 04_pcoa_amova/pcoa_amova.R
Rscript 05_mantel/mantel_scatter.R
Rscript 06_haplogroup/haplogroup_assignment.R
Rscript 06_haplogroup/fisher_discordance.R   # requires haplogroup_assignment.R output
Rscript 07_selection/per_gene_dnds.R
Rscript 08_figures/ne_ratio_plot.R
Rscript 08_figures/tajima_panel.R
Rscript 08_figures/sliding_window_plot.R     # fast: reads pre-computed CSVs
Rscript 08_figures/haplotype_network_mtdna.R
```

---

## Colour palette

All figures use the [Wong (2011)](https://www.nature.com/articles/nmeth.1618) seven-colour colorblind-safe palette:

| Population | Hex | Swatch |
|---|---|---|
| *L. g. cacsilensis* | `#000000` | Black |
| *L. g. guanicoe* | `#E69F00` | Orange |
| *L. glama* | `#009E73` | Bluish green |
| *V. v. mensalis* | `#0072B2` | Blue |
| *V. v. vicugna* | `#CC79A7` | Reddish purple |
| *V. pacos* | `#56B4E9` | Sky blue |

---

## Key results reproduced by these scripts

- **Haplogroups:** Two reciprocally monophyletic lineages on both markers — H-Lama (*Lama* spp.) and H-Vicugna (*Vicugna* spp.) — with 75% of *V. pacos* (alpaca) individuals carrying the H-Lama mtDNA haplogroup (Fisher's exact test, p < 0.001).
- **Sex-biased introgression:** π_mtDNA / π_MSY = 3.82 in *V. pacos* vs. 1.03–1.76 in wild species, consistent with female-biased introgression during domestication.
- **Purifying selection:** Per-gene ω = 0.20–0.24 across 22 MSY coding sequences. Formal BUSTED analysis identified episodic diversifying selection in 4 genes (CUL4BY, HSFY2, ZFY, EIF1AY; p < 0.05), all with roles in spermatogenesis or sex determination. HSFY1 (ω = 1.283) showed no significant BUSTED signal (p = 0.153).
- **Population structure:** MSY FST is uniformly higher than mtDNA FST across all pairwise comparisons, indicating stronger male-biased drift.

---

## Data availability

Raw sequencing reads and genome assemblies: NCBI BioProject **[XXXXXXX]** (to be released upon publication).

Alignments and intermediate outputs can be regenerated from the raw reads following the Methods section of the paper.

---

## Citation

If you use these scripts, please cite:

> Alemán-Gutiérrez A, *et al.* (2026) Sex-biased introgression during South American camelid domestication revealed by complete mitochondrial genomes and a 20.6 Mb Y-chromosome scaffold. *Heredity*. doi: *pending*.

---

## License

Scripts are released under the [MIT License](LICENSE).
