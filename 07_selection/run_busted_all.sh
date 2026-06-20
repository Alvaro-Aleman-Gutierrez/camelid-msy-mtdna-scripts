#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# run_busted_all.sh — BUSTED (HyPhy) en los 22 genes CDS del MSY
#
# Uso: bash run_busted_all.sh
# Tiempo estimado: 30–90 min en total (genes corren secuencialmente)
# Resultados: busted_results/json/<GEN>.BUSTED.json
# ─────────────────────────────────────────────────────────────────────────────

# HyPhy necesita libcrypto 1.1
export LD_LIBRARY_PATH=/home/aleman/anaconda3/pkgs/openssl-1.1.1w-hd590300_0/lib:$LD_LIBRARY_PATH

FASTA_DIR="/home/aleman/papers/paper_MSc/final_MSc/busted_results/fastas_clean"
TREE="/home/aleman/papers/paper_MSc/final_MSc/busted_results/msy_pruned.nwk"
JSON_DIR="/home/aleman/papers/paper_MSc/final_MSc/busted_results/json"
LOG_DIR="/home/aleman/papers/paper_MSc/final_MSc/busted_results/logs"

mkdir -p "$JSON_DIR" "$LOG_DIR"

GENES=(AMELY CLDN34Y CUL4BY DDX3Y EIF1AY EIF2S3Y HSFY1 HSFY2 KDM5D \
       NLGN4Y OFD1Y PRKY SRY STSY TBL1Y TMSB4Y TSPY1 TXNLGY \
       USP9Y UTY WWC3Y ZFY)

echo "Iniciando BUSTED en ${#GENES[@]} genes — $(date)"
echo "=================================================="

for GENE in "${GENES[@]}"; do
    FASTA="$FASTA_DIR/${GENE}_clean.fasta"
    JSON_OUT="$JSON_DIR/${GENE}.BUSTED.json"
    LOG="$LOG_DIR/${GENE}.log"

    if [[ -f "$JSON_OUT" ]]; then
        echo "  [SKIP] $GENE — resultado ya existe"
        continue
    fi

    if [[ ! -f "$FASTA" ]]; then
        echo "  [WARN] $GENE — fasta no encontrado, saltando"
        continue
    fi

    echo -n "  Procesando $GENE ... "
    START=$(date +%s)

    hyphy BUSTED \
        --alignment "$FASTA" \
        --tree "$TREE" \
        --branches All \
        --srv Yes \
        --output "$JSON_OUT" \
        > "$LOG" 2>&1

    END=$(date +%s)
    ELAPSED=$((END - START))

    if [[ -f "$JSON_OUT" ]]; then
        echo "OK (${ELAPSED}s)"
    else
        echo "ERROR — revisa $LOG"
    fi
done

echo ""
echo "=================================================="
echo "Completado — $(date)"
echo "Resultados en: $JSON_DIR"
echo ""
echo "Genera la tabla con:"
echo "  python3 parse_busted_results.py"
