#!/usr/bin/env python3
"""
parse_busted_results.py — Consolida resultados BUSTED en Tabla S3

Uso: python3 parse_busted_results.py
Salida: busted_summary_tableS3.csv
"""
import json, os, re

JSON_DIR = "/home/aleman/papers/paper_MSc/final_MSc/busted_results/json"
OUT_CSV  = "/home/aleman/papers/paper_MSc/final_MSc/busted_results/busted_summary_tableS3.csv"

GENES = ["AMELY","CLDN34Y","CUL4BY","DDX3Y","EIF1AY","EIF2S3Y","HSFY1","HSFY2",
         "KDM5D","NLGN4Y","OFD1Y","PRKY","SRY","STSY","TBL1Y","TMSB4Y","TSPY1",
         "TXNLGY","USP9Y","UTY","WWC3Y","ZFY"]

rows = []
for gene in GENES:
    json_path = os.path.join(JSON_DIR, f"{gene}.BUSTED.json")
    if not os.path.exists(json_path):
        rows.append({
            "Gene": gene, "n_sequences": "n.d.", "n_codons": "n.d.",
            "LRT": "n.d.", "p_value": "n.d.",
            "omega1": "n.d.", "prop1": "n.d.",
            "omega2": "n.d.", "prop2": "n.d.",
            "omega3": "n.d.", "prop3": "n.d.",
            "significant": "n.d."
        })
        print(f"  MISSING: {gene}")
        continue

    with open(json_path) as f:
        content = f.read().strip()
    if not content:
        rows.append({
            "Gene": gene, "n_sequences": "n.d.", "n_codons": "n.d.",
            "LRT": "n.d.", "p_value": "n.d.",
            "omega1": "n.d.", "prop1": "n.d.",
            "omega2": "n.d.", "prop2": "n.d.",
            "omega3": "n.d.", "prop3": "n.d.",
            "significant": "running"
        })
        print(f"  {gene}: JSON vacío (análisis en curso)")
        continue
    r = json.loads(content)

    test   = r.get("test results", {})
    pval   = test.get("p-value", None)
    lrt    = test.get("LRT", None)

    # Sequences and codons
    input_info = r.get("input", {})
    n_seq   = input_info.get("number of sequences", "?")
    n_sites = input_info.get("number of sites", "?")

    # Rate classes from unconstrained model
    fits = r.get("fits", {})
    unc  = fits.get("Unconstrained model", {})
    rd   = unc.get("Rate Distributions", {})
    test_rd = rd.get("Test", {})

    omegas = []
    for k in sorted(test_rd.keys()):
        v = test_rd[k]
        if isinstance(v, dict):
            omegas.append((float(v.get("omega", 0)), float(v.get("weight", 0))))

    while len(omegas) < 3:
        omegas.append((None, None))

    sig = "Yes" if (pval is not None and pval < 0.05) else "No"

    rows.append({
        "Gene": gene,
        "n_sequences": n_seq,
        "n_codons": n_sites,
        "LRT": f"{lrt:.4f}" if lrt is not None else "n.d.",
        "p_value": f"{pval:.4f}" if pval is not None else "n.d.",
        "omega1": f"{omegas[0][0]:.4f}" if omegas[0][0] is not None else "n.d.",
        "prop1":  f"{omegas[0][1]:.4f}" if omegas[0][1] is not None else "n.d.",
        "omega2": f"{omegas[1][0]:.4f}" if omegas[1][0] is not None else "n.d.",
        "prop2":  f"{omegas[1][1]:.4f}" if omegas[1][1] is not None else "n.d.",
        "omega3": f"{omegas[2][0]:.4f}" if omegas[2][0] is not None else "n.d.",
        "prop3":  f"{omegas[2][1]:.4f}" if omegas[2][1] is not None else "n.d.",
        "significant": sig
    })
    print(f"  {gene}: p={pval:.4f}, LRT={lrt:.3f}, ω₃={omegas[2][0]:.4f} "
          f"({float(omegas[2][1])*100:.1f}%)  → {sig}")

# Write CSV
header = ["Gene","n_sequences","n_codons","LRT","p_value",
          "omega1","prop_omega1","omega2","prop_omega2","omega3","prop_omega3","significant_p005"]
with open(OUT_CSV, 'w') as f:
    f.write(",".join(header) + "\n")
    for row in rows:
        line = [str(row.get(
            h.replace("prop_omega","prop").replace("significant_p005","significant"), "n.d."))
            for h in header]
        f.write(",".join(line) + "\n")

print(f"\nTabla guardada: {OUT_CSV}")
print(f"Genes significativos (p < 0.05): {sum(1 for r in rows if r['significant']=='Yes')}/{len(rows)}")
