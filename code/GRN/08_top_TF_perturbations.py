"""
08_top_TF_perturbations.py
Standalone replacement for 03a_top_TFs.py, updated for the per-(NK subtype
cluster x genetic_group) GRN units fitted in 03_grn_perturb.py.

Why a new script instead of patching 03a_top_TFs.py in place:
  * 03a assumed `links.filtered_links` was keyed by plain genetic group
    (e.g. "ATRXdel"). It's now keyed by "<cluster>_<genetic_group>"
    (e.g. "hNK_Bm2_ATRXdel"), so `links.filtered_links["ATRXdel"]` in 03a's
    "ATRXdel specific TFs" section would raise a KeyError as-is.
  * 03a relied on `oracle`/`links`/`GROUP`/`PATIENT`/`N_PROP` already being
    in scope from an interactive 03_grn_perturb.py session (same issue
    05_/06_/07_ hit before being made standalone).
  * 03a had accumulated duplication from iterative editing:
      - `links = co.load_hdf5(...)` loaded twice (top of file, and again
        before the CD6-pathway section)
      - `genes = list(oracle.adata.var_names)`, `cyto_in`, `idx`, and
        `base = oracle.adata.layers["imputed_count"]` recomputed 3-4x
        identically (none of these change between simulate_shift() calls --
        only the "simulated_count" layer changes)
      - four near-identical per-perturbation export helpers
        (export_module_and_readouts / export_module_and_runx3_only /
        export_module / export_module_and_tbx21) that differ only in which
        extra single-gene readouts get written alongside the cytotoxicity
        module -- consolidated into export_perturbation()/perturb_ko_oe() in
        grn_perturb_utils.py (shared with 10_additional_TF_perturbations.py).
  * Every output CSV here carries "group", "cluster", "patient" columns
    (matching 03_grn_perturb.py's current convention), not the old
    2-column "group"/"patient" from 03a.

Run any time after 03_grn_perturb.py (reloads its saved oracle/links).
"""
import celloracle as co
import pandas as pd
import numpy as np

from grn_perturb_utils import make_perturb_helpers

ORACLE_PATH = "oracle_fitted.celloracle.oracle"
LINKS_PATH  = "links_per_group.celloracle.links"
GENETIC_GROUP = "genetic_group"
CLUSTER       = "predicted_NK_type_de.n350"
PATIENT       = "patient_id"
N_PROP  = 3
REFERENCE_GROUP = "ATRXdel"

CYTO = ["GZMB", "PRF1", "IFNG", "GNLY", "NKG7", "KLRD1", "GZMA", "GZMH", "KLRF1", "FCGR3A"]

oracle = co.load_hdf5(ORACLE_PATH)
links  = co.load_hdf5(LINKS_PATH)

# computed once -- constant across every perturbation below (only
# oracle.adata.layers["simulated_count"] changes per simulate_shift() call)
genes    = list(oracle.adata.var_names)
cyto_in  = [g for g in CYTO if g in genes]
idx      = [genes.index(g) for g in cyto_in]
base     = oracle.adata.layers["imputed_count"]
tbx21_idx = genes.index("TBX21")
runx3_idx = genes.index("RUNX3")

clusters_present = sorted(oracle.adata.obs[CLUSTER].unique())

export_perturbation, perturb_ko_oe = make_perturb_helpers(
    oracle, genes, cyto_in, idx, base, GENETIC_GROUP, CLUSTER, PATIENT, N_PROP
)

# =============================================================================
# Pass 1: rank TFs by total strength of edges INTO the cytotoxicity genes,
# per (cluster x group) unit
# =============================================================================
rows = []
for unit_name, df in links.filtered_links.items():
    cyto_edges = df[df["target"].isin(CYTO)]
    ranked = (cyto_edges.groupby("source")
              .agg(n_cyto_targets=("target", "nunique"),
                   total_coef_abs=("coef_abs", "sum"))
              .sort_values("total_coef_abs", ascending=False))
    ranked["unit"] = unit_name
    rows.append(ranked.reset_index())
cyto_regulators = pd.concat(rows)
print("Top cytotoxicity-module regulators, summed across all (cluster x group) units:")
print(cyto_regulators.groupby("source")[["n_cyto_targets", "total_coef_abs"]].sum()
      .sort_values("total_coef_abs", ascending=False).head(15))

# =============================================================================
# Pass 2: of those, which also connect strongly to TBX21 or RUNX3?
# =============================================================================
top_tfs = cyto_regulators.groupby("source")["total_coef_abs"].sum().nlargest(15).index.tolist()

for unit_name, df in links.filtered_links.items():
    hits = df[(df["source"].isin(top_tfs)) & (df["target"].isin(["TBX21", "RUNX3"]))]
    if len(hits):
        print(f"\n[{unit_name}] cytotoxicity-regulators also linked to TBX21/RUNX3:")
        print(hits[["source", "target", "coef_mean", "coef_abs", "p"]].to_string(index=False))

# =============================================================================
# MEF2C KO/OE (edges to both TBX21 and RUNX3 per network scan)
# =============================================================================
perturb_ko_oe("MEF2C", "MEF2C", extra_readouts=("TBX21", "RUNX3"))

# =============================================================================
# KLF2 KO/OE (edge to RUNX3 only per network scan)
# =============================================================================
perturb_ko_oe("KLF2", "KLF2", extra_readouts=("RUNX3",))

# =============================================================================
# ATRXdel-specific TFs: rank using ONLY the ATRXdel-labeled units, pooled
# across every NK subtype cluster that has one (fixes 03a's
# links.filtered_links["ATRXdel"] KeyError under the new cluster x group keys)
# =============================================================================
ALREADY_TESTED = {"TBX21", "RUNX3", "MEF2C", "KLF2"}

atrxdel_units = [f"{cl}_{REFERENCE_GROUP}" for cl in clusters_present
                 if f"{cl}_{REFERENCE_GROUP}" in links.filtered_links]
assert atrxdel_units, f"No cluster has a {REFERENCE_GROUP} unit in links.filtered_links"

atrxdel_cyto_edges = pd.concat(
    [links.filtered_links[u][links.filtered_links[u]["target"].isin(CYTO)] for u in atrxdel_units]
)
ranked_atrxdel = (atrxdel_cyto_edges.groupby("source")
                  .agg(n_cyto_targets=("target", "nunique"),
                       total_coef_abs=("coef_abs", "sum"))
                  .sort_values("total_coef_abs", ascending=False))
ranked_atrxdel = ranked_atrxdel[~ranked_atrxdel.index.isin(ALREADY_TESTED)]

print(f"\nTop {REFERENCE_GROUP}-specific cytotoxicity-module regulators "
      f"(pooled across {atrxdel_units}, not yet tested):")
print(ranked_atrxdel.head(10))

top2_new_tfs = ranked_atrxdel.head(2).index.tolist()
print(f"\nSelected for simulation: {top2_new_tfs}")

for tf in top2_new_tfs:
    perturb_ko_oe(tf, tf)

# =============================================================================
# PRDM1 KO/OE -> effect on RUNX3 (RUNX3 is a known PRDM1 target)
# =============================================================================
perturb_ko_oe("PRDM1", "PRDM1", extra_readouts=("RUNX3",))

# =============================================================================
# FOS KO/OE -> effect on TBX21
# =============================================================================
perturb_ko_oe("FOS", "FOS", extra_readouts=("TBX21",))

# =============================================================================
# CD6 signalling pathway TFs (AP-1, NF-kB, NFAT) into TBX21/RUNX3/PRDM1/
# cytotoxicity genes -- edge table across every (cluster x group) unit
# =============================================================================
CD6_PATHWAY_TFS = [
    # AP-1 family (Ras/MAPK/ERK arm)
    "FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND",
    # NF-kB family (PKCtheta arm)
    "RELA", "RELB", "REL", "NFKB1", "NFKB2",
    # NFAT family (calcium influx arm)
    "NFATC1", "NFATC2", "NFATC3", "NFATC4",
]
TARGETS = ["TBX21", "RUNX3", "PRDM1"] + CYTO

rows = []
for unit_name, df in links.filtered_links.items():
    hits = df[(df["source"].isin(CD6_PATHWAY_TFS)) & (df["target"].isin(TARGETS))]
    for _, row in hits.iterrows():
        rows.append({
            "unit": unit_name, "source": row["source"], "target": row["target"],
            "coef_mean": row["coef_mean"], "coef_abs": row["coef_abs"], "p": row["p"],
        })

cd6_edge_table = pd.DataFrame(rows).sort_values(["target", "source", "unit"])
cd6_edge_table.to_csv("cd6_pathway_tf_edges.csv", index=False)
print(f"\nFound {cd6_edge_table.shape[0]} edges from CD6-pathway TFs into "
      f"TBX21/RUNX3/PRDM1/cytotoxicity genes")
print(cd6_edge_table.to_string(index=False))

pivot = cd6_edge_table.pivot_table(index=["source", "target"], columns="unit",
                                    values="coef_mean", aggfunc="first")
print("\n--- Presence/absence by (cluster x group) unit (NaN = no edge) ---")
print(pivot)

print("\n[top TF perturbations] done")
