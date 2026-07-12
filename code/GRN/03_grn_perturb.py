"""
03_grn_perturb.py
Fit ONE GRN per (NK subtype cluster x genetic group) from the shared base GRN,
then run an in silico TBX21 knockout and export the per-cell expression shift
for the cytotoxicity module. cluster_group is the GRN unit; PATIENT is the
inferential unit (handled in 04).

Design notes
  * Network unit = NK subtype cluster x genetic group (cluster_group), so that
    per-cluster GRNs can be compared ATRXdel-vs-other-groups within a matched
    NK subtype instead of pooling all subtypes into one group-level network.
    CLUSTER ("predicted_NK_type_de.n350") is the scRNA-side subtype label; it
    is independent of the ATAC-side de novo clustering (e.g. "Clusters_res1")
    used only when the base GRN itself is built from an ATAC subset.
  * Any cluster_group combination with < MIN_CELLS_PER_UNIT cells is dropped
    before fitting -- CellOracle recommends >=50 cells per GRN unit (Kamimoto
    et al., Nature 2023, Extended Data Fig. 3), and splitting the original
    genetic-group units (403-961 cells) further by subtype can easily push
    some combinations below that floor.
  * The KO is most informative FROM A HIGH-TBX21 BASELINE (control): does removing
    TBX21 reproduce the cytotoxicity-low state seen in ATRXdel? In ATRXmut, TBX21
    is already low, so its KO delta should be small -- that asymmetry IS the result.
  * The simulated delta is a model-internal extrapolation from the inferred GRN;
    it prioritizes/orders hypotheses, it does not prove causality. Frame results
    as "consistent with a causal role."
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co
import matplotlib.pyplot as plt

FIG_DIR = "/home/sara_wz/bioinf_isilon/Research/HALBRITTER/zHalbritter_TaschnerMandl/wernig-zorc.sara/out/NK_project/figures"
os.makedirs(FIG_DIR, exist_ok=True)

ADATA    = "/home/sara_wz/bioinf_isilon/Research/HALBRITTER/zHalbritter_TaschnerMandl/wernig-zorc.sara/out/NK_project/h5ad/nk_rna.h5ad"
BASEGRN  = "/home/sara_wz/bioinf_isilon/Research/HALBRITTER/zHalbritter_TaschnerMandl/wernig-zorc.sara/out/NK_project/GRN/base_GRN_hg38.parquet"
GENETIC_GROUP = "genetic_group"        # obs col -> ATRXdel / control / MYCNamp / ATRXwtMYCNwt
CLUSTER       = "predicted_NK_type_de.n350"  # obs col -> NK subtype (RNA-side label transferred to ATAC as predicted_NK_subtype)
GROUP    = "cluster_group"     # obs col -> GRN unit (cluster x genetic group)
PATIENT  = "patient_id"        # obs col -> replicate (exported for script 04)
EMB      = "X_umap"            # obsm key for the embedding
TF       = "TBX21"
KO_VALUE = 0.0                 # 0 = knockout
N_PROP   = 3                   # signal-propagation steps
ALPHA    = 10                  # ridge penalty for GRN fit
N_HVG    = 3000
MIN_CELLS_PER_UNIT = 50        # CellOracle's own floor for reliable GRN inference
OUT_DELTA = "perturb_deltas.csv"

# cytotoxicity / effector module 
CYTO = ["GZMB","PRF1","IFNG","GNLY","NKG7","KLRD1","GZMA","GZMH","KLRF1","FCGR3A"]
KEEP = list(set(CYTO + [TF, "RUNX3", "CD6"]))   # force-retain genes of interest

# read data
adata = sc.read_h5ad(ADATA)
adata.layers["raw_count"] = adata.X.copy()      # preserve raw counts

# ---- build the cluster x group GRN unit, drop underpowered combinations ---
assert CLUSTER in adata.obs.columns, \
    f"{CLUSTER} not found in adata.obs -- confirm it survived DietSeurat/sceasy export in 01_export_inputs.R"
adata.obs[GENETIC_GROUP] = adata.obs[GENETIC_GROUP].astype(str)
adata.obs[CLUSTER]       = adata.obs[CLUSTER].astype(str)
adata.obs[GROUP] = adata.obs[CLUSTER] + "_" + adata.obs[GENETIC_GROUP]

unit_counts = adata.obs[GROUP].value_counts()
print("[cluster_group] cells per (cluster x genetic_group) unit:")
print(unit_counts.to_string())

small_units = unit_counts[unit_counts < MIN_CELLS_PER_UNIT].index.tolist()
if small_units:
    print(f"[cluster_group] dropping {len(small_units)} unit(s) below "
          f"MIN_CELLS_PER_UNIT={MIN_CELLS_PER_UNIT}: {small_units}")
    adata = adata[~adata.obs[GROUP].isin(small_units)].copy()

# standard scanpy preprocessing (for HVG selection + embedding only)
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG)
for g in KEEP:
    if g in adata.var_names:
        adata.var.loc[g, "highly_variable"] = True
adata = adata[:, adata.var.highly_variable].copy()

# CellOracle re-imports RAW counts; restore them on the HVG-subset object
adata.X = adata.layers["raw_count"][:, adata.var["highly_variable"].values] \
    if adata.layers["raw_count"].shape[1] != adata.n_vars else adata.layers["raw_count"]
# (if the line above is awkward for your matrix type, simply re-subset raw counts:
#  adata.X = adata[:, adata.var_names].layers["raw_count"])

# ---- Oracle object ---------------------------------------------------------
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(
    adata=adata,
    cluster_column_name=GROUP,      # <-- cluster x genetic group is the GRN unit
    embedding_name=EMB,
)
base_GRN = pd.read_parquet(BASEGRN)
oracle.import_TF_data(TF_info_matrix=base_GRN)

# imputation (KNN over PCA). If donor batch structure remains, use a
# batch-corrected embedding upstream so this reflects state, not patient.
oracle.perform_PCA()

cum_var = np.cumsum(oracle.pca.explained_variance_ratio_) if hasattr(oracle, "pca") else None

if cum_var is not None:
    n_comps = int(np.searchsorted(cum_var, 0.90) + 1)   # first n_comps reaching 90% variance
    n_comps = min(n_comps, len(cum_var))                # can't exceed computed components
else:
    n_comps = 50

n_comps = min(max(n_comps, 20), 50)

# tutorial figure: cumulative explained-variance-ratio plot w/ chosen PC threshold
if cum_var is not None:
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(np.arange(1, len(cum_var) + 1), cum_var, marker="o", markersize=2)
    ax.axvline(n_comps, linestyle="--", color="grey")
    ax.axhline(0.90, linestyle="--", color="grey")
    ax.set_xlabel("# PCs"); ax.set_ylabel("cumulative explained variance ratio")
    ax.set_title(f"PCA variance (chosen n_comps={n_comps})")
    fig.savefig(f"{FIG_DIR}/pca_cumulative_variance.pdf"); plt.close(fig)

k = int(0.025 * oracle.adata.n_obs)
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True,b_sight=k*8, b_maxl=k*4, n_jobs=20)

# ---- per (cluster x group) GRN + fit for simulation ------------------------
# links.filtered_links is now keyed by "<cluster>_<genetic_group>", e.g. "NK1_ATRXdel"
links = oracle.get_links(cluster_name_for_GRN_unit=GROUP, alpha=ALPHA,verbose_level=1, test_mode=False)
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
links.to_hdf5(file_path="links_per_group.celloracle.links")

oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=ALPHA, use_cluster_specific_TFdict=True)

# confirm TBX21 is actually a regulator in the fitted model
assert TF in oracle.active_regulatory_genes, \
    f"{TF} is not an active regulator; check base GRN motif retention (script 02)."

# persist the fitted Oracle object so 05_/06_/07_ can reload it in a fresh
# process instead of requiring the same interactive session as this script
ORACLE_PATH = "oracle_fitted.celloracle.oracle"
oracle.to_hdf5(file_path=ORACLE_PATH)
print(f"[oracle] saved fitted Oracle object to {ORACLE_PATH}")

# ---- in silico TBX21 knockout ---------------------------------------------
oracle.simulate_shift(perturb_condition={TF: KO_VALUE}, n_propagation=N_PROP)

# per-cell delta = simulated - imputed (baseline), in normalized space
sim  = oracle.adata.layers["simulated_count"]
base = oracle.adata.layers["imputed_count"]
delta = np.asarray(sim - base)

genes = list(oracle.adata.var_names)
cyto_in = [g for g in CYTO if g in genes]
idx = [genes.index(g) for g in cyto_in]

df = pd.DataFrame(delta[:, idx], columns=cyto_in, index=oracle.adata.obs_names)
df.insert(0, "group",   oracle.adata.obs[GENETIC_GROUP].values)
df.insert(1, "cluster", oracle.adata.obs[CLUSTER].values)
df.insert(2, "patient", oracle.adata.obs[PATIENT].values)
df.to_csv(OUT_DELTA)
print(f"[perturb] wrote {OUT_DELTA}  cells={df.shape[0]}  cyto genes={len(cyto_in)}")

# complementary test: TBX21 over-expression "rescue" in ATRXdel 
OUT_DELTA_OE = "perturb_deltas_OE.csv"

tbx21_hi = np.quantile(base[:, genes.index(TF)], 0.9)
oracle.simulate_shift(perturb_condition={TF: tbx21_hi}, n_propagation=N_PROP)

sim_oe   = oracle.adata.layers["simulated_count"]
delta_oe = np.asarray(sim_oe - base)

df_oe = pd.DataFrame(delta_oe[:, idx], columns=cyto_in, index=oracle.adata.obs_names)
df_oe.insert(0, "group",   oracle.adata.obs[GENETIC_GROUP].values)
df_oe.insert(1, "cluster", oracle.adata.obs[CLUSTER].values)
df_oe.insert(2, "patient", oracle.adata.obs[PATIENT].values)
df_oe.to_csv(OUT_DELTA_OE)
print(f"[perturb-OE] wrote {OUT_DELTA_OE}  cells={df_oe.shape[0]}  cyto genes={len(cyto_in)}")

########

TF2 = "RUNX3"
assert TF2 in oracle.active_regulatory_genes, \
    f"{TF2} is not an active regulator; check base GRN motif retention (script 02)."

genes = list(oracle.adata.var_names)
cyto_in = [g for g in CYTO if g in genes]
idx = [genes.index(g) for g in cyto_in]
base = oracle.adata.layers["imputed_count"]   # re-grab if not already in scope

# ---- RUNX3 knockout ----
oracle.simulate_shift(perturb_condition={TF2: 0.0}, n_propagation=N_PROP)
sim_ko2 = oracle.adata.layers["simulated_count"]
delta_ko2 = np.asarray(sim_ko2 - base)
df_ko2 = pd.DataFrame(delta_ko2[:, idx], columns=cyto_in, index=oracle.adata.obs_names)
df_ko2.insert(0, "group",   oracle.adata.obs[GENETIC_GROUP].values)
df_ko2.insert(1, "cluster", oracle.adata.obs[CLUSTER].values)
df_ko2.insert(2, "patient", oracle.adata.obs[PATIENT].values)
df_ko2.to_csv("perturb_deltas_RUNX3_KO.csv")
print(f"[perturb-RUNX3-KO] wrote perturb_deltas_RUNX3_KO.csv  cells={df_ko2.shape[0]}")

# ---- RUNX3 over-expression ----
runx3_hi = np.quantile(base[:, genes.index(TF2)], 0.9)
oracle.simulate_shift(perturb_condition={TF2: runx3_hi}, n_propagation=N_PROP)
sim_oe2 = oracle.adata.layers["simulated_count"]
delta_oe2 = np.asarray(sim_oe2 - base)
df_oe2 = pd.DataFrame(delta_oe2[:, idx], columns=cyto_in, index=oracle.adata.obs_names)
df_oe2.insert(0, "group",   oracle.adata.obs[GENETIC_GROUP].values)
df_oe2.insert(1, "cluster", oracle.adata.obs[CLUSTER].values)
df_oe2.insert(2, "patient", oracle.adata.obs[PATIENT].values)
df_oe2.to_csv("perturb_deltas_RUNX3_OE.csv")
print(f"[perturb-RUNX3-OE] wrote perturb_deltas_RUNX3_OE.csv  cells={df_oe2.shape[0]}")


# how does RUNX3 KO/OE effect TBX21 expression? (sanity check: TBX21 is a known RUNX3 target)
tf_idx = genes.index("TBX21")

# RUNX3 KO -> effect on TBX21 expression
delta_ko2_tbx21 = np.asarray(sim_ko2 - base)[:, tf_idx]
df_ko2_tbx21 = pd.DataFrame({"TBX21": delta_ko2_tbx21}, index=oracle.adata.obs_names)
df_ko2_tbx21.insert(0, "group",   oracle.adata.obs[GENETIC_GROUP].values)
df_ko2_tbx21.insert(1, "cluster", oracle.adata.obs[CLUSTER].values)
df_ko2_tbx21.insert(2, "patient", oracle.adata.obs[PATIENT].values)
df_ko2_tbx21.to_csv("perturb_deltas_RUNX3_KO_TBX21.csv")
print(f"[RUNX3-KO->TBX21] wrote perturb_deltas_RUNX3_KO_TBX21.csv  cells={df_ko2_tbx21.shape[0]}")

# RUNX3 OE -> effect on TBX21 expression
delta_oe2_tbx21 = np.asarray(sim_oe2 - base)[:, tf_idx]
df_oe2_tbx21 = pd.DataFrame({"TBX21": delta_oe2_tbx21}, index=oracle.adata.obs_names)
df_oe2_tbx21.insert(0, "group",   oracle.adata.obs[GENETIC_GROUP].values)
df_oe2_tbx21.insert(1, "cluster", oracle.adata.obs[CLUSTER].values)
df_oe2_tbx21.insert(2, "patient", oracle.adata.obs[PATIENT].values)
df_oe2_tbx21.to_csv("perturb_deltas_RUNX3_OE_TBX21.csv")
print(f"[RUNX3-OE->TBX21] wrote perturb_deltas_RUNX3_OE_TBX21.csv  cells={df_oe2_tbx21.shape[0]}")

# =============================================================================
# Compare the fitted GRNs: ATRXdel vs every other genetic group, within each
# matched NK subtype cluster. Two complementary views:
#   (a) edge-level: which specific TF->target edges differ in strength/presence
#   (b) network-level: which genes' regulatory importance (centrality) shifts
# Units missing from links.filtered_links were dropped upstream for having
# fewer than MIN_CELLS_PER_UNIT cells -- those cluster/group combinations are
# silently skipped here rather than compared.
# =============================================================================
REFERENCE_GROUP = "ATRXdel"
clusters_present = sorted(adata.obs[CLUSTER].unique())
other_groups = [g for g in sorted(adata.obs[GENETIC_GROUP].unique()) if g != REFERENCE_GROUP]

# ---- (a) edge-level differences ---------------------------------------------
edge_diff_rows = []
for cl in clusters_present:
    unit_ref = f"{cl}_{REFERENCE_GROUP}"
    if unit_ref not in links.filtered_links:
        continue
    ref_edges = links.filtered_links[unit_ref][["source", "target", "coef_abs", "coef_mean"]]
    for grp in other_groups:
        unit_other = f"{cl}_{grp}"
        if unit_other not in links.filtered_links:
            continue
        other_edges = links.filtered_links[unit_other][["source", "target", "coef_abs", "coef_mean"]]
        merged = ref_edges.merge(
            other_edges, on=["source", "target"], how="outer",
            suffixes=(f"_{REFERENCE_GROUP}", f"_{grp}"),
        ).fillna(0)  # edge absent in one group's filtered network -> treated as coef 0
        merged.insert(0, "cluster", cl)
        merged.insert(1, "group_compared", grp)
        merged["delta_coef_abs"] = merged[f"coef_abs_{REFERENCE_GROUP}"] - merged[f"coef_abs_{grp}"]
        edge_diff_rows.append(merged)

if edge_diff_rows:
    edge_diff = pd.concat(edge_diff_rows, ignore_index=True)
    edge_diff = edge_diff.reindex(edge_diff["delta_coef_abs"].abs().sort_values(ascending=False).index)
    edge_diff.to_csv("grn_edge_diff_ATRXdel_vs_other_groups.csv", index=False)
    print(f"[edge diff] wrote grn_edge_diff_ATRXdel_vs_other_groups.csv "
          f"rows={edge_diff.shape[0]} (sorted by |delta_coef_abs|, largest first)")
else:
    print("[edge diff] no cluster had both ATRXdel and a comparison group above MIN_CELLS_PER_UNIT")

# ---- (b) network-score (centrality) differences ------------------------------
links.get_network_score()
# merged_score is gene-indexed per CellOracle docs; pull the gene name out as a column
merged_score = links.merged_score.reset_index().rename(columns={"index": "gene"})
score_cols = ["degree_centrality_all", "degree_centrality_in", "degree_centrality_out",
              "betweenness_centrality", "eigenvector_centrality"]
present_units = set(merged_score["cluster"])

score_diff_rows = []
for cl in clusters_present:
    unit_ref = f"{cl}_{REFERENCE_GROUP}"
    if unit_ref not in present_units:
        continue
    ref_score = merged_score[merged_score["cluster"] == unit_ref][["gene"] + score_cols]
    for grp in other_groups:
        unit_other = f"{cl}_{grp}"
        if unit_other not in present_units:
            continue
        other_score = merged_score[merged_score["cluster"] == unit_other][["gene"] + score_cols]
        merged = ref_score.merge(other_score, on="gene", suffixes=(f"_{REFERENCE_GROUP}", f"_{grp}"))
        for col in score_cols:
            merged[f"delta_{col}"] = merged[f"{col}_{REFERENCE_GROUP}"] - merged[f"{col}_{grp}"]
        merged.insert(0, "cluster", cl)
        merged.insert(1, "group_compared", grp)
        score_diff_rows.append(merged)

if score_diff_rows:
    score_diff = pd.concat(score_diff_rows, ignore_index=True)
    score_diff.to_csv("grn_network_score_diff_ATRXdel_vs_other_groups.csv", index=False)
    print(f"[network score diff] wrote grn_network_score_diff_ATRXdel_vs_other_groups.csv "
          f"rows={score_diff.shape[0]}")
else:
    print("[network score diff] no cluster had both ATRXdel and a comparison group above MIN_CELLS_PER_UNIT")