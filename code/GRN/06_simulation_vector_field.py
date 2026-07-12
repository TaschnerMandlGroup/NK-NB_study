"""
06_simulation_vector_field.py
Vector-field visualization from the CellOracle simulation tutorial
(Gata1_KO_simulation_with_Paul_etal_2015_data.html, section 4), applied to
each TF perturbation already run in 03_grn_perturb.py / 03a_top_TFs.py.

Standalone script -- reloads the fitted `oracle` from the file
03_grn_perturb.py saves (oracle.to_hdf5), and recomputes the TBX21/RUNX3
90th-percentile "high" values used for the OE perturbations from the
imputed baseline, so it does not need to run in the same process/session as
that script.

Every quiver/grid-flow plot below is generated TWICE per perturbation: once
for the actual simulated shift, once for a randomized control (CellOracle's
own recommended sanity check -- a real effect should look structured, the
random control should look like noise).

MIN_MASS: `oracle.suggest_mass_thresholds()` is written to file per
perturbation before `calculate_mass_filter()` runs, so you can override
MIN_MASS_OVERRIDE per tag after inspecting it if 0.01 filters too much/little
of a given perturbation's grid.

Refs: Kamimoto et al., Nature 2023 (PMID 36755098).
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import celloracle as co

FIG_DIR = "/home/sara_wz/bioinf_isilon/Research/HALBRITTER/zHalbritter_TaschnerMandl/wernig-zorc.sara/out/NK_project/figures/vector_field"
os.makedirs(FIG_DIR, exist_ok=True)

ORACLE_PATH = "oracle_fitted.celloracle.oracle"
CLUSTER       = "predicted_NK_type_de.n350"
GENETIC_GROUP = "genetic_group"
N_PROP = 3
N_GRID = 40
N_NEIGHBORS = 200
# scale=None -> matplotlib auto-computes arrow scale from the actual data
# range (the tutorial's literal scale=25/0.5 assumed its own embedding's
# coordinate range, which doesn't match ours -- that's what caused arrows to
# span the whole plot). Override to a fixed number here if auto-scale still
# looks off in either direction.
SCALE_QUIVER = None
SCALE_GRID = None
MIN_MASS_OVERRIDE = {}   # e.g. {"TBX21_KO": 80} to override just one perturbation
# 52 chosen from mass_threshold_suggestions.pdf: smallest value that already
# drops the empty-background grid points while still covering every branch/
# island of the actual cell cloud (gaps start appearing in sparser regions
# above ~100). calculate_p_mass only depends on embedding density, not on
# which TF is perturbed, so one value applies across all perturbations.
DEFAULT_MIN_MASS = 52

oracle = co.load_hdf5(ORACLE_PATH)

# 90th-percentile "high" expression used for the OE perturbations in
# 03_grn_perturb.py -- recomputed here from the reloaded imputed baseline
genes = list(oracle.adata.var_names)
base_imputed = oracle.adata.layers["imputed_count"]
tbx21_hi = np.quantile(base_imputed[:, genes.index("TBX21")], 0.9)
runx3_hi = np.quantile(base_imputed[:, genes.index("RUNX3")], 0.9)

PERTURBATIONS = {
    "TBX21_KO": {"TBX21": 0.0},
    "TBX21_OE": {"TBX21": tbx21_hi},
    "RUNX3_KO": {"RUNX3": 0.0},
    "RUNX3_OE": {"RUNX3": runx3_hi},
}


def plot_vector_field(oracle, tag):
    out_dir = f"{FIG_DIR}/{tag}"
    os.makedirs(out_dir, exist_ok=True)

    # ---- transition probability + embedding shift ---------------------------
    oracle.estimate_transition_prob(n_neighbors=N_NEIGHBORS, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    # ---- quiver: per-cell simulated shift vs randomized control -------------
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    oracle.plot_quiver(scale=SCALE_QUIVER, ax=ax[0]); ax[0].set_title(f"{tag}: simulated shift")
    oracle.plot_quiver_random(scale=SCALE_QUIVER, ax=ax[1]); ax[1].set_title(f"{tag}: randomized control")
    fig.savefig(f"{out_dir}/quiver.pdf"); plt.close(fig)

    # ---- grid mass + threshold selection -------------------------------------
    oracle.calculate_p_mass(smooth=0.8, n_grid=N_GRID, n_neighbors=N_NEIGHBORS)

    fig = plt.figure(figsize=(6, 5))
    oracle.suggest_mass_thresholds(n_suggestion=12)
    plt.gcf().savefig(f"{out_dir}/mass_threshold_suggestions.pdf"); plt.close("all")

    min_mass = MIN_MASS_OVERRIDE.get(tag, DEFAULT_MIN_MASS)
    oracle.calculate_mass_filter(min_mass=min_mass, plot=True)
    plt.gcf().savefig(f"{out_dir}/mass_filter_min-mass-{min_mass}.pdf"); plt.close("all")

    # ---- flow on grid: simulated vs randomized -------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    oracle.plot_simulation_flow_on_grid(scale=SCALE_GRID, ax=ax[0]); ax[0].set_title(f"{tag}: flow")
    oracle.plot_simulation_flow_random_on_grid(scale=SCALE_GRID, ax=ax[1]); ax[1].set_title(f"{tag}: random control")
    fig.savefig(f"{out_dir}/flow_on_grid.pdf"); plt.close(fig)

    # ---- flow overlaid on the (cluster x group) GRN unit used at import -----
    fig, ax = plt.subplots(figsize=(7, 7))
    oracle.plot_cluster_whole(ax=ax, s=10)
    oracle.plot_simulation_flow_on_grid(scale=SCALE_GRID, ax=ax, show_background=False)
    ax.set_title(f"{tag}: flow over cluster_group units")
    fig.savefig(f"{out_dir}/flow_over_cluster_group.pdf"); plt.close(fig)

    # ---- flow overlaid on NK subtype only (readable alternative -- the
    #      cluster_group coloring above has cluster x group categories, which
    #      gets crowded) --------------------------------------------------------
    emb = oracle.adata.obsm[oracle.embedding_name]   # EMB = "X_umap", set at import_anndata_as_raw_count()
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    for cat, color in zip(sorted(oracle.adata.obs[CLUSTER].unique()), plt.cm.tab20.colors):
        m = oracle.adata.obs[CLUSTER].values == cat
        ax[0].scatter(emb[m, 0], emb[m, 1], s=8, color=color, label=cat)
    ax[0].legend(fontsize=6, markerscale=2, title="NK subtype")
    ax[0].set_title(f"{tag}: NK subtype")
    for cat, color in zip(sorted(oracle.adata.obs[GENETIC_GROUP].unique()), plt.cm.Set2.colors):
        m = oracle.adata.obs[GENETIC_GROUP].values == cat
        ax[1].scatter(emb[m, 0], emb[m, 1], s=8, color=color, label=cat)
    ax[1].legend(fontsize=8, markerscale=2, title="genetic group")
    ax[1].set_title(f"{tag}: genetic group")
    for a in ax:
        oracle.plot_simulation_flow_on_grid(scale=SCALE_GRID, ax=a, show_background=False)
    fig.savefig(f"{out_dir}/flow_over_subtype_and_group.pdf"); plt.close(fig)

    print(f"[vector field] wrote figures under {out_dir}")


for tag, cond in PERTURBATIONS.items():
    oracle.simulate_shift(perturb_condition=cond, n_propagation=N_PROP)
    plot_vector_field(oracle, tag)

print("[vector field] done -- figures under", FIG_DIR)
