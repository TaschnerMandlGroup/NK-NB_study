"""
07_development_comparison.py
Compares each TF perturbation's simulated vector field against the NK cells'
own developmental (pseudotime) vector field, and computes the CellOracle
Perturbation Score (PS) -- from the simulation tutorial
(Gata1_KO_simulation_with_Paul_etal_2015_data.html), sections 5-6.

Pseudotime: loaded from tables/ALL_NK_TSCAN_pseudotime.csv, exported by the
new chunk added to 20_20250227_NEW_Trajectory_analysis.Rmd right after the
combined-ALL_NK TSCAN block (quickPseudotime(sce, use.dimred="UMAP",
start="hNK_Bm2") -> averagePseudotime()). ALL_NK is all four genetic groups
combined -- the same population as nk_rna.h5ad -- so cell barcodes should
match directly; the merge below reports how many cells matched as a sanity
check. This replaces an earlier scanpy-dpt placeholder: TSCAN's cluster-MST
ordering is both a more principled trajectory method for discrete NK
subtypes and already validated in your own analysis, rather than a fresh
guess at a root cell from this script.

Standalone script -- reloads the fitted `oracle` from the file
03_grn_perturb.py saves (oracle.to_hdf5), and recomputes the TBX21/RUNX3
90th-percentile "high" values used for the OE perturbations from the
imputed baseline, so it does not need to run in the same process/session as
that script (or as 06_simulation_vector_field.py).

Refs: Kamimoto et al., Nature 2023 (PMID 36755098).
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co
import matplotlib.pyplot as plt
from celloracle.applications import Gradient_calculator, Oracle_development_module
from celloracle.visualizations.config import CONFIG

FIG_DIR = "/home/sara_wz/bioinf_isilon/Research/HALBRITTER/zHalbritter_TaschnerMandl/wernig-zorc.sara/out/NK_project/figures/development_comparison"
os.makedirs(FIG_DIR, exist_ok=True)

ORACLE_PATH = "oracle_fitted.celloracle.oracle"
PSEUDOTIME_CSV = "/home/sara_wz/bioinf_isilon/Research/HALBRITTER/zHalbritter_TaschnerMandl/wernig-zorc.sara/out/NK_project/tables/ALL_NK_TSCAN_pseudotime.csv"  # host path -- /home/rstudio/mnt_out/... is the Docker-internal mount for the same file
CLUSTER       = "predicted_NK_type_de.n350"
N_PROP = 3
LINEAGE_CLUSTERS = ["hNK_Bm2", "hNK_Bm3", "hNK_Bm4", "hNK_Bm1"]  # developmental path: hNK_Bm2 -> hNK_Bm3 -> hNK_Bm4 -> hNK_Bm1
N_GRID = 40
N_NEIGHBORS = 200
# see SCALE_QUIVER/SCALE_GRID note in 06_simulation_vector_field.py -- the
# "Perturb simulation" panel used a hardcoded scale_for_simulation=0.5 tuned
# to the tutorial's own embedding, not ours, which is what made it look like
# chaotic long crossing arrows. scale_for_pseudotime is left at 40 since the
# "Development flow" panel it controls already renders cleanly.
SCALE_SIMULATION = None

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

# =============================================================================
# 1) pseudotime: load TSCAN's combined-ALL_NK ordering (skip if
#    oracle.adata.obs["Pseudotime"] was already set some other way)
# =============================================================================
if "Pseudotime" not in oracle.adata.obs.columns:
    pt = pd.read_csv(PSEUDOTIME_CSV, index_col="cell_id")["pseudotime"]

    n_matched = oracle.adata.obs_names.isin(pt.index).sum()
    print(f"[pseudotime] {n_matched}/{oracle.adata.n_obs} oracle cells matched "
          f"by barcode in {PSEUDOTIME_CSV}")
    assert n_matched > 0, (
        "No cell barcodes matched between oracle.adata and the TSCAN pseudotime "
        "CSV -- check that ALL_NK in 20_20250227_NEW_Trajectory_analysis.Rmd is "
        "the same object 01_export_inputs.R exported to nk_rna.h5ad."
    )
    if n_matched < oracle.adata.n_obs:
        print(f"[pseudotime] WARNING: {oracle.adata.n_obs - n_matched} cells have "
              f"no TSCAN pseudotime (dropped from the MST as NA, or barcode "
              f"mismatch) -- these cells will get NaN and be excluded downstream.")

    oracle.adata.obs["Pseudotime"] = oracle.adata.obs_names.map(pt)

fig, ax = plt.subplots(figsize=(7, 7))
sc.pl.embedding(adata=oracle.adata, basis=oracle.embedding_name, ax=ax,
                 cmap="rainbow", color=["Pseudotime"], show=False)
fig.savefig(f"{FIG_DIR}/pseudotime.pdf"); plt.close(fig)

# =============================================================================
# 2) developmental gradient field (built once, reused for every perturbation)
# =============================================================================
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=N_GRID, n_neighbors=N_NEIGHBORS)

gradient.calculate_mass_filter(min_mass=52, plot=True)  # see DEFAULT_MIN_MASS note in 06_simulation_vector_field.py
plt.gcf().savefig(f"{FIG_DIR}/gradient_mass_filter.pdf"); plt.close("all")

gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly": 3}, plot=True)
plt.gcf().savefig(f"{FIG_DIR}/pseudotime_on_grid.pdf"); plt.close("all")

gradient.calculate_gradient()

# visualize_results() creates its own figure (no ax= kwarg in this CellOracle
# version) -- same pattern as calculate_mass_filter/transfer_data_into_grid above.
# scale=None: see SCALE_QUIVER/SCALE_GRID note in 06_simulation_vector_field.py.
gradient.visualize_results(scale=None, s=5)
plt.gcf().savefig(f"{FIG_DIR}/developmental_flow.pdf"); plt.close("all")

fig, ax = plt.subplots(figsize=(7, 7))
gradient.plot_dev_flow_on_grid(scale=None, ax=ax)
fig.savefig(f"{FIG_DIR}/developmental_flow_on_grid.pdf"); plt.close(fig)

print("[development comparison] wrote pseudotime + gradient figures")

# =============================================================================
# 3) Perturbation Score (PS): does each KO/OE promote (+) or block (-) the
#    inferred NK developmental trajectory, cell by cell / on the grid?
# =============================================================================
for tag, cond in PERTURBATIONS.items():
    oracle.simulate_shift(perturb_condition=cond, n_propagation=N_PROP)
    oracle.estimate_transition_prob(n_neighbors=N_NEIGHBORS, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    dev = Oracle_development_module()
    dev.load_differentiation_reference_data(gradient_object=gradient)
    dev.load_perturb_simulation_data(oracle_object=oracle)
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=10)

    out_dir = f"{FIG_DIR}/{tag}"
    os.makedirs(out_dir, exist_ok=True)

    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    dev.plot_inner_product_on_grid(vm=0.02, s=50, ax=ax[0]); ax[0].set_title(f"{tag}: PS")
    dev.plot_inner_product_random_on_grid(vm=0.02, s=50, ax=ax[1]); ax[1].set_title(f"{tag}: PS (random control)")
    fig.savefig(f"{out_dir}/perturbation_score.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 7))
    dev.plot_inner_product_on_grid(vm=0.02, s=50, ax=ax)
    dev.plot_simulation_flow_on_grid(scale=SCALE_SIMULATION, show_background=False, ax=ax)
    fig.savefig(f"{out_dir}/perturbation_score_with_flow.pdf"); plt.close(fig)

    # visualize_development_module_layout_0() returns None in this CellOracle
    # version (draws internally, like visualize_results() above) -- grab the
    # figure it left open instead of using a return value.
    plt.close("all")
    dev.visualize_development_module_layout_0(
        s=5, scale_for_simulation=SCALE_SIMULATION, s_grid=50, scale_for_pseudotime=40, vm=0.02)
    plt.gcf().savefig(f"{out_dir}/development_module_layout.pdf"); plt.close("all")

    print(f"[development comparison] wrote PS figures under {out_dir}")

# =============================================================================
# 4) lineage-restricted view: PS computed only within LINEAGE_CLUSTERS
# =============================================================================
cell_idx_lineage = np.flatnonzero(oracle.adata.obs[CLUSTER].isin(LINEAGE_CLUSTERS).values)
assert len(cell_idx_lineage) > 0, f"No cells found in LINEAGE_CLUSTERS={LINEAGE_CLUSTERS}"

for tag, cond in {"TBX21_KO": {"TBX21": 0.0}, "RUNX3_KO": {"RUNX3": 0.0}}.items():
    oracle.simulate_shift(perturb_condition=cond, n_propagation=N_PROP)
    oracle.estimate_transition_prob(n_neighbors=N_NEIGHBORS, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    dev = Oracle_development_module()
    dev.load_differentiation_reference_data(gradient_object=gradient)
    dev.load_perturb_simulation_data(
        oracle_object=oracle, cell_idx_use=cell_idx_lineage,
        name=f"Lineage_{'_'.join(LINEAGE_CLUSTERS)}",
    )
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=10)

    out_dir = f"{FIG_DIR}/{tag}_lineage"
    os.makedirs(out_dir, exist_ok=True)

    plt.close("all")
    dev.visualize_development_module_layout_0(
        s=5, scale_for_simulation=SCALE_SIMULATION, s_grid=50, scale_for_pseudotime=40, vm=0.02)
    plt.gcf().savefig(f"{out_dir}/development_module_layout.pdf"); plt.close("all")

    # coolwarm variant, as shown in the tutorial's final section
    CONFIG["cmap_ps"] = "coolwarm"
    dev.visualize_development_module_layout_0(
        s=5, scale_for_simulation=SCALE_SIMULATION, s_grid=50, scale_for_pseudotime=40, vm=0.02)
    plt.gcf().savefig(f"{out_dir}/development_module_layout_coolwarm.pdf"); plt.close("all")
    CONFIG["cmap_ps"] = "PiYG"  # CellOracle's actual built-in default (not the literal string "default")

    print(f"[development comparison] wrote lineage-restricted figures under {out_dir}")

print("[development comparison] done -- figures under", FIG_DIR)
