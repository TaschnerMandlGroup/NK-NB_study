"""
05_network_score_analysis.py
Full network-score analysis from the CellOracle network-analysis tutorial
(Network_analysis_with_Paul_etal_2015_data.html, sections 5-8), applied to
the per-(NK subtype cluster x genetic_group) links fitted in 03_grn_perturb.py.

Standalone script -- reloads `links` and the fitted `oracle` (for its
.obs[CLUSTER]/.obs[GENETIC_GROUP] columns) from the files 03_grn_perturb.py
saves, so it does not need to run in the same process/session as that script.

Notes on API stability: plot_scores_as_rank / plot_score_comparison_2D /
plot_score_per_cluster / plot_score_discributions / plot_network_entropy_distributions
are the names documented in the tutorial notebook (celloracle 0.10.x). Confirm
exact keyword arguments against `help(links.<method>)` in your installed
version before a long batch run -- signatures have drifted across releases.

Refs: Kamimoto et al., Nature 2023 (PMID 36755098).
"""
import os
import matplotlib.pyplot as plt
import celloracle as co

FIG_DIR = "/home/sara_wz/bioinf_isilon/Research/HALBRITTER/zHalbritter_TaschnerMandl/wernig-zorc.sara/out/NK_project/figures/network_score"
os.makedirs(FIG_DIR, exist_ok=True)

# Most of the plotting calls below (plot_degree_distributions,
# plot_scores_as_rank, plot_score_per_cluster, plot_score_discributions,
# plot_network_entropy_distributions) let CellOracle create + save its own
# figure internally rather than us passing an ax we control. With 12
# cluster_group units and long labels (e.g. "hNK_Bm1_ATRXwtMYCNwt"),
# CellOracle's default figure size crushes/overlaps everything. Raising
# matplotlib's default figure size and font here widens every figure
# CellOracle creates for the rest of this script.
plt.rcParams["figure.figsize"] = (16, 9)
plt.rcParams["figure.dpi"] = 110
plt.rcParams["font.size"] = 9
plt.rcParams["savefig.bbox"] = "tight"

LINKS_PATH  = "links_per_group.celloracle.links"
ORACLE_PATH = "oracle_fitted.celloracle.oracle"
CLUSTER       = "predicted_NK_type_de.n350"
GENETIC_GROUP = "genetic_group"

links  = co.load_hdf5(LINKS_PATH)
oracle = co.load_hdf5(ORACLE_PATH)
adata  = oracle.adata

# TFs of interest already perturbed in 03_grn_perturb.py / 03a_top_TFs.py
GOI = ["TBX21", "RUNX3", "MEF2C", "KLF2", "POU2F2", "PRDM1", "FOS"]
REFERENCE_GROUP = "ATRXdel"
UNITS = list(links.filtered_links.keys())   # e.g. "NK1_ATRXdel", "NK1_control", ...

# =============================================================================
# 1) degree distributions per unit, with power-law fit overlay
# =============================================================================
os.makedirs(f"{FIG_DIR}/degree_distribution", exist_ok=True)
links.plot_degree_distributions(plot_model=True,
                                 save=f"{FIG_DIR}/degree_distribution/")
plt.close("all")
print("[network score] wrote degree distributions")

# =============================================================================
# 2) network score calculation (degree/betweenness/eigenvector centrality)
# =============================================================================
links.get_network_score()
links.merged_score.to_csv(f"{FIG_DIR}/network_score_merged.csv")
print(f"[network score] wrote network_score_merged.csv  units={len(UNITS)}")

# =============================================================================
# 3) ranked bar chart of top genes by centrality, per unit
# =============================================================================
os.makedirs(f"{FIG_DIR}/ranked_score", exist_ok=True)
for unit in UNITS:
    try:
        links.plot_scores_as_rank(cluster=unit, n_gene=30,
                                   save=f"{FIG_DIR}/ranked_score/{unit}")
    except Exception as e:
        print(f"[network score] plot_scores_as_rank failed for {unit}: {e}")
    finally:
        plt.close("all")

# =============================================================================
# 4) 2D score comparison: ATRXdel vs each other group, within matched cluster
#    (same cluster-matching logic as the edge/network-score diff in
#    03_grn_perturb.py, but using CellOracle's own built-in comparison plot)
# =============================================================================
os.makedirs(f"{FIG_DIR}/score_comparison", exist_ok=True)
clusters_present = sorted(adata.obs[CLUSTER].unique())
other_groups = [g for g in sorted(adata.obs[GENETIC_GROUP].unique()) if g != REFERENCE_GROUP]
SCORE_VALUES = ["degree_centrality_all", "betweenness_centrality", "eigenvector_centrality"]

for cl in clusters_present:
    unit_ref = f"{cl}_{REFERENCE_GROUP}"
    if unit_ref not in UNITS:
        continue
    for grp in other_groups:
        unit_other = f"{cl}_{grp}"
        if unit_other not in UNITS:
            continue
        for val in SCORE_VALUES:
            try:
                links.plot_score_comparison_2D(
                    value=val, cluster1=unit_ref, cluster2=unit_other,
                    percentile=98,
                    save=f"{FIG_DIR}/score_comparison/{cl}_{val}_{REFERENCE_GROUP}-vs-{grp}",
                )
            except Exception as e:
                print(f"[network score] plot_score_comparison_2D failed "
                      f"({unit_ref} vs {unit_other}, {val}): {e}")
            finally:
                plt.close("all")

# =============================================================================
# 5) network score dynamics for genes of interest, across all units
# =============================================================================
os.makedirs(f"{FIG_DIR}/score_per_gene", exist_ok=True)
for gene in GOI:
    try:
        links.plot_score_per_cluster(goi=gene, save=f"{FIG_DIR}/score_per_gene/{gene}")
    except Exception as e:
        print(f"[network score] plot_score_per_cluster failed for {gene}: {e}")
    finally:
        plt.close("all")

# =============================================================================
# 6) score distributions across all units (boxplots)
# =============================================================================
# plt.close("all") both before and after: these two calls (6 and 7) previously
# ran back-to-back with no figure cleanup in between, and produced visually
# identical output (7 silently reused/half-overwrote 6's still-open figure
# instead of drawing its own) -- closing figures immediately before each call
# guarantees it starts from a clean matplotlib state.
os.makedirs(f"{FIG_DIR}/score_distribution", exist_ok=True)
plt.close("all")
links.plot_score_discributions(values=["degree_centrality_all", "eigenvector_centrality"],
                                method="boxplot",
                                save=f"{FIG_DIR}/score_distribution/")
plt.close("all")

# =============================================================================
# 7) network entropy per unit (boxplot)
# =============================================================================
os.makedirs(f"{FIG_DIR}/network_entropy", exist_ok=True)
print("[network score] merged_score columns:", list(links.merged_score.columns))
plt.close("all")
try:
    links.plot_network_entropy_distributions(save=f"{FIG_DIR}/network_entropy/")
except Exception as e:
    print(f"[network score] plot_network_entropy_distributions failed: {e}")
finally:
    plt.close("all")

print("[network score] done -- figures under", FIG_DIR)
