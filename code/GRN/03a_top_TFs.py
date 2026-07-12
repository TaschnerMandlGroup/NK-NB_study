# top candidates

import celloracle as co
import pandas as pd
import numpy as np

links = co.load_hdf5("links_per_group.celloracle.links")
CYTO = ["GZMB","PRF1","IFNG","GNLY","NKG7","KLRD1","GZMA","GZMH","KLRF1","FCGR3A"]

# ---- Pass 1: rank TFs by total strength of edges INTO the cytotoxicity genes ----
rows = []
for group_name, df in links.filtered_links.items():
    cyto_edges = df[df["target"].isin(CYTO)]
    ranked = (cyto_edges.groupby("source")
              .agg(n_cyto_targets=("target","nunique"),
                   total_coef_abs=("coef_abs","sum"))
              .sort_values("total_coef_abs", ascending=False))
    ranked["group"] = group_name
    rows.append(ranked.reset_index())
cyto_regulators = pd.concat(rows)
print(cyto_regulators.groupby("source")[["n_cyto_targets","total_coef_abs"]].sum()
      .sort_values("total_coef_abs", ascending=False).head(15))

# ---- Pass 2: of those, which also connect strongly to TBX21 or RUNX3? ----
top_tfs = cyto_regulators.groupby("source")["total_coef_abs"].sum().nlargest(15).index.tolist()

for group_name, df in links.filtered_links.items():
    hits = df[(df["source"].isin(top_tfs)) & (df["target"].isin(["TBX21","RUNX3"]))]
    if len(hits):
        print(f"\n[{group_name}] cytotoxicity-regulators also linked to TBX21/RUNX3:")
        print(hits[["source","target","coef_mean","coef_abs","p"]].to_string(index=False))
        

# preturb top TFs

# ============================================================
# MEF2C KO/OE  (edges to both TBX21 and RUNX3 per network scan)
# ============================================================
TF3 = "MEF2C"
assert TF3 in oracle.active_regulatory_genes, \
    f"{TF3} is not an active regulator; check base GRN motif retention (script 02)."

genes = list(oracle.adata.var_names)
cyto_in = [g for g in CYTO if g in genes]
idx = [genes.index(g) for g in cyto_in]
base = oracle.adata.layers["imputed_count"]
tbx21_idx = genes.index("TBX21")
runx3_idx = genes.index("RUNX3")

def export_module_and_readouts(sim_layer, base_layer, tag):
    delta = np.asarray(sim_layer - base_layer)
    df_mod = pd.DataFrame(delta[:, idx], columns=cyto_in, index=oracle.adata.obs_names)
    df_mod.insert(0, "group",   oracle.adata.obs[GROUP].values)
    df_mod.insert(1, "patient", oracle.adata.obs[PATIENT].values)
    df_mod.to_csv(f"perturb_deltas_{tag}.csv")
    print(f"[{tag}] wrote perturb_deltas_{tag}.csv  cells={df_mod.shape[0]}")
    for readout_name, readout_idx in [("TBX21", tbx21_idx), ("RUNX3", runx3_idx)]:
        df_gene = pd.DataFrame({readout_name: delta[:, readout_idx]}, index=oracle.adata.obs_names)
        df_gene.insert(0, "group",   oracle.adata.obs[GROUP].values)
        df_gene.insert(1, "patient", oracle.adata.obs[PATIENT].values)
        df_gene.to_csv(f"perturb_deltas_{tag}_{readout_name}.csv")
        print(f"[{tag}->{readout_name}] wrote perturb_deltas_{tag}_{readout_name}.csv")

# ---- MEF2C KO ----
oracle.simulate_shift(perturb_condition={TF3: 0.0}, n_propagation=N_PROP)
sim_mef2c_ko = oracle.adata.layers["simulated_count"]
export_module_and_readouts(sim_mef2c_ko, base, "MEF2C_KO")

# ---- MEF2C OE ----
mef2c_hi = np.quantile(base[:, genes.index(TF3)], 0.9)
oracle.simulate_shift(perturb_condition={TF3: mef2c_hi}, n_propagation=N_PROP)
sim_mef2c_oe = oracle.adata.layers["simulated_count"]
export_module_and_readouts(sim_mef2c_oe, base, "MEF2C_OE")

# ============================================================
# KLF2 KO/OE  (edge to RUNX3 only per network scan)
# ============================================================
TF4 = "KLF2"
assert TF4 in oracle.active_regulatory_genes, \
    f"{TF4} is not an active regulator; check base GRN motif retention (script 02)."

def export_module_and_runx3_only(sim_layer, base_layer, tag):
    delta = np.asarray(sim_layer - base_layer)
    df_mod = pd.DataFrame(delta[:, idx], columns=cyto_in, index=oracle.adata.obs_names)
    df_mod.insert(0, "group",   oracle.adata.obs[GROUP].values)
    df_mod.insert(1, "patient", oracle.adata.obs[PATIENT].values)
    df_mod.to_csv(f"perturb_deltas_{tag}.csv")
    print(f"[{tag}] wrote perturb_deltas_{tag}.csv  cells={df_mod.shape[0]}")
    df_gene = pd.DataFrame({"RUNX3": delta[:, runx3_idx]}, index=oracle.adata.obs_names)
    df_gene.insert(0, "group",   oracle.adata.obs[GROUP].values)
    df_gene.insert(1, "patient", oracle.adata.obs[PATIENT].values)
    df_gene.to_csv(f"perturb_deltas_{tag}_RUNX3.csv")
    print(f"[{tag}->RUNX3] wrote perturb_deltas_{tag}_RUNX3.csv")

# ---- KLF2 KO ----
oracle.simulate_shift(perturb_condition={TF4: 0.0}, n_propagation=N_PROP)
sim_klf2_ko = oracle.adata.layers["simulated_count"]
export_module_and_runx3_only(sim_klf2_ko, base, "KLF2_KO")

# ---- KLF2 OE ----
klf2_hi = np.quantile(base[:, genes.index(TF4)], 0.9)
oracle.simulate_shift(perturb_condition={TF4: klf2_hi}, n_propagation=N_PROP)
sim_klf2_oe = oracle.adata.layers["simulated_count"]
export_module_and_runx3_only(sim_klf2_oe, base, "KLF2_OE")

##################################
### ATRXdel specific TFs only ####
##################################

ALREADY_TESTED = {"TBX21", "RUNX3", "MEF2C", "KLF2"}

# ---- rank TFs using ONLY the ATRXdel cluster's fitted network --------------
atrxdel_links = links.filtered_links["ATRXdel"]
cyto_edges = atrxdel_links[atrxdel_links["target"].isin(CYTO)]

ranked_atrxdel = (cyto_edges.groupby("source")
                    .agg(n_cyto_targets=("target", "nunique"),
                        total_coef_abs=("coef_abs", "sum"))
                    .sort_values("total_coef_abs", ascending=False))
ranked_atrxdel = ranked_atrxdel[~ranked_atrxdel.index.isin(ALREADY_TESTED)]

print("Top ATRXdel-specific cytotoxicity-module regulators (not yet tested):")
print(ranked_atrxdel.head(10))

top2_new_tfs = ranked_atrxdel.head(2).index.tolist()
print(f"\nSelected for simulation: {top2_new_tfs}")

# ---- simulate KO/OE for the top 2 new candidates ---------------------------
genes = list(oracle.adata.var_names)
cyto_in = [g for g in CYTO if g in genes]
idx = [genes.index(g) for g in cyto_in]
base = oracle.adata.layers["imputed_count"]

def export_module(sim_layer, base_layer, tag):
    delta = np.asarray(sim_layer - base_layer)
    df_mod = pd.DataFrame(delta[:, idx], columns=cyto_in, index=oracle.adata.obs_names)
    df_mod.insert(0, "group",   oracle.adata.obs[GROUP].values)
    df_mod.insert(1, "patient", oracle.adata.obs[PATIENT].values)
    df_mod.to_csv(f"perturb_deltas_{tag}.csv")
    print(f"[{tag}] wrote perturb_deltas_{tag}.csv  cells={df_mod.shape[0]}")

for tf in top2_new_tfs:
    assert tf in oracle.active_regulatory_genes, \
        f"{tf} is not an active regulator; check base GRN motif retention (script 02)."
    # KO
    oracle.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=N_PROP)
    sim_ko = oracle.adata.layers["simulated_count"]
    export_module(sim_ko, base, f"{tf}_KO")
    # OE
    tf_hi = np.quantile(base[:, genes.index(tf)], 0.9)
    oracle.simulate_shift(perturb_condition={tf: tf_hi}, n_propagation=N_PROP)
    sim_oe = oracle.adata.layers["simulated_count"]
    export_module(sim_oe, base, f"{tf}_OE")
    


# 
TF5 = "PRDM1"
runx3_idx = genes.index("RUNX3")  # reuse if already defined; otherwise: genes.index("RUNX3")

# ---- PRDM1 KO -> RUNX3 ----
oracle.simulate_shift(perturb_condition={TF5: 0.0}, n_propagation=N_PROP)
sim_prdm1_ko = oracle.adata.layers["simulated_count"]
delta_prdm1_ko_runx3 = np.asarray(sim_prdm1_ko - base)[:, runx3_idx]
df1 = pd.DataFrame({"RUNX3": delta_prdm1_ko_runx3}, index=oracle.adata.obs_names)
df1.insert(0, "group",   oracle.adata.obs[GROUP].values)
df1.insert(1, "patient", oracle.adata.obs[PATIENT].values)
df1.to_csv("perturb_deltas_PRDM1_KO_RUNX3.csv")
print(f"[PRDM1-KO->RUNX3] wrote perturb_deltas_PRDM1_KO_RUNX3.csv  cells={df1.shape[0]}")

# ---- PRDM1 OE -> RUNX3 ----
prdm1_hi = np.quantile(base[:, genes.index(TF5)], 0.9)
oracle.simulate_shift(perturb_condition={TF5: prdm1_hi}, n_propagation=N_PROP)
sim_prdm1_oe = oracle.adata.layers["simulated_count"]
delta_prdm1_oe_runx3 = np.asarray(sim_prdm1_oe - base)[:, runx3_idx]
df2 = pd.DataFrame({"RUNX3": delta_prdm1_oe_runx3}, index=oracle.adata.obs_names)
df2.insert(0, "group",   oracle.adata.obs[GROUP].values)
df2.insert(1, "patient", oracle.adata.obs[PATIENT].values)
df2.to_csv("perturb_deltas_PRDM1_OE_RUNX3.csv")
print(f"[PRDM1-OE->RUNX3] wrote perturb_deltas_PRDM1_OE_RUNX3.csv  cells={df2.shape[0]}")


# CD6 signalling pathway TFs (AP-1, NF-kB, NFAT) into TBX21/RUNX3/PRDM1/cytotoxicity genes

import celloracle as co
import pandas as pd

links = co.load_hdf5("links_per_group.celloracle.links")

CD6_PATHWAY_TFS = [
    # AP-1 family (Ras/MAPK/ERK arm)
    "FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND",
    # NF-kB family (PKCtheta arm)
    "RELA", "RELB", "REL", "NFKB1", "NFKB2",
    # NFAT family (calcium influx arm)
    "NFATC1", "NFATC2", "NFATC3", "NFATC4",
]
TARGETS = ["TBX21", "RUNX3", "PRDM1"] + \
    ["GZMB","PRF1","IFNG","GNLY","NKG7","KLRD1","GZMA","GZMH","KLRF1","FCGR3A"]

rows = []
for group_name in links.filtered_links.keys():
    df = links.filtered_links[group_name]
    hits = df[(df["source"].isin(CD6_PATHWAY_TFS)) & (df["target"].isin(TARGETS))]
    for _, row in hits.iterrows():
        rows.append({
            "group": group_name, "source": row["source"], "target": row["target"],
            "coef_mean": row["coef_mean"], "coef_abs": row["coef_abs"], "p": row["p"]
        })

cd6_edge_table = pd.DataFrame(rows).sort_values(["target","source","group"])
cd6_edge_table.to_csv("cd6_pathway_tf_edges.csv", index=False)
print(f"Found {cd6_edge_table.shape[0]} edges from CD6-pathway TFs into TBX21/RUNX3/PRDM1/cytotoxicity genes")
print(cd6_edge_table.to_string(index=False))

# ---- flag which are ATRXdel-specific (present in ATRXdel, absent elsewhere) or
#      ATRXdel-absent (present elsewhere, absent specifically in ATRXdel) ------
pivot = cd6_edge_table.pivot_table(index=["source","target"], columns="group",
                                     values="coef_mean", aggfunc="first")
print("\n--- Presence/absence by group (NaN = no edge) ---")
print(pivot)


# TBX21 repression by FOS TF

TF6 = "FOS"
assert TF6 in oracle.active_regulatory_genes, \
    f"{TF6} is not an active regulator; check base GRN motif retention (script 02)."

genes = list(oracle.adata.var_names)
cyto_in = [g for g in CYTO if g in genes]
idx = [genes.index(g) for g in cyto_in]
base = oracle.adata.layers["imputed_count"]
tbx21_idx = genes.index("TBX21")

def export_module_and_tbx21(sim_layer, base_layer, tag):
    delta = np.asarray(sim_layer - base_layer)
    df_mod = pd.DataFrame(delta[:, idx], columns=cyto_in, index=oracle.adata.obs_names)
    df_mod.insert(0, "group",   oracle.adata.obs[GROUP].values)
    df_mod.insert(1, "patient", oracle.adata.obs[PATIENT].values)
    df_mod.to_csv(f"perturb_deltas_{tag}.csv")
    print(f"[{tag}] wrote perturb_deltas_{tag}.csv  cells={df_mod.shape[0]}")
    df_gene = pd.DataFrame({"TBX21": delta[:, tbx21_idx]}, index=oracle.adata.obs_names)
    df_gene.insert(0, "group",   oracle.adata.obs[GROUP].values)
    df_gene.insert(1, "patient", oracle.adata.obs[PATIENT].values)
    df_gene.to_csv(f"perturb_deltas_{tag}_TBX21.csv")
    print(f"[{tag}->TBX21] wrote perturb_deltas_{tag}_TBX21.csv")

# ---- FOS KO ----
oracle.simulate_shift(perturb_condition={TF6: 0.0}, n_propagation=N_PROP)
sim_fos_ko = oracle.adata.layers["simulated_count"]
export_module_and_tbx21(sim_fos_ko, base, "FOS_KO")

# ---- FOS OE ----
fos_hi = np.quantile(base[:, genes.index(TF6)], 0.9)
oracle.simulate_shift(perturb_condition={TF6: fos_hi}, n_propagation=N_PROP)
sim_fos_oe = oracle.adata.layers["simulated_count"]
export_module_and_tbx21(sim_fos_oe, base, "FOS_OE")