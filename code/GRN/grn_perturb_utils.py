"""
grn_perturb_utils.py
Shared CSV-export / KO-OE helpers for the TF perturbation scripts
(08_top_TF_perturbations.py, 10_additional_TF_perturbations.py). Both reload
the same fitted Oracle object and need identical export logic; factored out
here so it isn't duplicated verbatim across scripts (same rationale as
pseudobulk_helpers.R on the R side).
"""
import numpy as np
import pandas as pd


def make_perturb_helpers(oracle, genes, cyto_in, idx, base,
                          GENETIC_GROUP, CLUSTER, PATIENT, N_PROP):
    """Returns (export_perturbation, perturb_ko_oe) closures bound to the
    given oracle/context, so callers don't need to thread every constant
    through every call."""

    def export_perturbation(sim_layer, tag, extra_readouts=()):
        """Write the cytotoxicity-module delta CSV, plus one CSV per extra
        single-gene readout (e.g. TBX21, RUNX3), each with
        group/cluster/patient columns."""
        delta = np.asarray(sim_layer - base)

        df_mod = pd.DataFrame(delta[:, idx], columns=cyto_in, index=oracle.adata.obs_names)
        df_mod.insert(0, "group",   oracle.adata.obs[GENETIC_GROUP].values)
        df_mod.insert(1, "cluster", oracle.adata.obs[CLUSTER].values)
        df_mod.insert(2, "patient", oracle.adata.obs[PATIENT].values)
        df_mod.to_csv(f"perturb_deltas_{tag}.csv")
        print(f"[{tag}] wrote perturb_deltas_{tag}.csv  cells={df_mod.shape[0]}")

        for readout_name in extra_readouts:
            readout_idx = genes.index(readout_name)
            df_gene = pd.DataFrame({readout_name: delta[:, readout_idx]}, index=oracle.adata.obs_names)
            df_gene.insert(0, "group",   oracle.adata.obs[GENETIC_GROUP].values)
            df_gene.insert(1, "cluster", oracle.adata.obs[CLUSTER].values)
            df_gene.insert(2, "patient", oracle.adata.obs[PATIENT].values)
            df_gene.to_csv(f"perturb_deltas_{tag}_{readout_name}.csv")
            print(f"[{tag}->{readout_name}] wrote perturb_deltas_{tag}_{readout_name}.csv")

    def perturb_ko_oe(tf, tag_prefix, extra_readouts=()):
        """Run KO (set to 0) then OE (90th-percentile baseline) for one TF,
        exporting both through export_perturbation()."""
        assert tf in oracle.active_regulatory_genes, \
            f"{tf} is not an active regulator; check base GRN motif retention (script 02)."

        oracle.simulate_shift(perturb_condition={tf: 0.0}, n_propagation=N_PROP)
        export_perturbation(oracle.adata.layers["simulated_count"], f"{tag_prefix}_KO", extra_readouts)

        tf_hi = np.quantile(base[:, genes.index(tf)], 0.9)
        oracle.simulate_shift(perturb_condition={tf: tf_hi}, n_propagation=N_PROP)
        export_perturbation(oracle.adata.layers["simulated_count"], f"{tag_prefix}_OE", extra_readouts)

    return export_perturbation, perturb_ko_oe
