"""
10_additional_TF_perturbations.py
KO/OE simulations for three TFs flagged from 08_top_TF_perturbations.py's
output as showing an ATRXdel-vs-control-specific signal into RUNX3 (not just
"top ranked overall" -- these were picked because their edge coefficient is
distinctly different in the ATRXdel-labeled units vs every other group,
within the same NK subtype cluster):

  * JUND  -> RUNX3: highest coefficient of all 4 genetic groups in BOTH
    ATRXdel-bearing clusters (hNK_Bm2: 0.197 vs control 0.157; hNK_Bm3: 0.189
    vs control 0.134) -- a reproducible ATRXdel > control gain, not just
    ATRXdel presence.
  * KLF6 -> RUNX3: repressive edge (-0.047 hNK_Bm2_ATRXdel, -0.044
    hNK_Bm3_ATRXdel) that does not survive the p<0.001/top-2000 filter in
    ANY other group in either cluster -- an ATRXdel-gained repressive edge.
  * STAT1 -> RUNX3: weaker case (also appears in hNK_Bm3_MYCNamp, so not as
    cleanly ATRXdel-restricted as the two above) but STAT1 is the
    highest-ranked untested TF overall by total cytotoxicity-edge strength
    (rank 5, 33 targets) -- included as a general-interest third candidate.
  * ZNF281 -> RUNX3: mirror image of JUND -- ATRXdel has the WEAKEST
    ZNF281->RUNX3 activation of all 4 groups in BOTH clusters (hNK_Bm2:
    ATRXdel 0.015 vs control 0.053, ATRXwtMYCNwt 0.085, MYCNamp 0.055;
    hNK_Bm3: ATRXdel 0.026 vs control 0.049, ATRXwtMYCNwt 0.090, MYCNamp
    0.052) -- also the single strongest untested TF overall by
    total_coef_abs (rank 6, 5.36).
  * TAL1: broadest untested direct cytotoxicity-gene connectivity (53
    targets -- more than KLF2's 25, close to POU2F2's 47) but never appears
    linked to TBX21/RUNX3 in any unit, so its cytotoxicity-gene effect (if
    any) runs through a route none of the other tested TFs touch -- module
    score only, no RUNX3/TBX21 readout to chase.

Standalone script -- reloads `oracle`/`links` the same way as
08_top_TF_perturbations.py, and shares its export/perturb helpers via
grn_perturb_utils.py rather than redefining them.

Run any time after 03_grn_perturb.py.
"""
import celloracle as co

from grn_perturb_utils import make_perturb_helpers

ORACLE_PATH = "oracle_fitted.celloracle.oracle"
GENETIC_GROUP = "genetic_group"
CLUSTER       = "predicted_NK_type_de.n350"
PATIENT       = "patient_id"
N_PROP  = 3

CYTO = ["GZMB", "PRF1", "IFNG", "GNLY", "NKG7", "KLRD1", "GZMA", "GZMH", "KLRF1", "FCGR3A"]

oracle = co.load_hdf5(ORACLE_PATH)

genes   = list(oracle.adata.var_names)
cyto_in = [g for g in CYTO if g in genes]
idx     = [genes.index(g) for g in cyto_in]
base    = oracle.adata.layers["imputed_count"]

export_perturbation, perturb_ko_oe = make_perturb_helpers(
    oracle, genes, cyto_in, idx, base, GENETIC_GROUP, CLUSTER, PATIENT, N_PROP
)

# =============================================================================
# JUND KO/OE -> effect on RUNX3 (direct edge) and TBX21 (downstream)
# =============================================================================
perturb_ko_oe("JUND", "JUND", extra_readouts=("RUNX3", "TBX21"))

# =============================================================================
# KLF6 KO/OE -> effect on RUNX3
# KO is the direction of interest: if the ATRXdel-gained repressive edge is
# real, KO should de-repress RUNX3 specifically in ATRXdel cells with little
# effect in control (same asymmetric-baseline logic as the TBX21 KO design
# in 03_grn_perturb.py).
# =============================================================================
perturb_ko_oe("KLF6", "KLF6", extra_readouts=("RUNX3",))

# =============================================================================
# STAT1 KO/OE -> effect on RUNX3 and TBX21
# =============================================================================
perturb_ko_oe("STAT1", "STAT1", extra_readouts=("RUNX3", "TBX21"))

# =============================================================================
# ZNF281 KO/OE -> effect on RUNX3
# OE is the direction of interest here: ATRXdel shows the weakest
# ZNF281->RUNX3 activation of all groups, so OE (pushing ZNF281 to its 90th
# percentile) is the more informative "rescue-style" test -- does restoring
# strong ZNF281 activity push ATRXdel's RUNX3 signal back toward control?
# =============================================================================
perturb_ko_oe("ZNF281", "ZNF281", extra_readouts=("RUNX3",))

# =============================================================================
# TAL1 KO/OE -> cytotoxicity module only (no RUNX3/TBX21 edge in any unit)
# =============================================================================
perturb_ko_oe("TAL1", "TAL1")

print("\n[additional TF perturbations] done")
