"""
02_base_grn.py
Build the CellOracle BASE GRN (candidate TF -> target edges) from snATAC-seq data
peaks + Cicero co-accessibility. One shared backbone for all groups, so that the
per-group differences later reflect expression-driven edge activation, not
different peak sets.

Refs:
  CellOracle ......... Kamimoto et al., Nature 2023 (PMID 36755098)
  Motif scan ......... gimmemotifs (van Heeringen lab), default vertebrate set

Prior env requirements:
  pip install celloracle
  python -c "import genomepy; genomepy.install_genome('hg38','UCSC')"
"""
import pandas as pd
import celloracle as co
from celloracle import motif_analysis as ma

IN_DIR          = "/home/sara_wz/bioinf_isilon/Research/HALBRITTER/zHalbritter_TaschnerMandl/wernig-zorc.sara/out/NK_project/h5ad"
OUT             = "/home/sara_wz/bioinf_isilon/Research/HALBRITTER/zHalbritter_TaschnerMandl/wernig-zorc.sara/out/NK_project/GRN/base_GRN_hg38.parquet"
REF             = "hg38"
COACCESS        = 0.8     # keep strong cis-connections only
FPR             = 0.02    # motif-scan false-positive rate (CellOracle default)
MOTIF_SCORE_THR = 10      # motif score cutoff (CellOracle default)

# candidate accessible peaks (chr_start_end) + cicero connections (Peak1,Peak2,coaccess)
peaks  = pd.read_csv(f"{IN_DIR}/all_peaks.csv")["peak"].tolist()
cicero = pd.read_csv(f"{IN_DIR}/cicero_connections.csv")

# 1) annotate TSS, integrate with co-accessibility -> cis-regulatory peaks per gene
tss = ma.get_tss_info(peak_str_list=peaks, ref_genome=REF)
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss, cicero_connections=cicero)
peak_df = (integrated[integrated.coaccess >= COACCESS]
          [["peak_id", "gene_short_name"]]
          .reset_index(drop=True))
print(f"[base GRN] cis-regulatory peaks retained: {peak_df.shape[0]}")

# 2) scan those peaks for TF motifs -> candidate TF->target edges
tfi = ma.TFinfo(peak_data_frame=peak_df, ref_genome=REF)
tfi.scan(fpr=FPR, motifs=None, verbose=True)      # None = default motif set
tfi.reset_filtering()
tfi.filter_motifs_by_score(threshold=MOTIF_SCORE_THR)
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

base_GRN = tfi.to_dataframe()
base_GRN.to_parquet(OUT)
print(f"[base GRN] saved {OUT}  shape={base_GRN.shape}")

# sanity: confirm TBX21 and RUNX3 survive as regulators in the backbone
for tf in ("TBX21", "RUNX3"):
    present = tf in base_GRN.columns
    print(f"[base GRN] {tf} present as regulator: {present}")
    # If TBX21 is absent, its motif wasn't found under any retained peak -> the
    # perturbation can't act on it -> Loosen MOTIF_SCORE_THR or COACCESS if so.
