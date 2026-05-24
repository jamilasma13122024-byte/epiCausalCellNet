
# epiCausalCellNet: Granger Causality-based
# Cell-Cell Communication Network Inference
# Author: [Your Name]
# Date: 2024

# ── Dependencies ──
# pip install scanpy statsmodels gseapy scipy pandas numpy matplotlib

import scanpy as sc
import pandas as pd
import numpy as np
from statsmodels.tsa.stattools import grangercausalitytests
import warnings
warnings.filterwarnings("ignore")

# ══════════════════════════════════════════
# STEP 1: Load Data
# ══════════════════════════════════════════
def load_data(h5_path, meta_path):
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    meta = pd.read_csv(meta_path)
    meta.index = meta["Barcode"]
    adata.obs = adata.obs.join(meta, how="left")
    return adata

# ══════════════════════════════════════════
# STEP 2: Quality Control
# ══════════════════════════════════════════
def quality_control(adata, min_genes=200, max_genes=6000, max_mt=20):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"],
        percent_top=None, log1p=False, inplace=True)
    before = adata.n_obs
    adata = adata[
        (adata.obs["n_genes_by_counts"] > min_genes) &
        (adata.obs["n_genes_by_counts"] < max_genes) &
        (adata.obs["pct_counts_mt"] < max_mt)
    ].copy()
    print(f"QC: {before} -> {adata.n_obs} cells retained")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata

# ══════════════════════════════════════════
# STEP 3: Granger Causality
# ══════════════════════════════════════════
def granger_causality(adata, lr_pairs, ct_pairs, 
                      diagnosis_col="Diagnosis",
                      celltype_col="Cell.Type",
                      sample_col="SampleID",
                      diagnosis_val="AD",
                      maxlag=2):
    mask = adata.obs[diagnosis_col] == diagnosis_val
    adata_d = adata[mask]
    samples = adata_d.obs[sample_col].unique()
    results = []
    
    for ligand, receptor in lr_pairs:
        if ligand not in adata.var_names: continue
        if receptor not in adata.var_names: continue
        
        for sender, receiver in ct_pairs:
            lig_series, rec_series = [], []
            
            for samp in samples:
                mask_samp = adata_d.obs[sample_col] == samp
                mask_s = (adata_d.obs[celltype_col]==sender) & mask_samp
                mask_r = (adata_d.obs[celltype_col]==receiver) & mask_samp
                
                if mask_s.sum() > 0 and mask_r.sum() > 0:
                    lig_series.append(float(
                        adata_d[mask_s, ligand].X.mean()))
                    rec_series.append(float(
                        adata_d[mask_r, receptor].X.mean()))
            
            min_len = min(len(lig_series), len(rec_series))
            if min_len < 5: continue
            
            try:
                df_gc = pd.DataFrame({
                    "Ligand":   lig_series[:min_len],
                    "Receptor": rec_series[:min_len]
                })
                test = grangercausalitytests(
                    df_gc[["Receptor","Ligand"]],
                    maxlag=maxlag, verbose=False
                )
                p_val = test[1][0]["ssr_ftest"][1]
                results.append({
                    "Diagnosis": diagnosis_val,
                    "Sender":    sender,
                    "Receiver":  receiver,
                    "Ligand":    ligand,
                    "Receptor":  receptor,
                    "P_value":   round(p_val, 4),
                    "Causal":    "YES" if p_val < 0.05 else "NO"
                })
            except: continue
    
    return pd.DataFrame(results)

# ══════════════════════════════════════════
# STEP 4: Permutation Validation
# ══════════════════════════════════════════
def permutation_granger(lig_vals, rec_vals, 
                        n_perm=1000, maxlag=2):
    df_real = pd.DataFrame({
        "Ligand": lig_vals, "Receptor": rec_vals})
    try:
        real_test = grangercausalitytests(
            df_real[["Receptor","Ligand"]],
            maxlag=maxlag, verbose=False)
        real_f = real_test[1][0]["ssr_ftest"][0]
        real_p = real_test[1][0]["ssr_ftest"][1]
    except: return None, None
    
    perm_f_stats = []
    for _ in range(n_perm):
        lig_perm = np.random.permutation(lig_vals)
        df_perm = pd.DataFrame({
            "Ligand": lig_perm, "Receptor": rec_vals})
        try:
            perm_test = grangercausalitytests(
                df_perm[["Receptor","Ligand"]],
                maxlag=maxlag, verbose=False)
            perm_f_stats.append(
                perm_test[1][0]["ssr_ftest"][0])
        except: continue
    
    if not perm_f_stats: return real_p, None
    perm_p = np.mean(np.array(perm_f_stats) >= real_f)
    return real_p, perm_p

# ══════════════════════════════════════════
# STEP 5: Differential Expression
# ══════════════════════════════════════════
def differential_expression(adata, genes,
                             diagnosis_col="Diagnosis"):
    from scipy import stats
    mask_ad   = adata.obs[diagnosis_col] == "AD"
    mask_ctrl = adata.obs[diagnosis_col] == "Control"
    results = []
    
    for gene in genes:
        if gene not in adata.var_names: continue
        ad_expr   = adata[mask_ad,   gene].X.toarray().flatten()
        ctrl_expr = adata[mask_ctrl, gene].X.toarray().flatten()
        log2fc = np.log2((ad_expr.mean()+1e-9) / 
                         (ctrl_expr.mean()+1e-9))
        _, pval = stats.mannwhitneyu(ad_expr, ctrl_expr,
                                     alternative="two-sided")
        results.append({
            "Gene": gene,
            "AD_mean":   round(ad_expr.mean(), 4),
            "Ctrl_mean": round(ctrl_expr.mean(), 4),
            "Log2FC":    round(log2fc, 4),
            "P_value":   round(pval, 6),
            "Regulation": "UP" if log2fc > 0 else "DOWN"
        })
    return pd.DataFrame(results)

# ══════════════════════════════════════════
# MAIN — Example Usage
# ══════════════════════════════════════════
if __name__ == "__main__":
    # L-R pairs from literature
    LR_PAIRS = [
        ("CCL2",   "CCR2"),
        ("APOE",   "LRP1"),
        ("CX3CL1", "CX3CR1"),
        ("NRG1",   "ERBB4"),
        ("TLR1",   "MYD88"),
        ("APP",    "CD74"),
        ("PROS1",  "AXL"),
        ("PTPRC",  "MRC1"),
    ]
    CT_PAIRS = [
        ("ASC", "MG"),
        ("ODC", "MG"),
        ("MG",  "EX"),
        ("ODC", "ASC"),
        ("MG",  "INH"),
    ]
    
    print("epiCausalCellNet — Granger Causality Framework")
    print("Usage:")
    print("  1. adata = load_data(h5_path, meta_path)")
    print("  2. adata = quality_control(adata)")
    print("  3. results = granger_causality(adata, LR_PAIRS, CT_PAIRS)")
    print("  4. Validate: permutation_granger(lig_vals, rec_vals)")
    print("  5. de = differential_expression(adata, genes)")
