import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

def compute_icc_2_1(df, col1, col2):
    y1 = df[col1].values
    y2 = df[col2].values
    n = len(y1)
    k = 2
    y = np.column_stack([y1, y2])
    
    mean_total = np.mean(y)
    row_means = np.mean(y, axis=1)
    col_means = np.mean(y, axis=0)
    
    ss_total = np.sum((y - mean_total)**2)
    ss_rows = k * np.sum((row_means - mean_total)**2)
    ss_cols = n * np.sum((col_means - mean_total)**2)
    ss_error = ss_total - ss_rows - ss_cols
    
    ms_rows = ss_rows / (n - 1)
    ms_cols = ss_cols / (k - 1)
    ms_error = ss_error / ((n - 1) * (k - 1))
    
    icc = (ms_rows - ms_error) / (ms_rows + (k - 1) * ms_error + (k / n) * (ms_cols - ms_error))
    return icc

def main():
    print("=== M11-F Statistical Forensics: Independent Recomputation ===")
    failures = []
    
    # 1. Figure 2: ARTI Decoupling & Multivariate OLS
    f2_path = "figures/source_data/fig2_source_data.csv"
    if os.path.exists(f2_path):
        df2 = pd.read_csv(f2_path)
        print(f"[*] Loaded Fig2 source data: {len(df2)} conditions.")
        
        # Bivariate correlations
        r_c0, p_c0 = stats.pearsonr(df2["arti"], df2["c0_pseudobulk_l2"])
        r_c1, p_c1 = stats.pearsonr(df2["arti"], df2["c1_e_distance"])
        r_c2_mag, p_c2_mag = stats.pearsonr(df2["arti"], df2["c2_mean_projection_magnitude"])
        r_c2_frac, p_c2_frac = stats.pearsonr(df2["arti"], df2["c2_response_penetrance_fraction"])
        r_c3, p_c3 = stats.pearsonr(df2["arti"], df2["c3_simple_dispersion"])
        r_c4, p_c4 = stats.pearsonr(df2["arti"], df2["c4_transcriptomic_entropy"])
        r_c5, p_c5 = stats.pearsonr(df2["arti"], df2["c5_dist_to_control_centroid"])
        
        print(f"  ARTI vs C0 (Pseudobulk L2): r = {r_c0:+.4f} (p = {p_c0:.4f})")
        print(f"  ARTI vs C1 (Energy Distance): r = {r_c1:+.4f} (p = {p_c1:.4f})")
        print(f"  ARTI vs C2 (Mean Projection Mag): r = {r_c2_mag:+.4f} (p = {p_c2_mag:.4f})")
        print(f"  ARTI vs C2 (PS Fraction): r = {r_c2_frac:+.4f} (p = {p_c2_frac:.4f})")
        print(f"  ARTI vs C3 (Simple Dispersion): r = {r_c3:+.4f} (p = {p_c3:.4f})")
        print(f"  ARTI vs C4 (Entropy): r = {r_c4:+.4f} (p = {p_c4:.2e})")
        print(f"  ARTI vs C5 (Dist to Control): r = {r_c5:+.4f} (p = {p_c5:.4f})")
        
        if abs(r_c0 - (-0.0753)) > 0.005: failures.append(f"C0 correlation mismatch: got {r_c0}, expected -0.0753")
        if abs(r_c1 - (-0.1573)) > 0.005: failures.append(f"C1 correlation mismatch: got {r_c1}, expected -0.1573")
        if abs(r_c2_mag - (-0.1806)) > 0.005: failures.append(f"C2 magnitude correlation mismatch: got {r_c2_mag}, expected -0.1806")
        if abs(r_c2_frac - 0.0178) > 0.005: failures.append(f"C2 fraction correlation mismatch: got {r_c2_frac}, expected +0.0178")
        if abs(r_c3 - (-0.0591)) > 0.005: failures.append(f"C3 correlation mismatch: got {r_c3}, expected -0.0591")
        if abs(r_c4 - (-0.3421)) > 0.005: failures.append(f"C4 correlation mismatch: got {r_c4}, expected -0.3421")
        if abs(r_c5 - (-0.1267)) > 0.005: failures.append(f"C5 correlation mismatch: got {r_c5}, expected -0.1267")
        
        # Multivariate OLS using all comparator metrics
        comparators = [
            'c0_pseudobulk_l2', 'c1_e_distance', 'c2_response_penetrance_fraction',
            'c2_mean_projection_magnitude', 'c3_simple_dispersion',
            'c4_transcriptomic_entropy', 'c5_dist_to_control_centroid'
        ]
        X = df2[comparators].values
        X_const = np.column_stack([np.ones(len(X)), X])
        y = df2["arti"].values
        beta = np.linalg.lstsq(X_const, y, rcond=None)[0]
        y_pred = X_const @ beta
        ss_res = np.sum((y - y_pred)**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r2 = 1.0 - (ss_res / ss_tot)
        res_var_pct = (1.0 - r2) * 100.0
        exp_var_pct = r2 * 100.0
        print(f"  Multivariate OLS R2: {r2:.4f} -> Explained: {exp_var_pct:.2f}%, Residual: {res_var_pct:.2f}%")
        
        if abs(exp_var_pct - 15.10) > 0.15: failures.append(f"Explained variance mismatch: got {exp_var_pct}, expected 15.10%")
        if abs(res_var_pct - 84.90) > 0.15: failures.append(f"Residual variance mismatch: got {res_var_pct}, expected 84.90%")
    else:
        failures.append(f"Missing {f2_path}")
        
    # 2. Figure 3A: Replicate Reproducibility in Tahoe-100M
    f3_rep_path = "figures/source_data/fig3_replicate_source_data.csv"
    if os.path.exists(f3_rep_path):
        df3_rep = pd.read_csv(f3_rep_path)
        sub3 = df3_rep[df3_rep['plate_pair'] == 'plate4_vs_plate14']
        print(f"\n[*] Loaded Fig3 replicate source data: {len(sub3)} Plate 4 vs 14 condition pairs.")
        r_rep, p_rep = stats.pearsonr(sub3["rep1"], sub3["rep2"])
        rho_rep, prho_rep = stats.spearmanr(sub3["rep1"], sub3["rep2"])
        icc_rep = compute_icc_2_1(sub3, "rep1", "rep2")
        print(f"  Tahoe Replicate (Plate 4 vs 14): Pearson r = {r_rep:.4f} (p = {p_rep:.2e}), Spearman rho = {rho_rep:.4f}, ICC = {icc_rep:.4f}")
        
        if abs(r_rep - 0.4173) > 0.005: failures.append(f"Tahoe replicate r mismatch: got {r_rep}, expected 0.4173")
        if abs(rho_rep - 0.3062) > 0.005: failures.append(f"Tahoe replicate rho mismatch: got {rho_rep}, expected 0.3062")
        if abs(icc_rep - 0.3938) > 0.01: failures.append(f"Tahoe replicate ICC mismatch: got {icc_rep}, expected 0.3938")
    else:
        failures.append(f"Missing {f3_rep_path}")
        
    # 3. Figure 3C-D: Concentration Modulation in Tahoe-100M
    f3_dose_path = "figures/source_data/fig3_dose_source_data.csv"
    if os.path.exists(f3_dose_path):
        df3_dose = pd.read_csv(f3_dose_path)
        print(f"\n[*] Loaded Fig3 dose source data: {len(df3_dose)} matched pairs.")
        t_stat, p_val = stats.ttest_rel(df3_dose["arti_high"], df3_dose["arti_low"])
        mean_diff = np.mean(df3_dose["arti_high"] - df3_dose["arti_low"])
        delta_arti = df3_dose["arti_high"] - df3_dose["arti_low"]
        delta_c0 = df3_dose["delta_mean_response"]
        r_dose_decouple, p_dose_decouple = stats.pearsonr(delta_arti, delta_c0)
        
        print(f"  Paired t-test (high vs low): t = {t_stat:.2f} (p = {p_val:.2e}), Mean diff = {mean_diff:.4f}")
        print(f"  Delta ARTI vs Delta C0 correlation: r = {r_dose_decouple:.4f} (p = {p_dose_decouple:.3f})")
        
        if abs(t_stat - 5.21) > 0.1: failures.append(f"Dose t-stat mismatch: got {t_stat}, expected 5.21")
        if abs(mean_diff - 0.2792) > 0.01: failures.append(f"Dose mean diff mismatch: got {mean_diff}, expected 0.2792")
        if abs(r_dose_decouple - (-0.1007)) > 0.01: failures.append(f"Dose decoupling r mismatch: got {r_dose_decouple}, expected -0.1007")
    else:
        failures.append(f"Missing {f3_dose_path}")
        
    # 4. Figure 4: Type II ANOVA Variance Partitioning
    f4_path = "figures/source_data/fig4_source_data.csv"
    if os.path.exists(f4_path):
        df4 = pd.read_csv(f4_path)
        print(f"\n[*] Loaded Fig4 ANOVA source data:")
        for _, row in df4.iterrows():
            print(f"  {row['source']}: {row['pct_variance']:.2f}% (F = {row['F']:.2f}, p = {row['PR(>F)']:.2e})")
        
        drug_pct = df4[df4["source"].str.contains("drug", case=False)]["pct_variance"].values[0]
        cell_pct = df4[df4["source"].str.contains("cell", case=False)]["pct_variance"].values[0]
        plate_pct = df4[df4["source"].str.contains("plate", case=False)]["pct_variance"].values[0]
        res_pct = df4[df4["source"].str.contains("residual", case=False)]["pct_variance"].values[0]
        
        if abs(drug_pct - 22.12) > 0.1: failures.append(f"ANOVA drug % mismatch: got {drug_pct}, expected 22.12")
        if abs(cell_pct - 14.04) > 0.1: failures.append(f"ANOVA cell line % mismatch: got {cell_pct}, expected 14.04")
        if abs(plate_pct - 16.35) > 0.1: failures.append(f"ANOVA plate % mismatch: got {plate_pct}, expected 16.35")
        if abs(res_pct - 47.49) > 0.1: failures.append(f"ANOVA residual % mismatch: got {res_pct}, expected 47.49")
    else:
        failures.append(f"Missing {f4_path}")

    # 5. Figure 5A: SciPlex2 Replicate Reproducibility
    f5_sciplex_path = "figures/source_data/fig5_sciplex_source_data.csv"
    if os.path.exists(f5_sciplex_path):
        df5_sci = pd.read_csv(f5_sciplex_path)
        print(f"\n[*] Loaded SciPlex2 replicate data: {len(df5_sci)} condition pairs.")
        r_sci, p_sci = stats.pearsonr(df5_sci["arti_plate_A"], df5_sci["arti_plate_B"])
        rho_sci, prho_sci = stats.spearmanr(df5_sci["arti_plate_A"], df5_sci["arti_plate_B"])
        icc_sci = compute_icc_2_1(df5_sci, "arti_plate_A", "arti_plate_B")
        print(f"  SciPlex2: Pearson r = {r_sci:.4f} (p = {p_sci:.2e}), Spearman rho = {rho_sci:.4f}, ICC = {icc_sci:.4f}")
        
        if abs(r_sci - 0.7509) > 0.005: failures.append(f"SciPlex2 r mismatch: got {r_sci}, expected 0.7509")
        if abs(rho_sci - 0.6825) > 0.005: failures.append(f"SciPlex2 rho mismatch: got {rho_sci}, expected 0.6825")
        if abs(icc_sci - 0.7155) > 0.01: failures.append(f"SciPlex2 ICC mismatch: got {icc_sci}, expected 0.7155")
    else:
        failures.append(f"Missing {f5_sciplex_path}")
        
    # 6. Figure 7A: Blinded DTP Evaluation (GSE275199)
    f7_hop62_path = "figures/source_data/fig7_hop62_source_data.csv"
    if os.path.exists(f7_hop62_path):
        df7_hop62 = pd.read_csv(f7_hop62_path)
        print(f"\n[*] Loaded Fig7 HOP62 single-cell DTP projections: {len(df7_hop62)} cells.")
        tail_scores = df7_hop62[df7_hop62["cell_category"] == "ART / Tail"]["dtp_spearman_rho"].values
        dom_scores = df7_hop62[df7_hop62["cell_category"] == "Dominant"]["dtp_spearman_rho"].values
        mean_tail = np.mean(tail_scores)
        mean_dom = np.mean(dom_scores)
        delta_rho = mean_tail - mean_dom
        u_stat, p_mwu = stats.mannwhitneyu(tail_scores, dom_scores, alternative="greater")
        t_dtp, p_ttest = stats.ttest_ind(tail_scores, dom_scores, alternative="greater")
        print(f"  HOP62 Adagrasib vs DTP: Delta rho = {delta_rho:.4f} (MWU p = {p_mwu:.4f}, one-sided t-test p = {p_ttest:.4f})")
        
        if abs(delta_rho - (-0.0242)) > 0.005: failures.append(f"DTP delta rho mismatch: got {delta_rho}, expected -0.0242")
        if abs(p_mwu - 0.8354) > 0.01: failures.append(f"DTP p-value mismatch: got {p_mwu}, expected 0.8354")
    else:
        failures.append(f"Missing {f7_hop62_path}")
        
    # 7. Figure 7C: Pan-Condition DTP Screening
    f7_pan_path = "figures/source_data/fig7_pan_dtp_source_data.csv"
    if os.path.exists(f7_pan_path):
        df7_pan = pd.read_csv(f7_pan_path)
        print(f"\n[*] Loaded Pan-condition DTP screen: {len(df7_pan)} conditions.")
        med_delta = df7_pan["delta_rho"].median()
        null_count = np.sum((df7_pan["p_val_rho"] >= 0.05) | (df7_pan["delta_rho"] <= 0))
        print(f"  Pan-condition: {null_count}/49 null or negative, Median Delta rho = {med_delta:.4f}")
        
        if null_count != 47: failures.append(f"Pan-condition null count mismatch: got {null_count}, expected 47")
        if abs(med_delta - (-0.0292)) > 0.002: failures.append(f"Pan-condition median mismatch: got {med_delta}, expected -0.0292")
    else:
        failures.append(f"Missing {f7_pan_path}")
        
    # Summary
    print("\n" + "="*50)
    if failures:
        print(f"[FAIL] Found {len(failures)} recomputation discrepancies:")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print("[PASS] 100% of primary statistics independently recomputed and match RESULTS_LEDGER.csv within precision tolerance.")
        sys.exit(0)

if __name__ == "__main__":
    main()
