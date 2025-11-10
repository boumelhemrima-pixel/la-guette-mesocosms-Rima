# -*- coding: utf-8 -*-
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import t
import matplotlib.dates as mdates
from datetime import datetime
from matplotlib.lines import Line2D

# ======================================
# USER SETTINGS — EDIT THESE TWO PATHS
# ======================================
data_file = r"C:/Users/rboumelh/Desktop/Mesocosms/Modeling GPP/GPPmax_by campaign/GPPvsPAR_modified.xlsx"
output_directory = r"C:/Users/rboumelh/Desktop/Mesocosms/Modeling GPP/GPPmax_by campaign/by replicate/GPPmax/bounded/for git"
os.makedirs(output_directory, exist_ok=True)

# ======================================
# Model & helpers
# ======================================
def gpp_model(PAR, GPP_max, k):
    """Rectangular hyperbola"""
    return (GPP_max * PAR) / (k + PAR)

def se_calc(vals):
    """Standard error; SE=0 when n<=1."""
    vals = np.asarray(vals, dtype=float)
    vals = vals[~np.isnan(vals)]
    n = len(vals)
    if n <= 1:
        return 0.0
    return float(np.std(vals, ddof=1) / np.sqrt(n))

# Mesocosm order + group mapping
mesocosm_order = [
    "B1","B2","B3","BG1","BG2","BG3","BGS1","BGS2","BGS3",
    "GS1","GS2","GS3","G1","G2","G3"
]
meso_to_group = {
    "B1": "B", "B2": "B", "B3": "B",
    "BG1": "BG", "BG2": "BG", "BG3": "BG",
    "BGS1": "BGS", "BGS2": "BGS", "BGS3": "BGS",
    "GS1": "GS", "GS2": "GS", "GS3": "GS",
    "G1": "G", "G2": "G", "G3": "G"
}

# ======================================
# Load raw data (PAR, GPP)
# ======================================
df = pd.read_excel(data_file)
df["Date"] = pd.to_datetime(df["Date"]).dt.floor("D")

# ======================================
# Fit GPP vs PAR per (Date, Mesocosm)
# ======================================
results = []
for (date, meso), g in df.groupby(["Date", "Mesocosm"], sort=False):
    PAR = g["PAR"].to_numpy()
    GPP = g["GPP"].to_numpy()

    m = PAR <= 2200
    PAR = PAR[m]; GPP = GPP[m]
    if np.sum(~np.isnan(PAR)) <= 2:
        continue

    try:
        # bounds: GPPmax ≤ 0 ;  10 ≤ k ≤ 2200
        p0 = (-0.5, 200.0)
        bounds = ([-30, 10.0], [0.0, 2200.0])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            params, cov = curve_fit(gpp_model, PAR, GPP, p0=p0, bounds=bounds, maxfev=10000)

        gppmax_opt, k_opt = params

        # Wald p-values
        dof = max(1, len(PAR) - len(params))
        if cov is None or np.any(~np.isfinite(np.diag(cov))):
            p_gpp = np.nan; p_k = np.nan
        else:
            se = np.sqrt(np.diag(cov))
            tvals = params / se
            p_gpp = 2 * (1 - t.cdf(abs(tvals[0]), dof))
            p_k   = 2 * (1 - t.cdf(abs(tvals[1]), dof))

        # R^2 for info
        pred = gpp_model(PAR, gppmax_opt, k_opt)
        resid = GPP - pred
        ss_res = np.nansum(resid**2)
        ss_tot = np.nansum((GPP - np.nanmean(GPP))**2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else np.nan

        results.append({
            "Date": pd.Timestamp(date),
            "Mesocosm": meso,
            "Group": meso_to_group.get(meso, np.nan),
            "GPP_max": float(gppmax_opt),
            "k": float(k_opt),
            "p_value_GPP_max": float(p_gpp),
            "p_value_k": float(p_k),
            "R_squared": float(r2),
        })

    except Exception as e:
        print(f"⚠️ Fit failed for {meso} on {date}: {e}")

results_df = pd.DataFrame(results)
if results_df.empty:
    raise SystemExit("No fits produced.")
results_df["Date"] = pd.to_datetime(results_df["Date"]).dt.floor("D")

# Save raw per-campaign fits
raw_xlsx = os.path.join(output_directory, "GPPmax_k_PER_CAMPAIGN_RAW.xlsx")
results_df.to_excel(raw_xlsx, index=False)
print(f"✅ Saved raw fits: {raw_xlsx}")

# ======================================
# Significance masks & exports
# ======================================
mask_either = (results_df["p_value_GPP_max"] < 0.05) | (results_df["p_value_k"] < 0.05)
df_either = results_df[mask_either].copy()   # used for per-replicate plots & raw scatter

sig_gpp = results_df[results_df["p_value_GPP_max"] < 0.05].copy()
sig_k   = results_df[results_df["p_value_k"] < 0.05].copy()

sig_gpp.to_excel(os.path.join(output_directory, "Significant_GPPmax_ONLY.xlsx"), index=False)
sig_k.to_excel(os.path.join(output_directory, "Significant_k_ONLY.xlsx"), index=False)

# ======================================
# Per-replicate plots (free y-axis)
# ======================================
for meso in [m for m in mesocosm_order if m in df_either["Mesocosm"].unique()]:
    sub = df_either[df_either["Mesocosm"] == meso].sort_values("Date")
    if sub.empty:
        continue
    colors = plt.cm.tab20(np.linspace(0, 1, max(1, len(sub))))

    fig, ax = plt.subplots(figsize=(12, 7))

    # collect all y-values that will be plotted (for autoscale)
    y_all = []

    for j, (_, row) in enumerate(sub.iterrows()):
        d = row["Date"]
        raw = df[(df["Mesocosm"] == meso) & (df["Date"] == d)]
        PAR = raw["PAR"].to_numpy(); GPP = raw["GPP"].to_numpy()
        m = PAR <= 2200
        PAR = PAR[m]; GPP = GPP[m]

        # keep y for autoscale
        if GPP.size:
            y_all.extend(GPP.tolist())

        PAR_fit = np.linspace(0, 2200, 200)
        ax.scatter(PAR, GPP, s=40, color=colors[j % len(colors)], alpha=0.85)
        ax.plot(PAR_fit, gpp_model(PAR_fit, row["GPP_max"], row["k"]),
                color=colors[j % len(colors)],
                label=f"{d:%Y-%m-%d}: GPPmax={row['GPP_max']:.1f}{'*' if row['p_value_GPP_max']<0.05 else ''}, "
                      f"k={row['k']:.1f}{'*' if row['p_value_k']<0.05 else ''}")

    ax.set_title(meso, fontsize=16, fontweight="bold")
    ax.set_xlabel("PPFD (µmol m$^{-2}$ s$^{-1}$)", fontsize=13)
    ax.set_ylabel("GPP (µmol m$^{-2}$ s$^{-1}$)", fontsize=13)
    ax.grid(True, alpha=0.4)
    ax.legend(fontsize=12, ncol=1, frameon=False, loc="lower right",
              bbox_to_anchor=(1.15, 0.02), handlelength=1.6, columnspacing=0.8)

    # --- FREE Y-AXIS: autoscale with padding ---
    if len(y_all) > 0:
        y_all = np.asarray(y_all, dtype=float)
        y_all = y_all[np.isfinite(y_all)]
        if y_all.size:
            y_min, y_max = np.nanmin(y_all), np.nanmax(y_all)
            if y_min == y_max:
                pad = 1.0 if y_max == 0 else 0.05 * abs(y_max)
                ax.set_ylim(y_min - pad, y_max + pad)
            else:
                pad = 0.05 * (y_max - y_min)
                ax.set_ylim(y_min - pad, y_max + pad)

    out_png = os.path.join(output_directory, f"REPLICATE_{meso}_eitherSIG.png")
    plt.savefig(out_png, dpi=300, bbox_inches="tight"); plt.close()
    print(f"📊 Saved replicate plot → {out_png}")

# ======================================
# Significance-aware GROUP AVERAGING
# ======================================
def build_sigaware_panel(df_in):
    rows = []
    for (d, grp), sub in df_in.groupby(["Date", "Group"]):
        # GPPmax
        gpp_sig_vals = sub.loc[sub["p_value_GPP_max"] < 0.05, "GPP_max"].to_numpy(dtype=float)
        gpp_mean = float(np.nanmean(gpp_sig_vals)) if gpp_sig_vals.size > 0 else np.nan
        gpp_se   = se_calc(gpp_sig_vals) if gpp_sig_vals.size > 0 else np.nan

        # k
        k_sig_vals = sub.loc[sub["p_value_k"] < 0.05, "k"].to_numpy(dtype=float)
        k_mean = float(np.nanmean(k_sig_vals)) if k_sig_vals.size > 0 else np.nan
        k_se   = se_calc(k_sig_vals) if k_sig_vals.size > 0 else np.nan

        rows.append({
            "Date": pd.Timestamp(d), "Group": grp,
            "GPPmax_mean": gpp_mean, "GPPmax_se": gpp_se, "GPPmax_n": int(len(gpp_sig_vals)),
            "k_mean": k_mean, "k_se": k_se, "k_n": int(len(k_sig_vals))
        })
    out = pd.DataFrame(rows).sort_values(["Group", "Date"]).reset_index(drop=True)
    return out

panel = build_sigaware_panel(results_df)

# Save the exact panel used for all group/time-series plots
panel_xlsx = os.path.join(output_directory, "GROUP_SIGAWARE_PANEL.xlsx")
panel.to_excel(panel_xlsx, index=False)
print(f"✅ Saved significance-aware panel → {panel_xlsx}")

# ======================================
# GROUP CURVES (2×3) — both significant, show ±SE
# ======================================
dates_sorted = sorted(panel["Date"].unique())
date_colors = {pd.Timestamp(d): plt.cm.tab20(float(i % 20) / 20.0)
               for i, d in enumerate(dates_sorted)}

def plot_group_curves(ax, grp):
    sub = panel[panel["Group"] == grp].sort_values("Date")
    if sub.empty:
        ax.axis("off")
        return

    sig_mesos = [m for m, g in meso_to_group.items() if g == grp]
    df_both = results_df[
        (results_df["p_value_GPP_max"] < 0.05) &
        (results_df["p_value_k"] < 0.05) &
        (results_df["Mesocosm"].isin(sig_mesos))
    ]

    y_vals, handles, labels = [], [], []
    for _, row in sub.iterrows():
        d = pd.Timestamp(row["Date"])
        color = date_colors.get(d, "black")

        # raw scatter only when both-significant for that mesocosm/date
        for meso in sig_mesos:
            if not ((df_both["Mesocosm"] == meso) & (df_both["Date"] == d)).any():
                continue
            raw = df[(df["Mesocosm"] == meso) & (df["Date"] == d)]
            if raw.empty:
                continue
            PAR = raw["PAR"].to_numpy()
            GPP = raw["GPP"].to_numpy()
            m = PAR <= 2200
            if np.sum(m) > 2:
                ax.scatter(PAR[m], GPP[m], s=12, alpha=0.4, color=color)
                y_vals.extend(GPP[m].tolist())

        # averaged fitted curve using group mean ± SE
        if np.isfinite(row["GPPmax_mean"]) and np.isfinite(row["k_mean"]):
            PAR_fit = np.linspace(0, 2200, 200)
            fit_vals = gpp_model(PAR_fit, row["GPPmax_mean"], row["k_mean"])
            ax.plot(PAR_fit, fit_vals, color=color, lw=1.4)
            y_vals.extend(fit_vals.tolist())

            gpp_str = (f"GPPmax={row['GPPmax_mean']:.1f}±{row['GPPmax_se']:.1f}"
                       if np.isfinite(row['GPPmax_se']) else f"GPPmax={row['GPPmax_mean']:.1f}")
            k_str = (f"k={row['k_mean']:.1f}±{row['k_se']:.1f}"
                     if np.isfinite(row['k_se']) else f"k={row['k_mean']:.1f}")
            handles.append(Line2D([0],[0], color=color, lw=1.6))
            labels.append(f"{d:%b-%Y} | {gpp_str}, {k_str}")

    ax.set_title(grp, fontsize=18, fontweight='bold')
    ax.set_xlabel("PPFD (µmol m$^{-2}$ s$^{-1}$)", fontsize=14)
    ax.set_ylabel("GPP (µmol m$^{-2}$ s$^{-1}$)", fontsize=14)
    ax.set_xlim(0, 2200)

    # free y-axis autoscaling
    if y_vals:
        y_vals = np.asarray(y_vals, dtype=float)
        y_vals = y_vals[np.isfinite(y_vals)]
        if y_vals.size:
            y_min, y_max = np.nanmin(y_vals), np.nanmax(y_vals)
            pad = 0.05 * (y_max - y_min) if y_max != y_min else 1.0
            ax.set_ylim(y_min - pad, y_max + pad)

    ax.grid(True, alpha=0.4)
    ax.tick_params(axis='both', labelsize=13)
    if handles:
        ax.legend(handles, labels, fontsize=12, ncol=1, frameon=False,
                  loc='lower right', bbox_to_anchor=(1.15, 0.02),
                  handlelength=1.6, columnspacing=0.8)

# --- Draw 2×3 grid (last panel empty) ---
fig, axes = plt.subplots(3, 2, figsize=(18, 18))
fig.subplots_adjust(hspace=0.35, wspace=0.25, left=0.10, right=0.86, top=0.95, bottom=0.06)

for ax, grp in zip(axes.flat, ["B", "BG", "BGS", "GS", "G", None]):
    if grp is None:
        ax.axis("off")
        continue
    plot_group_curves(ax, grp)

out_groups = os.path.join(output_directory, "GROUP_CURVES_SIGAWARE_BOTH_SIGNIFICANT_2x3.png")
plt.savefig(out_groups, dpi=300, bbox_inches="tight")
plt.close()
print(f"✅ Saved group curves (2×3 layout, both significant, ±SE shown) → {out_groups}")

# ======================================
# REPLICATE CURVES (3×5) — both significant, show ±SE
# ======================================
dates_sorted = sorted(results_df["Date"].unique())
date_colors = {pd.Timestamp(d): plt.cm.tab20(float(i % 20) / 20.0)
               for i, d in enumerate(dates_sorted)}

def plot_replicate_curves(ax, meso):
    sub = results_df[
        (results_df["Mesocosm"] == meso) &
        (results_df["p_value_GPP_max"] < 0.05) &
        (results_df["p_value_k"] < 0.05)
    ].sort_values("Date")
    if sub.empty:
        ax.axis("off")
        return

    y_vals, handles, labels = [], [], []

    for _, row in sub.iterrows():
        d = pd.Timestamp(row["Date"])
        color = date_colors.get(d, "black")

        # --- raw data for this date ---
        raw = df[(df["Mesocosm"] == meso) & (df["Date"] == d)]
        if raw.empty:
            continue
        PAR = raw["PAR"].to_numpy()
        GPP = raw["GPP"].to_numpy()
        m = PAR <= 2200
        if np.sum(m) > 2:
            ax.scatter(PAR[m], GPP[m], s=12, alpha=0.4, color=color)
            y_vals.extend(GPP[m].tolist())

        # --- fitted curve for this replicate/date ---
        if np.isfinite(row["GPP_max"]) and np.isfinite(row["k"]):
            PAR_fit = np.linspace(0, 2200, 200)
            fit_vals = gpp_model(PAR_fit, row["GPP_max"], row["k"])
            ax.plot(PAR_fit, fit_vals, color=color, lw=1.4)
            y_vals.extend(fit_vals.tolist())

            # build label (average ± SE not applicable here → single fit)
            label_text = (f"{d:%b-%Y} | "
                          f"GPPmax={row['GPP_max']:.1f}, "
                          f"k={row['k']:.1f}")
            handles.append(Line2D([0],[0], color=color, lw=1.6))
            labels.append(label_text)

    ax.set_title(meso, fontsize=14, fontweight='bold')
    ax.set_xlabel("PPFD (µmol m$^{-2}$ s$^{-1}$)", fontsize=11)
    ax.set_ylabel("GPP (µmol m$^{-2}$ s$^{-1}$)", fontsize=11)
    ax.set_xlim(0, 2200)

    # --- free y-axis autoscaling ---
    if y_vals:
        y_vals = np.asarray(y_vals, dtype=float)
        y_vals = y_vals[np.isfinite(y_vals)]
        if y_vals.size:
            y_min, y_max = np.nanmin(y_vals), np.nanmax(y_vals)
            pad = 0.05 * (y_max - y_min) if y_max != y_min else 1.0
            ax.set_ylim(y_min - pad, y_max + pad)

    ax.grid(True, alpha=0.4)
    ax.tick_params(axis='both', labelsize=9)
    if handles:
        ax.legend(handles, labels, fontsize=8, ncol=1, frameon=False,
                  loc='lower right', bbox_to_anchor=(1.12, 0.02),
                  handlelength=1.3, columnspacing=0.8)

# --- Layout: 3 columns × 5 rows (15 replicates) ---
replicate_order = [
    "B1","B2","B3",
    "BG1","BG2","BG3",
    "BGS1","BGS2","BGS3",
    "GS1","GS2","GS3",
    "G1","G2","G3"
]

fig, axes = plt.subplots(5, 3, figsize=(18, 25))
fig.subplots_adjust(hspace=0.35, wspace=0.25,
                    left=0.08, right=0.9, top=0.95, bottom=0.06)

for ax, meso in zip(axes.flat, replicate_order):
    plot_replicate_curves(ax, meso)

out_reps = os.path.join(output_directory, "REPLICATE_CURVES_BOTH_SIGNIFICANT_3x5.png")
plt.savefig(out_reps, dpi=300, bbox_inches="tight")
plt.close()
print(f"✅ Saved replicate curves (3×5 layout, both significant) → {out_reps}")


# ======================================
# TIME-SERIES (2×3) — unified y-axis across panels
# ======================================
def ts_panel(df_panel, value_col, se_col, title, ylabel, fname):
    # --- compute global y-limits from ALL groups using value ± SE ---
    df_valid = df_panel[np.isfinite(df_panel[value_col])].copy()
    if df_valid.empty:
        print(f"⚠️ No valid data for {value_col}; skipping {fname}")
        return
    se_all = df_valid[se_col].to_numpy(dtype=float)
    se_all = np.where(np.isfinite(se_all), se_all, 0.0)

    y_all = df_valid[value_col].to_numpy(dtype=float)
    y_low  = y_all - se_all
    y_high = y_all + se_all
    y_min  = float(np.nanmin(y_low))
    y_max  = float(np.nanmax(y_high))
    pad = 0.05 * (y_max - y_min) if y_max != y_min else 1.0
    y_limits = (y_min - pad, y_max + pad)

    # --- plot with unified y-limits on every panel ---
    fig, axes = plt.subplots(2, 3, figsize=(30, 12))
    fig.subplots_adjust(hspace=0.35, wspace=0.25, left=0.06, right=0.95, top=0.90, bottom=0.10)
    start, end = datetime(2023, 2, 1), datetime(2024, 10, 31)
    order = ["B","BG","BGS","GS","G",None]

    for ax, grp in zip(axes.flat, order):
        if grp is None:
            ax.axis("off"); continue

        sub = df_panel[df_panel["Group"] == grp].sort_values("Date")
        ax.set_title(grp, fontsize=20, fontweight='bold')

        sub = sub[np.isfinite(sub[value_col])]
        if not sub.empty:
            y = sub[value_col].to_numpy(dtype=float)
            se = sub[se_col].to_numpy(dtype=float)
            se = np.where(np.isfinite(se), se, 0.0)

            ax.errorbar(sub["Date"], y, yerr=se,
                        fmt='o', linestyle='none', color='black', ecolor='black',
                        capsize=3, elinewidth=1.5, markeredgewidth=1.5)
        else:
            ax.text(0.5, 0.5, "No significant data", ha='center', va='center', transform=ax.transAxes)

        ax.set_xlim(start, end)
        ax.set_ylim(*y_limits)  # << unified y-scale applied here
        ax.grid(True, alpha=0.4)
        ax.set_ylabel(ylabel, fontsize=20); ax.set_xlabel("Date", fontsize=20)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(45); lbl.set_ha('right')

    outpath = os.path.join(output_directory, fname)
    plt.savefig(outpath, dpi=300, bbox_inches="tight"); plt.close()
    print(f"📈 Saved {title} (UNIFIED y-axis) → {outpath}")


ts_panel(panel, "GPPmax_mean", "GPPmax_se",
         "Averaged GPPmax Time Series (significance-aware)",
         "GPPmax (µmol m$^{-2}$ s$^{-1}$)",
         "TS_GPPmax_SIGAWARE.png")

ts_panel(panel, "k_mean", "k_se",
         "Averaged k Time Series (significance-aware)",
         "k (µmol m$^{-2}$ s$^{-1}$)",
         "TS_k_SIGAWARE.png")

print("🎉 Done. Y-axes are now free (autoscaled with padding) for replicate plots, group curves, and time-series.")

# ======================================
#  REPLICATE-LEVEL TIME SERIES (3×5) — UNIFIED SCALE
# ======================================

print("📈 Building replicate-level GPPmax and k time series (3×5 layout)...")

# Prepare replicate panel (use all available data)
rep_panel = results_df.copy()
rep_panel = rep_panel.sort_values(["Mesocosm", "Date"]).reset_index(drop=True)

# Define replicate and layout order
replicate_order = [
    "B1","B2","B3",
    "BG1","BG2","BG3",
    "BGS1","BGS2","BGS3",
    "GS1","GS2","GS3",
    "G1","G2","G3"
]
rep_panel["Mesocosm"] = pd.Categorical(rep_panel["Mesocosm"],
                                       categories=replicate_order,
                                       ordered=True)

# === Helper: replicate time-series plot (POINTS ONLY, SIG-ONLY, CLEAN AXIS)
def plot_replicate_ts(df, value_col, title, ylabel, fname):
    # Keep only significant rows for the chosen value
    if value_col == "GPP_max":
        df_sig = df[df["p_value_GPP_max"] < 0.05].copy()
    elif value_col == "k":
        df_sig = df[df["p_value_k"] < 0.05].copy()
        # EXTRA FILTER: drop very large k values
        df_sig = df_sig[df_sig["k"] <= 1500].copy()
    else:
        df_sig = df.copy()

    if df_sig.empty or not np.isfinite(df_sig[value_col]).any():
        print(f"⚠️ No usable data for {value_col} after filtering; skipping {fname}")
        return

    # Unified Y limits (from filtered data)
    y_all = df_sig[value_col].to_numpy(dtype=float)
    y_all = y_all[np.isfinite(y_all)]
    y_min, y_max = np.nanmin(y_all), np.nanmax(y_all)
    pad = 0.05 * (y_max - y_min) if y_max != y_min else 1.0
    y_limits = (y_min - pad, y_max + pad)

    # Unified X limits (from filtered data)
    start = pd.to_datetime(df_sig["Date"].min()) - pd.Timedelta(days=7)
    end   = pd.to_datetime(df_sig["Date"].max()) + pd.Timedelta(days=7)

    # Date ticks: every 2 months, YYYY-MM
    locator = mdates.MonthLocator(interval=2)
    formatter = mdates.DateFormatter('%Y-%m')

    nrows, ncols = 5, 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 18))
    fig.subplots_adjust(hspace=0.22, wspace=0.18,
                        left=0.08, right=0.96, top=0.90, bottom=0.08)

    for idx, (ax, meso) in enumerate(zip(axes.flat, replicate_order)):
        sub = df_sig[df_sig["Mesocosm"] == meso].sort_values("Date")
        if sub.empty:
            ax.axis("off")
            continue

        # Points only
        ax.scatter(sub["Date"], sub[value_col], color="black", s=35, alpha=0.9)

        ax.set_title(meso, fontsize=12, fontweight="bold", pad=4)
        ax.set_ylim(*y_limits)
        ax.set_xlim(start, end)
        ax.grid(True, alpha=0.35)

        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)

        is_last_row = idx >= (nrows - 1) * ncols
        is_first_col = (idx % ncols == 0)

        if is_last_row:
            for lbl in ax.get_xticklabels():
                lbl.set_rotation(45); lbl.set_ha('right')
            ax.set_xlabel("Date", fontsize=10)
        else:
            ax.tick_params(axis="x", labelbottom=False)

        if is_first_col:
            ax.set_ylabel(ylabel, fontsize=10)
            ax.tick_params(axis="y", labelsize=9)
        else:
            ax.tick_params(axis="y", labelleft=False)

    fig.suptitle(title, fontsize=18, fontweight="bold")
    out_path = os.path.join(output_directory, fname)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"✅ Saved replicate time series (filtered & significant only) → {out_path}")

# === Plot 1: GPPmax (GPP_max column) ===
plot_replicate_ts(rep_panel, "GPP_max",
                  "GPPmax Time Series — per Replicate (3×5, unified scale)",
                  "GPPmax (µmol m$^{-2}$ s$^{-1}$)",
                  "TS_GPPmax_REPLICATES_UNIFIED.png")

# === Plot 2: k (k column) ===
plot_replicate_ts(rep_panel, "k",
                  "k Time Series — per Replicate (3×5, unified scale)",
                  "k (µmol m$^{-2}$ s$^{-1}$)",
                  "TS_k_REPLICATES_UNIFIED.png")

print("🎉 Replicate-level time series completed with unified y-axis scale for all 15 panels.")

# ======================================
# INDIVIDUAL CAMPAIGN PLOTS — BGS3 only (both significant, custom colors)
# ======================================
meso_target = "BGS3"
output_subdir = os.path.join(output_directory, f"{meso_target}_campaigns")
os.makedirs(output_subdir, exist_ok=True)

# Filter both-significant rows for BGS3
bgs3_df = results_df[
    (results_df["Mesocosm"] == meso_target) &
    (results_df["p_value_GPP_max"] < 0.05) &
    (results_df["p_value_k"] < 0.05)
].sort_values("Date")

if bgs3_df.empty:
    print(f"⚠️ No both-significant fits found for {meso_target}.")
else:
    print(f"✅ {len(bgs3_df)} significant campaign fits found for {meso_target}.")
    for _, row in bgs3_df.iterrows():
        date = pd.Timestamp(row["Date"])
        date_str = date.strftime("%Y-%m-%d")

        # --- extract measured data for that date ---
        raw = df[(df["Mesocosm"] == meso_target) & (df["Date"] == date)]
        if raw.empty:
            continue
        PAR = raw["PAR"].to_numpy()
        GPP = raw["GPP"].to_numpy()
        m = PAR <= 2200
        PAR = PAR[m]; GPP = GPP[m]

        # --- fitted curve ---
        PAR_fit = np.linspace(0, 2200, 200)
        fit_vals = gpp_model(PAR_fit, row["GPP_max"], row["k"])

        # --- custom colors ---
        peat_brown = "black"     # curve color
        sphagnum_green = "black" # points color

        # --- plot ---
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(PAR_fit, fit_vals, color=peat_brown, lw=2.2, label=f"Fit (R²={row['R_squared']:.2f})")
        ax.scatter(PAR, GPP, color=sphagnum_green, s=40, alpha=1, label="Measured GPP")
        
        # --- annotations ---
        ax.set_title(f"{meso_target} | {date_str}", fontsize=16, fontweight="bold")
        ax.set_xlabel("PPFD (µmol m$^{-2}$ s$^{-1}$)", fontsize=13)
        ax.set_ylabel("GPP (µmol m$^{-2}$ s$^{-1}$)", fontsize=13)
        ax.text(0.05, 0.93,
                f"GPPmax = {row['GPP_max']:.2f}\n"
                f"k = {row['k']:.1f}",
                transform=ax.transAxes, fontsize=12, va='top', ha='left',
                bbox=dict(facecolor='white', edgecolor='gray', alpha=0.7))

        ax.grid(True, alpha=0.4)
        ax.legend(frameon=False, fontsize=12)
        ax.set_xlim(0, 2200)

        # --- autoscale Y-axis ---
        y_all = np.concatenate([GPP, fit_vals])
        y_min, y_max = np.nanmin(y_all), np.nanmax(y_all)
        pad = 0.05 * (y_max - y_min) if y_max != y_min else 1.0
        ax.set_ylim(y_min - pad, y_max + pad)

        # --- save figure ---
        out_file = os.path.join(output_subdir, f"{meso_target}_{date_str}_fit.png")
        plt.savefig(out_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"📊 Saved → {out_file}")

# ======================================
# 2×5 PANEL — Column 1: GPPmax, Column 2: k; Rows: B, BG, BGS, GS, G
# significant-only, k<=1500, fixed axes, legends top-right, circle markers
# ======================================
def plot_treatment_gpp_k_panel_2x5(df, fname="TS_GPPmax_k_2x5_TREATMENTS_FIXED.png"):
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import pandas as pd
    import numpy as np
    import os

    # --- keep only significant points ---
    gpp_sig = df[df["p_value_GPP_max"] < 0.05].copy()
    k_sig   = df[(df["p_value_k"] < 0.05) & (df["k"] <= 1500)].copy()

    if gpp_sig.empty and k_sig.empty:
        print("⚠️ No significant data available.")
        return

    treatments = ["B", "BG", "BGS", "GS", "G"]
    replicates = ["1", "2", "3"]

    # =====================================================
    #  GREYSCALE + SHAPES (replicate ID)
    # =====================================================
    rep_colors = {
        "1": "black",       # darkest
        "2": "black",     # medium grey
        "3": "black",   # light grey
    }
    rep_markers = {
        "1": "o",  # circle
        "2": "s",  # square
        "3": "^",  # triangle
    }

    # --- fixed X axis: Feb-2023 .. Oct-2024, ticks every 5 months ---
    x_start = pd.Timestamp(2023, 2, 1)
    x_end   = pd.Timestamp(2024, 10, 31)
    xticks  = pd.date_range(start=x_start, end=x_end, freq="5MS")  # every 5 months
    xfmt    = mdates.DateFormatter('%Y-%m')

    # --- fixed Y limits ---
    gpp_ylim = (-32, 2)     # for GPPmax
    k_ylim   = (-50, 1550)  # for k

    # --- figure/grid (rows = treatments) ---
    fig, axes = plt.subplots(5, 2, figsize=(16, 22), sharex=False)
    fig.subplots_adjust(hspace=0.35, wspace=0.25,
                        left=0.08, right=0.97, top=0.93, bottom=0.07)

    for i, trt in enumerate(treatments):
        # ------------------------------------------------
        # LEFT column: GPPmax
        # ------------------------------------------------
        ax_g = axes[i, 0]
        sub_gpp = gpp_sig[gpp_sig["Group"] == trt].copy()
        has_any = False

        if not sub_gpp.empty:
            for rep in replicates:
                meso = f"{trt}{rep}"
                rep_data = sub_gpp[sub_gpp["Mesocosm"] == meso].sort_values("Date")
                if not rep_data.empty:
                    ax_g.scatter(
                        rep_data["Date"],
                        rep_data["GPP_max"],
                        color=rep_colors.get(rep, "gray"),
                        s=46,
                        alpha=0.9,
                        marker=rep_markers.get(rep, "o"),
                        label=meso,
                        edgecolor="black" if rep != "3" else "none"  # tiny contrast
                    )
                    has_any = True

        if not has_any:
            ax_g.text(0.5, 0.5, "No significant data",
                      ha="center", va="center",
                      transform=ax_g.transAxes, fontsize=10)

        # treatment label (left top)
        ax_g.text(0.02, 0.94, trt,
                  transform=ax_g.transAxes,
                  fontsize=13, fontweight="bold",
                  va="top", ha="left")

        ax_g.set_ylim(*gpp_ylim)
        ax_g.set_xlim(x_start, x_end)
        ax_g.set_xticks(xticks)
        ax_g.xaxis.set_major_formatter(xfmt)
        ax_g.tick_params(axis="x", labelrotation=45)
        ax_g.grid(True, alpha=0.35)

        # legend: show only used ones, no frame
        handles_g, labels_g = ax_g.get_legend_handles_labels()
        if handles_g:
            ax_g.legend(handles_g, labels_g,
                        fontsize=9, frameon=False,
                        loc="upper right", ncol=1)

        ax_g.set_ylabel("GPPmax (µmol m$^{-2}$ s$^{-1}$)", fontsize=11)

        # ------------------------------------------------
        # RIGHT column: k
        # ------------------------------------------------
        ax_k = axes[i, 1]
        sub_k = k_sig[k_sig["Group"] == trt].copy()
        has_any = False

        if not sub_k.empty:
            for rep in replicates:
                meso = f"{trt}{rep}"
                rep_data = sub_k[sub_k["Mesocosm"] == meso].sort_values("Date")
                if not rep_data.empty:
                    ax_k.scatter(
                        rep_data["Date"],
                        rep_data["k"],
                        color=rep_colors.get(rep, "gray"),
                        s=42,
                        alpha=0.9,
                        marker=rep_markers.get(rep, "o"),
                        label=meso,
                        edgecolor="black" if rep != "3" else "none"
                    )
                    has_any = True

        if not has_any:
            ax_k.text(0.5, 0.5, "No significant data",
                      ha="center", va="center",
                      transform=ax_k.transAxes, fontsize=10)

        ax_k.set_ylim(*k_ylim)
        ax_k.set_xlim(x_start, x_end)
        ax_k.set_xticks(xticks)
        ax_k.xaxis.set_major_formatter(xfmt)
        ax_k.tick_params(axis="x", labelrotation=45)
        ax_k.grid(True, alpha=0.35)

        handles_k, labels_k = ax_k.get_legend_handles_labels()
        if handles_k:
            ax_k.legend(handles_k, labels_k,
                        fontsize=9, frameon=False,
                        loc="upper right", ncol=1)

        ax_k.set_ylabel("k (µmol m$^{-2}$ s$^{-1}$)", fontsize=11)

    # save
    out_path = os.path.join(output_directory, fname) if 'output_directory' in globals() else fname
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"✅ Saved 2×5 treatment-level panel (greyscale + shapes) → {out_path}")

plot_treatment_gpp_k_panel_2x5(rep_panel)
