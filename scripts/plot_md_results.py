#!/usr/bin/env python3
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set publication style font size
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14
})

# Paths
BASE_DIR = "/media/hp/Data/Hosni/openmm_windows_setup"
MD_SIM_DIR = os.path.join(BASE_DIR, "AMR_Work/results/md_simulation")
TIMESERIES_FILE = os.path.join(MD_SIM_DIR, "md_analysis_timeseries.csv")
RMSF_FILE = os.path.join(MD_SIM_DIR, "md_analysis_rmsf.csv")

# Output dirs
PLOTS_DIRS = [
    os.path.join(BASE_DIR, "AMR_Work/paper1/plots"),
    os.path.join(BASE_DIR, "AMR_Work/manuscript/plots"),
    os.path.join(BASE_DIR, "AMR_Work/results/figures")
]

for d in PLOTS_DIRS:
    os.makedirs(d, exist_ok=True)

if not os.path.exists(TIMESERIES_FILE) or not os.path.exists(RMSF_FILE):
    print(f"[ERROR] Required CSV files not found in {MD_SIM_DIR}")
    sys.exit(1)

# Load data
df_ts = pd.read_csv(TIMESERIES_FILE)
df_rmsf = pd.read_csv(RMSF_FILE)

# 1. Figure 10a: Backbone RMSD and Rg Dual Plot
print("Plotting Backbone RMSD and Radius of Gyration...")
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)

# RMSD Plot
ax1.plot(df_ts['Time_ns'], df_ts['RMSD_Angstrom'], color='#2b5c8f', linewidth=1.2, label='Backbone RMSD')
# Smooth curve
if len(df_ts) > 100:
    smooth_rmsd = df_ts['RMSD_Angstrom'].rolling(window=100, center=True).mean()
    ax1.plot(df_ts['Time_ns'], smooth_rmsd, color='#e05a47', linewidth=1.5, label='Running Average (1 ns)')
ax1.set_ylabel('RMSD (Å)', fontweight='bold')
ax1.grid(True, linestyle='--', alpha=0.5)
ax1.legend(loc='lower right')
ax1.set_title('Protein structural evolution & dynamics', fontweight='bold', fontsize=12, loc='left')

# Rg Plot
ax2.plot(df_ts['Time_ns'], df_ts['Rg_Angstrom'], color='#108040', linewidth=1.2, label='Radius of Gyration ($R_g$)')
if len(df_ts) > 100:
    smooth_rg = df_ts['Rg_Angstrom'].rolling(window=100, center=True).mean()
    ax2.plot(df_ts['Time_ns'], smooth_rg, color='#d0a010', linewidth=1.5, label='Running Average (1 ns)')
ax2.set_xlabel('Time (ns)', fontweight='bold')
ax2.set_ylabel('Radius of Gyration $R_g$ (Å)', fontweight='bold')
ax2.grid(True, linestyle='--', alpha=0.5)
ax2.legend(loc='upper right')

plt.tight_layout()

for d in PLOTS_DIRS:
    plt.savefig(os.path.join(d, "figure10_md_rmsd_rg.png"), dpi=300)
    plt.savefig(os.path.join(d, "figure10_md_rmsd_rg.pdf"))
plt.close()

# 2. Figure 10b: C-alpha RMSF Plot
print("Plotting Residue-wise C-alpha RMSF...")
plt.figure(figsize=(9, 4.5))
plt.plot(df_rmsf['Residue'], df_rmsf['RMSF_Angstrom'], color='#503080', linewidth=1.5, label='C-alpha RMSF')
plt.grid(True, linestyle='--', alpha=0.5)

# Highlight active site catalytic residues: His125, Asp187, Asp188
catalytic_res = [125, 187, 188]
catalytic_labels = ["His125", "Asp187", "Asp188"]
colors = ['#d01010', '#10b010', '#1010d0']

for res, label, color in zip(catalytic_res, catalytic_labels, colors):
    val = df_rmsf.loc[df_rmsf['Residue'] == res, 'RMSF_Angstrom'].values[0]
    plt.scatter(res, val, color=color, s=50, zorder=5)
    plt.annotate(f"{label} ({val:.2f} Å)", (res, val), textcoords="offset points", 
                 xytext=(0,10), ha='center', fontweight='bold', color=color,
                 arrowprops=dict(arrowstyle="->", color=color, lw=1.0))

plt.xlabel('Residue Number', fontweight='bold')
plt.ylabel('Fluctuation RMSF (Å)', fontweight='bold')
plt.title('Local Residue-wise Fluctuations (C-alpha RMSF) with Catalytic Triad stability', fontweight='bold', fontsize=12, loc='left')
plt.xlim(1, 351)
plt.ylim(0, max(df_rmsf['RMSF_Angstrom']) + 0.5)

plt.tight_layout()

for d in PLOTS_DIRS:
    plt.savefig(os.path.join(d, "figure10_md_rmsf.png"), dpi=300)
    plt.savefig(os.path.join(d, "figure10_md_rmsf.pdf"))
plt.close()

print("MD figures generated successfully in paper1, manuscript, and results/figures directories.")
