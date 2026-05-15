#!/usr/bin/env python3
"""
3-way comparison: opt_cluster_1 vs no_opt_1 vs literature
SISO (spin-orbit coupled) 4f energy levels of Er3+ in LiYF4
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

base = "/data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+"

# ============================================================
# Load computed data
# ============================================================
opt_data = np.loadtxt(f"{base}/opt_cluster_1/rohf_casscf_nevpt2_siso/LiYF4:Er3+_mag.txt")
opt_unique = np.unique(np.round(opt_data, 4))

noopt_data = np.loadtxt(f"{base}/no_opt_1/cahf_casscf_nevpt2_siso/LiYF4:Er3+_mag.txt")
noopt_unique = np.unique(np.round(noopt_data, 4))

# ============================================================
# Load literature data
# ============================================================
lit_file = f"{base}/literature_Er_LiYF4.txt"
lit_levels = []
lit_labels = []
with open(lit_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        manifold = parts[0]
        energy = float(parts[2])
        lit_levels.append(energy)
        # Short label
        label = manifold.replace('4', '⁴').replace('2', '²')
        lit_labels.append(label)
lit_levels = np.array(lit_levels)

# ============================================================
# Plot
# ============================================================
fig, ax = plt.subplots(1, 1, figsize=(12, 14))

x_opt, x_noopt, x_lit = 1.0, 2.0, 3.0
colors = {'opt': '#E74C3C', 'noopt': '#3498DB', 'lit': '#2ECC71'}

# --- opt_cluster_1 ---
for e in opt_unique:
    if e <= 21000:
        ax.hlines(e, x_opt - 0.32, x_opt + 0.32, colors=colors['opt'], linewidth=2.5, zorder=4)

# --- no_opt_1 ---
for e in noopt_unique:
    if e <= 21000:
        ax.hlines(e, x_noopt - 0.32, x_noopt + 0.32, colors=colors['noopt'], linewidth=2.5, zorder=4)

# --- Literature ---
for e, lbl in zip(lit_levels, lit_labels):
    ax.hlines(e, x_lit - 0.32, x_lit + 0.32, colors=colors['lit'], linewidth=2.5, zorder=4)

# --- Connect corresponding levels ---
# Match by sorting all levels and connecting nearest neighbors
for i in range(len(opt_unique)):
    if opt_unique[i] > 21000: break
    # Find nearest no-opt level
    idx_noopt = np.argmin(np.abs(noopt_unique - opt_unique[i]))
    if noopt_unique[idx_noopt] <= 21000:
        ax.plot([x_opt + 0.32, x_noopt - 0.32], [opt_unique[i], noopt_unique[idx_noopt]],
                '-', color='gray', alpha=0.25, linewidth=0.5, zorder=1)

for i in range(len(noopt_unique)):
    if noopt_unique[i] > 21000: break
    idx_lit = np.argmin(np.abs(lit_levels - noopt_unique[i]))
    if lit_levels[idx_lit] <= 21000:
        ax.plot([x_noopt + 0.32, x_lit - 0.32], [noopt_unique[i], lit_levels[idx_lit]],
                '-', color='gray', alpha=0.2, linewidth=0.5, zorder=1)

# --- J-manifold labels on the right ---
manifold_boundaries = {
    '⁴I₁₅/₂': (0, 500),
    '⁴I₁₃/₂': (6400, 6800),
    '⁴I₁₁/₂': (10100, 10400),
    '⁴I₉/₂': (12200, 12800),
    '⁴F₉/₂': (15200, 15600),
    '⁴S₃/₂': (18300, 18600),
    '²H₁₁/₂': (19000, 19500),
    '⁴F₇/₂': (20400, 20800),
}
for label, (ylo, yhi) in manifold_boundaries.items():
    ax.text(x_lit + 0.55, (ylo + yhi) / 2, label, fontsize=9, va='center',
            fontstyle='italic', color='#555555')

# --- Styling ---
ax.set_xlim(0.3, 4.2)
ax.set_xticks([x_opt, x_noopt, x_lit])
ax.set_xticklabels([
    'opt_1 (ErF₈, ROHF, opt)',
    'no_opt_1 (ErF₈, CAHF, pre-opt)',
    'Experiment\n(Literature)'
], fontsize=11)
ax.set_ylabel('Energy (cm⁻¹)', fontsize=13)
ax.set_title('Er³⁺ 4f Crystal-Field Energy Levels in LiYF₄\nNEVPT2+CASSCF+SISO vs Experiment', 
             fontsize=14, fontweight='bold')
ax.grid(axis='y', linestyle='--', alpha=0.3)
ax.set_ylim(-50, 21000)

# Legend
from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([0], [0], color=colors['opt'], lw=3, label='opt_1 (ErF₈, ROHF, opt)'),
    Line2D([0], [0], color=colors['noopt'], lw=3, label='no_opt_1 (ErF₈, CAHF, pre-opt)'),
    Line2D([0], [0], color=colors['lit'], lw=3, label='Experiment'),
], loc='upper right', fontsize=10)

plt.tight_layout()
outfile = f'{base}/opt_cluster_1/rohf_casscf_nevpt2_siso/4f_levels_3way_comparison.png'
plt.savefig(outfile, dpi=200, bbox_inches='tight')
print(f"Saved: {outfile}")

# ============================================================
# RMS error table
# ============================================================
print(f"\n{'='*70}")
print(f"  RMS Error vs Experiment (cm⁻¹)")
print(f"{'='*70}")
print(f"{'J-manifold':>12s}  {'opt-cluster_1':>14s}  {'no-opt_1':>14s}")

manifold_groups = {}
for i, lbl in enumerate(lit_labels):
    m = lbl.replace('⁴', '4').replace('²', '2')
    if m not in manifold_groups:
        manifold_groups[m] = []
    manifold_groups[m].append(lit_levels[i])

for m_name, lit_vals in manifold_groups.items():
    lit_arr = np.array(lit_vals)
    # Find nearest computed levels
    opt_match = np.array([opt_unique[np.argmin(np.abs(opt_unique - v))] for v in lit_arr])
    noopt_match = np.array([noopt_unique[np.argmin(np.abs(noopt_unique - v))] for v in lit_arr])
    rms_opt = np.sqrt(np.mean((opt_match - lit_arr)**2))
    rms_noopt = np.sqrt(np.mean((noopt_match - lit_arr)**2))
    print(f"{m_name:>12s}  {rms_opt:14.1f}  {rms_noopt:14.1f}")

# Overall RMS (within 21000 cm-1)
opt_all = np.array([opt_unique[np.argmin(np.abs(opt_unique - v))] for v in lit_levels])
noopt_all = np.array([noopt_unique[np.argmin(np.abs(noopt_unique - v))] for v in lit_levels])
rms_opt_all = np.sqrt(np.mean((opt_all - lit_levels)**2))
rms_noopt_all = np.sqrt(np.mean((noopt_all - lit_levels)**2))
print(f"{'─'*40}")
print(f"{'OVERALL':>12s}  {rms_opt_all:14.1f}  {rms_noopt_all:14.1f}")
