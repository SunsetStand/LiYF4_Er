#!/usr/bin/env python3
"""
4-way comparison: opt_1 vs no_opt_1 vs no_opt_2 vs experiment
Er3+ 4f SISO energy levels in LiYF4
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

base = "/data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+"

# ============================================================
# Load data
# ============================================================
def load_unique(path):
    data = np.loadtxt(path)
    return np.unique(np.round(data, 4))

opt = load_unique(f"{base}/opt_cluster_1/rohf_casscf_nevpt2_siso/LiYF4:Er3+_mag.txt")
noopt1 = load_unique(f"{base}/no_opt_1/cahf_casscf_nevpt2_siso/LiYF4:Er3+_mag.txt")
noopt2 = load_unique(f"{base}/no_opt_2/rohf_casscf/LiYF4:Er3+_mag.txt")

lit_manifolds, lit_levels = [], []
with open(f"{base}/literature_Er_LiYF4.txt") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'): continue
        parts = line.split()
        lit_manifolds.append(parts[0])
        lit_levels.append(float(parts[2]))
lit_levels = np.array(lit_levels)

# ============================================================
# Per-manifold detail chart
# ============================================================
manifold_order = ['4I15/2', '4I13/2', '4I11/2', '4I9/2', '4F9/2', '4S3/2', '2H11/2', '4F7/2']
manifold_label = {'4I15/2': '⁴I₁₅/₂', '4I13/2': '⁴I₁₃/₂', '4I11/2': '⁴I₁₁/₂',
                  '4I9/2': '⁴I₉/₂', '4F9/2': '⁴F₉/₂', '4S3/2': '⁴S₃/₂',
                  '2H11/2': '²H₁₁/₂', '4F7/2': '⁴F₇/₂'}
colors = {'opt': '#E74C3C', 'noopt1': '#3498DB', 'noopt2': '#F39C12'}

n_manifolds = len(manifold_order)
fig, axes = plt.subplots(n_manifolds, 1, figsize=(11, 2.2 * n_manifolds),
                         gridspec_kw={'hspace': 0.18})
fig.suptitle('Er³⁺ in LiYF₄ — Crystal-Field Levels vs Experiment\n'
             'opt_1 (ErF₈+ROHF+opt) | no_opt_1 (ErF₈+CAHF) | no_opt_2 (25-atom+ROHF)',
             fontsize=14, fontweight='bold', y=0.996)

for ax_idx, m_name in enumerate(manifold_order):
    ax = axes[ax_idx]
    mask = np.array([m == m_name for m in lit_manifolds])
    lit_m = lit_levels[mask]
    n_levels = len(lit_m)
    
    opt_match = np.array([opt[np.argmin(np.abs(opt - v))] for v in lit_m])
    noopt1_match = np.array([noopt1[np.argmin(np.abs(noopt1 - v))] for v in lit_m])
    noopt2_match = np.array([noopt2[np.argmin(np.abs(noopt2 - v))] for v in lit_m])
    
    x_positions = np.arange(n_levels)
    width = 0.22
    
    diff_opt = opt_match - lit_m
    diff_no1 = noopt1_match - lit_m
    diff_no2 = noopt2_match - lit_m
    
    ax.bar(x_positions - width, diff_opt, width, color=colors['opt'], alpha=0.85, zorder=3)
    ax.bar(x_positions, diff_no1, width, color=colors['noopt1'], alpha=0.85, zorder=3)
    ax.bar(x_positions + width, diff_no2, width, color=colors['noopt2'], alpha=0.85, zorder=3)
    ax.axhline(0, color='black', linewidth=0.8, zorder=2)
    
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f'{v:.0f}' for v in lit_m], fontsize=7.5)
    ax.set_ylabel('Δ cm⁻¹', fontsize=9)
    ax.text(-0.08, 0.5, manifold_label[m_name], transform=ax.transAxes, fontsize=11,
            va='center', ha='right', fontweight='bold', fontstyle='italic')
    
    rms_o = np.sqrt(np.mean(diff_opt**2))
    rms_n1 = np.sqrt(np.mean(diff_no1**2))
    rms_n2 = np.sqrt(np.mean(diff_no2**2))
    ax.text(0.98, 0.95, f'RMS: opt₁={rms_o:.0f}  no₁={rms_n1:.0f}  no₂={rms_n2:.0f}',
            transform=ax.transAxes, fontsize=7.5, ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    
    ymax = max(abs(diff_opt).max(), abs(diff_no1).max(), abs(diff_no2).max()) * 1.3
    ymax = max(ymax, 25)
    ax.set_ylim(-ymax, ymax)

axes[-1].set_xlabel('Experimental level energy (cm⁻¹)', fontsize=10)

from matplotlib.lines import Line2D
fig.legend(handles=[
    Line2D([0],[0], color=colors['opt'], lw=4, label='opt_1 (ErF₈, ROHF, optimized)'),
    Line2D([0],[0], color=colors['noopt1'], lw=4, label='no_opt_1 (ErF₈, CAHF, pre-opt)'),
    Line2D([0],[0], color=colors['noopt2'], lw=4, label='no_opt_2 (25-atom, ROHF, pre-opt)'),
], loc='lower center', ncol=3, fontsize=9.5, bbox_to_anchor=(0.5, -0.02))

plt.tight_layout(rect=[0.08, 0.04, 1, 0.99])
outfile = f'{base}/opt_cluster_1/rohf_casscf_nevpt2_siso/4f_levels_4way_detail.png'
plt.savefig(outfile, dpi=200, bbox_inches='tight')
print(f"Saved: {outfile}")

# ============================================================
# RMS comparison table
# ============================================================
print(f"\n{'='*85}")
print(f"  RMS Error vs Experiment (cm⁻¹)")
print(f"{'='*85}")
print(f"{'J-manifold':>12s}  {'opt_1':>10s}  {'no_opt_1':>10s}  {'no_opt_2':>10s}  {'best':>8s}")

rms_all = {'opt_1': [], 'no_opt_1': [], 'no_opt_2': []}
for m_name in manifold_order:
    mask = np.array([m == m_name for m in lit_manifolds])
    lit_m = lit_levels[mask]
    rms = {}
    for label, data in [('opt_1', opt), ('no_opt_1', noopt1), ('no_opt_2', noopt2)]:
        match = np.array([data[np.argmin(np.abs(data - v))] for v in lit_m])
        r = np.sqrt(np.mean((match - lit_m)**2))
        rms[label] = r
        rms_all[label].append(r)
    best_label = min(rms, key=rms.get)
    print(f"{m_name:>12s}  {rms['opt_1']:10.1f}  {rms['no_opt_1']:10.1f}  {rms['no_opt_2']:10.1f}  {best_label:>8s}")

print(f"{'─'*55}")
for label in ['opt_1', 'no_opt_1', 'no_opt_2']:
    overall = np.sqrt(np.mean(np.array(rms_all[label])**2))
    print(f"  {label} OVERALL = {overall:.1f}")

# Count "best" per manifold
print(f"\n  Best per manifold: ", end='')
best_counts = {'opt_1': 0, 'no_opt_1': 0, 'no_opt_2': 0}
for m_name in manifold_order:
    mask = np.array([m == m_name for m in lit_manifolds])
    lit_m = lit_levels[mask]
    best = min([
        ('opt_1', np.sqrt(np.mean((np.array([opt[np.argmin(np.abs(opt-v))] for v in lit_m])-lit_m)**2))),
        ('no_opt_1', np.sqrt(np.mean((np.array([noopt1[np.argmin(np.abs(noopt1-v))] for v in lit_m])-lit_m)**2))),
        ('no_opt_2', np.sqrt(np.mean((np.array([noopt2[np.argmin(np.abs(noopt2-v))] for v in lit_m])-lit_m)**2)))
    ], key=lambda x: x[1])[0]
    best_counts[best] += 1
print(f"opt_1={best_counts['opt_1']}, no_opt_1={best_counts['no_opt_1']}, no_opt_2={best_counts['no_opt_2']}")
