#!/usr/bin/env python3
"""
Per-manifold detailed comparison: opt vs no-opt vs experiment
Side-by-side level-by-level within each J-manifold
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

base = "/data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+"

# Load computed data
opt_raw = np.loadtxt(f"{base}/opt_cluster_1/rohf_casscf_nevpt2_siso/LiYF4:Er3+_mag.txt")
opt = np.unique(np.round(opt_raw, 4))

noopt_raw = np.loadtxt(f"{base}/no_opt_1/cahf_casscf_nevpt2_siso/LiYF4:Er3+_mag.txt")
noopt = np.unique(np.round(noopt_raw, 4))

# Load literature with manifold labels
lit_levels_all = []
lit_labels_all = []
lit_manifolds = []
with open(f"{base}/literature_Er_LiYF4.txt") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'): continue
        parts = line.split()
        lit_manifolds.append(parts[0])
        lit_levels_all.append(float(parts[2]))
        lit_labels_all.append(f"{parts[0]}({parts[1]})")
lit_levels_all = np.array(lit_levels_all)

# Group by manifold
manifold_order = ['4I15/2', '4I13/2', '4I11/2', '4I9/2', '4F9/2', '4S3/2', '2H11/2', '4F7/2']
manifold_label = {'4I15/2': '⁴I₁₅/₂', '4I13/2': '⁴I₁₃/₂', '4I11/2': '⁴I₁₁/₂',
                  '4I9/2': '⁴I₉/₂', '4F9/2': '⁴F₉/₂', '4S3/2': '⁴S₃/₂',
                  '2H11/2': '²H₁₁/₂', '4F7/2': '⁴F₇/₂'}

manifold_ranges = {}
for m_name in manifold_order:
    mask = np.array([m == m_name for m in lit_manifolds])
    levels = lit_levels_all[mask]
    manifold_ranges[m_name] = (levels.min() - 30, levels.max() + 30)

# ============================================================
# Plot: subplot per manifold
# ============================================================
n_manifolds = len(manifold_order)
fig, axes = plt.subplots(n_manifolds, 1, figsize=(10, 2.0 * n_manifolds),
                         gridspec_kw={'hspace': 0.15})
fig.suptitle('Er³⁺ in LiYF₄ — Crystal-Field Levels: opt vs no-opt vs Experiment',
             fontsize=15, fontweight='bold', y=0.995)

for ax_idx, m_name in enumerate(manifold_order):
    ax = axes[ax_idx]
    mask = np.array([m == m_name for m in lit_manifolds])
    lit_m = lit_levels_all[mask]
    n_levels = len(lit_m)
    
    # For each experimental level, find nearest computed level
    x_positions = np.arange(n_levels)
    width = 0.25
    
    # --- Bar chart: difference from experiment ---
    opt_match = np.array([opt[np.argmin(np.abs(opt - v))] for v in lit_m])
    noopt_match = np.array([noopt[np.argmin(np.abs(noopt - v))] for v in lit_m])
    
    diff_opt = opt_match - lit_m
    diff_noopt = noopt_match - lit_m
    
    bars1 = ax.bar(x_positions - width, diff_opt, width, color='#E74C3C', alpha=0.85,
                   label='opt_1 − exp', zorder=3)
    bars2 = ax.bar(x_positions, diff_noopt, width, color='#3498DB', alpha=0.85,
                   label='no-opt_1 − exp', zorder=3)
    ax.axhline(0, color='black', linewidth=0.8, zorder=2)
    
    # Value labels on bars
    for bar in bars1:
        h = bar.get_height()
        va = 'bottom' if h >= 0 else 'top'
        ax.text(bar.get_x() + bar.get_width()/2, h + (2 if h>=0 else -2),
                f'{h:+.0f}', ha='center', va=va, fontsize=7, color='#C0392B')
    for bar in bars2:
        h = bar.get_height()
        va = 'bottom' if h >= 0 else 'top'
        ax.text(bar.get_x() + bar.get_width()/2, h + (2 if h>=0 else -2),
                f'{h:+.0f}', ha='center', va=va, fontsize=7, color='#2471A3')
    
    # X-axis labels: exp energy
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f'{v:.0f}' for v in lit_m], fontsize=8)
    
    # Manifold label on left
    ax.text(-0.5, 0.5, manifold_label[m_name], transform=ax.transAxes, fontsize=12,
            va='center', ha='right', fontweight='bold', fontstyle='italic')
    
    # y-label
    ax.set_ylabel('Δ cm⁻¹', fontsize=9)
    
    # RMS annotation
    rms_opt = np.sqrt(np.mean(diff_opt**2))
    rms_noopt = np.sqrt(np.mean(diff_noopt**2))
    ax.text(0.98, 0.95, f'RMS opt={rms_opt:.0f} | no-opt={rms_noopt:.0f}',
            transform=ax.transAxes, fontsize=8, ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))
    
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    
    # Auto y-lim
    ymax = max(abs(diff_opt).max(), abs(diff_noopt).max()) * 1.3
    ymax = max(ymax, 20)
    ax.set_ylim(-ymax, ymax)

# Shared legend at bottom
handles = [plt.Rectangle((0,0),1,1,fc='#E74C3C', alpha=0.85),
           plt.Rectangle((0,0),1,1,fc='#3498DB', alpha=0.85)]
labels = ['opt_1 (ErF₈, ROHF, opt)', 'no_opt_1 (ErF₈, CAHF, pre-opt)']
fig.legend(handles, labels, loc='lower center', ncol=2, fontsize=10, bbox_to_anchor=(0.5, 0.0))

ax_xlabel = axes[-1]
ax_xlabel.set_xlabel('Experimental level energy (cm⁻¹)', fontsize=10)

plt.tight_layout(rect=[0.08, 0.03, 1, 0.99])
outfile = f'{base}/opt_cluster_1/rohf_casscf_nevpt2_siso/4f_levels_detail_comparison.png'
plt.savefig(outfile, dpi=200, bbox_inches='tight')
print(f"Saved: {outfile}")

# Also print per-manifold detailed table
print(f"\n{'='*80}")
print(f"  Per-Manifold Detailed Comparison (cm⁻¹)")
print(f"{'='*80}")
for m_name in manifold_order:
    mask = np.array([m == m_name for m in lit_manifolds])
    lit_m = lit_levels_all[mask]
    opt_m = np.array([opt[np.argmin(np.abs(opt - v))] for v in lit_m])
    noopt_m = np.array([noopt[np.argmin(np.abs(noopt - v))] for v in lit_m])
    print(f"\n  {manifold_label[m_name]}:")
    print(f"  {'Exp':>8s}  {'opt':>8s}  {'Δopt':>8s}  {'no-opt':>8s}  {'Δno':>8s}")
    for e, o, n in zip(lit_m, opt_m, noopt_m):
        print(f"  {e:8.0f}  {o:8.1f}  {o-e:+8.1f}  {n:8.1f}  {n-e:+8.1f}")
