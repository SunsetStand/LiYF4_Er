 import numpy as np
from ase.io import read
import os

# Read the optimized CONTCAR
contcar = read("/data/home/wangcx/LiYF4_Er3+/optv2/CONTCAR")
print(f"CONTCAR: {len(contcar)} atoms, cell: {contcar.get_cell()}")

# Get the Er atom in the center — shift without wrapping
er_idx = [i for i, s in enumerate(contcar.get_chemical_symbols()) if s == "Er"]
er_pos_orig = contcar.get_positions()[er_idx[0]]
print(f"Er original pos: {er_pos_orig}")

# Shift all atoms relative to Er (NO wrapping)
all_pos = contcar.get_positions() - er_pos_orig
symbols = contcar.get_chemical_symbols()

# For atoms that are far, apply fractional wrapping to get nearest image
cell = contcar.get_cell()
for i in range(len(all_pos)):
    # Convert to fractional, wrap to [-0.5, 0.5)
    frac = np.linalg.solve(cell.T, all_pos[i].T).T
    frac = frac - np.round(frac)
    all_pos[i] = cell.T @ frac

distances = np.linalg.norm(all_pos, axis=1)

# Count atoms by distance
for r in [2.5, 3.0, 3.5, 4.0, 5.0, 8.0]:
    n = np.sum(distances < r)
    print(f"  r < {r:.1f}: {n} atoms")

# Use rCluster = 3.8 to match no_opt_2 pattern (same 25-atom selection)
rCluster = 3.8
rAIMP = 8.0

cluster_mask = distances < rCluster
aimp_mask = (distances >= rCluster) & (distances < rAIMP)

n_cluster = np.sum(cluster_mask)
n_aimp = np.sum(aimp_mask)
print(f"\nCluster: {n_cluster} atoms (r < {rCluster})")
print(f"AIMP: {n_aimp} atoms ({rCluster} < r < {rAIMP})")

# Sort cluster atoms by distance
c_idx = np.where(cluster_mask)[0]
c_dist = distances[c_idx]
c_sort = np.argsort(c_dist)
c_idx = c_idx[c_sort]

outdir = "/data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+/opt_cluster_2"
os.makedirs(outdir, exist_ok=True)

# Write cluster.xyz
with open(f"{outdir}/cluster.xyz", "w") as f:
    f.write(f"{n_cluster}\n")
    f.write(" Optimized big cluster (r < 3.8 A) from CONTCAR\n")
    for idx in c_idx:
        sym = symbols[idx]
        pos = all_pos[idx]
        f.write(f"  {sym:3s} {pos[0]:14.7f} {pos[1]:14.7f} {pos[2]:14.7f}\n")
print(f"\nWrote {outdir}/cluster.xyz ({n_cluster} atoms)")

# Write aimp.xyz
with open(f"{outdir}/aimp.xyz", "w") as f:
    a_idx = np.where(aimp_mask)[0]
    a_dist = distances[a_idx]
    a_sort = np.argsort(a_dist)
    a_idx = a_idx[a_sort]
    f.write(f"{len(a_idx)}\n")
    f.write(" AIMP environment from optimized CONTCAR\n")
    for idx in a_idx:
        sym = symbols[idx]
        pos = all_pos[idx]
        f.write(f"  {sym:3s} {pos[0]:14.7f} {pos[1]:14.7f} {pos[2]:14.7f}\n")
print(f"Wrote {outdir}/aimp.xyz ({len(a_idx)} atoms)")

# Show cluster composition
print("\nCluster composition:")
from collections import Counter
c_syms = [symbols[i] for i in c_idx]
for sym, count in Counter(c_syms).items():
    print(f"  {sym}: {count}")
