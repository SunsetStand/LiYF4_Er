#!/usr/bin/env python3
"""
Build merged 4-layer embedded cluster:
- Inner region (r < d/2): optimized CONTCAR coordinates
- Outer region (r > d/2): pre-optimization POSCAR coordinates
- cluster.xyz: ErF8 only (r < rCluster)
- Re-fit surface charges using SVD
"""

import sys, os
import numpy as np

# Add path for absolute imports (same as main.py)
sys.path.insert(0, '/data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main')
sys.path.insert(0, '/data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/src')
from fitting_ewald.neighborTools import neighbors
from fitting_ewald.exact_potential import formal_charges
from fitting_ewald.potential_fitting import PotentialFitOnlyCharges
from ase.io import read

# === Parameters ===
d_half = 5.1361       # half of minimum Er-Er distance in optimized supercell
rCluster = 2.993       # for ErF8 (8 nearest F atoms at ~2.2-2.3A, next F at ~3.7A)
rAIMP = 8.0
rChgs = 12.0
rSurface = 14.0
num_sites = 800

# File paths
pre_opt_poscar = "/data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+/no_opt_2/LiYF4 (1).poscar"
post_opt_poscar = "/data/home/wangcx/LiYF4_Er3+/optv2/CONTCAR"
workdir = "/data/home/wangcx/LiYF4_Er3+/AIMPModelGenerator-main/LiYF4_Er3+/opt_cluster_2"

os.makedirs(workdir, exist_ok=True)
os.chdir(workdir)

print("=" * 60)
print("  Building Merged 4-Layer Embedded Cluster")
print("=" * 60)
print(f"  d/2 = {d_half:.4f} A")
print(f"  rCluster (ErF8) = {rCluster:.3f} A")
print(f"  rAIMP = {rAIMP} A, rChgs = {rChgs} A, rSurface = {rSurface} A")
print()

# ============================================================
# Step 1: Build neighbors from both structures
# ============================================================
print("--- Step 1: Building neighbor lists ---")

# Pre-opt: undoped primitive cell, centered on first Y atom
pre_nbs = neighbors(pre_opt_poscar, cAtom='Y', cAtomIndex=1, rCut=rSurface, sort=True)
n_pre = pre_nbs.get_number_of_neighbors()
print(f"  Pre-opt  (LiYF4 prim, Y-center): {n_pre} atoms within {rSurface} A")

# Post-opt: optimized 2x2x2 supercell with Er, centered on Er
post_nbs = neighbors(post_opt_poscar, cAtom='Er', cAtomIndex=1, rCut=max(rSurface, d_half), sort=True)
n_post = post_nbs.get_number_of_neighbors()
print(f"  Post-opt (CONTCAR, Er-center): {n_post} atoms within {max(rSurface, d_half)} A")

# ============================================================
# Step 2: Extract atom data
# ============================================================
print("\n--- Step 2: Extracting atom data ---")

def get_atom_data(nbs):
    indices, offsets = nbs.get_neighbors()
    distances = nbs.get_distances()
    coords = nbs.get_cartesian_coordinates(origin_shifted=True)
    symbols = np.array([nbs.mol.get_chemical_symbols()[idx] for idx in indices])
    return symbols, coords, distances

post_syms, post_coords, post_dists = get_atom_data(post_nbs)
pre_syms, pre_coords, pre_dists = get_atom_data(pre_nbs)

# ============================================================
# Step 3: Merge atoms
# ============================================================
print("\n--- Step 3: Merging at d/2 boundary ---")

# Strategy:
# 1. All post-opt atoms within d/2 go into merged list
# 2. Pre-opt atoms beyond d/2 go into merged list
# 3. Skip pre-opt atoms that are very close to any post-opt atom (boundary overlap)

merged_syms = []
merged_coords = []
merged_dists = []
merged_layer = []  # 'cluster', 'aimp', 'rawChgs', 'surfChgs'

# Post-opt inner atoms
post_coords_inner = post_coords[post_dists <= d_half]
post_inner_count = 0
for s, c, d in zip(post_syms, post_coords, post_dists):
    if d <= d_half:
        if d <= rCluster:
            layer = 'cluster'
        elif d <= rAIMP:
            layer = 'aimp'
        else:
            layer = 'rawChgs'
        merged_syms.append(s)
        merged_coords.append(c)
        merged_dists.append(d)
        merged_layer.append(layer)
        post_inner_count += 1

print(f"  Post-opt atoms within d/2 (={d_half:.2f}A): {post_inner_count}")

# Pre-opt outer atoms (beyond d/2)
pre_outer_count = 0
overlap_skipped = 0
for s, c, d in zip(pre_syms, pre_coords, pre_dists):
    if d > d_half:
        # Check overlap with post-opt inner atoms
        min_dist_to_post = np.min(np.linalg.norm(post_coords_inner - c, axis=1)) if len(post_coords_inner) > 0 else 1e10
        if min_dist_to_post < 0.3:  # 0.3A tolerance
            overlap_skipped += 1
            continue
        
        if d <= rAIMP:
            layer = 'aimp'
        elif d <= rChgs:
            layer = 'rawChgs'
        elif d <= rSurface:
            layer = 'surfChgs'
        else:
            continue
        
        merged_syms.append(s)
        merged_coords.append(c)
        merged_dists.append(d)
        merged_layer.append(layer)
        pre_outer_count += 1

print(f"  Pre-opt atoms beyond d/2: {pre_outer_count}")
print(f"  Overlap atoms skipped: {overlap_skipped}")

merged_syms = np.array(merged_syms)
merged_coords = np.array(merged_coords)
merged_dists = np.array(merged_dists)
merged_layer = np.array(merged_layer)

# Summary by layer
for layer_name in ['cluster', 'aimp', 'rawChgs', 'surfChgs']:
    mask = merged_layer == layer_name
    count = np.sum(mask)
    if count > 0:
        uniq, cnts = np.unique(merged_syms[mask], return_counts=True)
        comp = ", ".join([f"{u}:{c}" for u, c in zip(uniq, cnts)])
        print(f"  {layer_name}: {count} atoms ({comp})")

# ============================================================
# Step 4: Write output files (preliminary, before re-fitting)
# ============================================================
print("\n--- Step 4: Writing output files ---")

def write_xyz(filename, symbols, coords):
    with open(filename, 'w') as f:
        f.write(f"  {len(symbols)}\n")
        f.write(f" structure file generated manually.\n")
        for s, c in zip(symbols, coords):
            f.write(f"  {s:3s}  {c[0]:15.7f} {c[1]:15.7f} {c[2]:15.7f}\n")

for layer_name in ['cluster', 'aimp', 'rawChgs', 'surfChgs']:
    mask = merged_layer == layer_name
    fname = f"{layer_name}.xyz"
    write_xyz(fname, merged_syms[mask], merged_coords[mask])
    print(f"  Wrote {fname}: {np.sum(mask)} atoms")

# rawCharges.dat: formal charges for rawChgs layer
raw_mask = merged_layer == 'rawChgs'
raw_charges = np.array([formal_charges[s] for s in merged_syms[raw_mask]])
np.savetxt('rawCharges.dat', raw_charges)
print(f"  Wrote rawCharges.dat: {len(raw_charges)} values")

# temp surface charges (will be replaced by fitted ones)
surf_mask = merged_layer == 'surfChgs'
surf_raw = np.array([formal_charges[s] for s in merged_syms[surf_mask]])
np.savetxt('surfaceCharges.dat', surf_raw)

# ============================================================
# Step 5: Generate pre-opt surface charges via PotentialFitOnlyCharges
#   The surface charges come from the pre-opt structure (same lattice as
#   post-opt, ISIF=2). Inner-region relaxation has negligible effect on
#   far-field Madelung potential, so pre-opt surface charges are sufficient.
#   This way, adjusting cluster radii auto-updates the fitting.
# ============================================================
print("\n--- Step 5: Fitting surface charges on pre-opt structure ---")
print("  Using PotentialFitOnlyCharges on pre-opt crystal...")

pot_fit = PotentialFitOnlyCharges(
    pre_opt_poscar,
    cAtom='Y', cAtomIndex=1,
    rCut=rChgs, rSurface=rSurface,
    num_sites=num_sites
)
pot_fit.run_fit()
pot_fit.show_res()

# Write pre-opt surface charge .xyz and fitted .dat
pre_surf_nbs = pot_fit.surface_neighbors
pre_surf_coords = pre_surf_nbs.get_cartesian_coordinates(origin_shifted=True)
pre_surf_syms = [pre_surf_nbs.mol.get_chemical_symbols()[idx] 
                 for idx in pre_surf_nbs.get_neighbors()[0]]
write_xyz('surfChgs.xyz', pre_surf_syms, pre_surf_coords)
np.savetxt('surfaceCharges.dat', pot_fit.surf_chgs)

print(f"\n  Wrote surfChgs.xyz: {len(pre_surf_syms)} atoms")
print(f"  Wrote surfaceCharges.dat: {len(pot_fit.surf_chgs)} fitted charges")

# Copy CONTCAR for reference
import shutil
shutil.copy(post_opt_poscar, 'CONTCAR')

print("\n" + "=" * 60)
print("  DONE! All files written to:")
print(f"  {workdir}")
print("  Surface charges: fitted on pre-opt crystal via PotentialFitOnlyCharges")
print("=" * 60)
print("\n  Output files:")
for f in ['cluster.xyz', 'aimp.xyz', 'rawChgs.xyz', 'rawCharges.dat',
           'surfChgs.xyz', 'surfaceCharges.dat', 'CONTCAR']:
    fpath = os.path.join(workdir, f)
    if os.path.exists(fpath):
        size = os.path.getsize(fpath)
        if f.endswith('.xyz'):
            with open(fpath) as fh:
                n_atoms = fh.readline().strip()
            print(f"  {f}: {n_atoms} atoms, {size} bytes")
        else:
            print(f"  {f}: {size} bytes")
