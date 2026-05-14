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
from fitting_ewald.exact_potential import get_exact_potential, formal_charges
from fitting_ewald.finite_potential import get_finite_potential, coef_pot
from fitting_ewald.genSites import gen_random_sites
from ase.io import read
from ase import Atoms

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
# Step 5: Re-fit surface charges
# ============================================================
print("\n--- Step 5: Re-fitting surface charges ---")
print("  (Computing Ewald potential via calcmad...)")

# Surface charge atoms (from pre-opt, since rChgs > d/2)
surf_indices = np.where(surf_mask)[0]
surf_coords_only = merged_coords[surf_indices]
surf_syms_only = merged_syms[surf_indices]
n_surf = len(surf_indices)
print(f"  Number of surface charges to fit: {n_surf}")

# Strategy:
# 1. Use pre-opt poscar for gen_random_sites + get_exact_potential (infinite lattice)
# 2. Use get_finite_potential from pre-opt poscar (properly handles periodic images)
# 3. ADD correction: (post-opt inner contribution) - (pre-opt inner contribution)
#    This accounts for the relaxed inner region
# 4. Build A matrix from merged surface charge coordinates
# 5. SVD fit

pre_mol = read(pre_opt_poscar)
center_index = pre_nbs.center_index
n_sites = max(num_sites, 2 * n_surf)

# Generate random evaluation sites
sites = gen_random_sites(pre_opt_poscar, num=n_sites, deps=0.1, atom=center_index)
n_sites_actual = len(sites)
print(f"  Number of evaluation sites: {n_sites_actual}")

# Step 5a: Exact potential from infinite pre-opt crystal
print("  Computing exact Madelung potential via calcmad...")
sys.stdout.flush()
exact_pot = get_exact_potential(pre_opt_poscar, sites)
print(f"  Exact potential computed.")

# Step 5b: Finite potential from pre-opt crystal (proper periodic treatment)
print("  Computing finite potential from pre-opt crystal (periodic)...")
sys.stdout.flush()
finite_pot_pre = get_finite_potential(pre_opt_poscar, sites, rCut=rChgs, atom=center_index, 
                                        cAtom='Y', cAtomIndex=1)

# Step 5c: Compute correction from inner region relaxation
# Correction = sum_{post inner} q/r - sum_{pre inner} q/r
# Inner atoms within d/2 from the center, in the CENTRAL replica only
print("  Computing inner-region correction...")
sys.stdout.flush()

# Get pre-opt atoms within d/2 (central replica only)
pre_inner_mask = pre_dists <= d_half
pre_inner_syms = pre_syms[pre_inner_mask]
pre_inner_coords = pre_coords[pre_inner_mask]

# Get post-opt atoms within d/2 (central replica only)
post_inner_mask = post_dists <= d_half
post_inner_syms = post_syms[post_inner_mask]
post_inner_coords = post_coords[post_inner_mask]

correction = np.zeros(n_sites_actual)
for i, site in enumerate(sites):
    # Post contribution
    for sym, coord in zip(post_inner_syms, post_inner_coords):
        d = np.linalg.norm(coord - site)
        if d > 1e-10:
            correction[i] += formal_charges[sym] / d
    # Subtract pre contribution
    for sym, coord in zip(pre_inner_syms, pre_inner_coords):
        d = np.linalg.norm(coord - site)
        if d > 1e-10:
            correction[i] -= formal_charges[sym] / d
correction *= coef_pot

# Apply correction: finite_pot = finite_pot_pre + correction
finite_pot = finite_pot_pre + correction
print(f"  Correction applied (inner region relaxation).")
print(f"  Max correction: {np.max(np.abs(correction)):.5f} V")

# Step 5d: Build A matrix from merged surface charge coordinates
print("  Building fitting matrix and solving SVD...")
sys.stdout.flush()

# Convert sites from fractional to Cartesian for A matrix computation
site_coords = np.zeros((n_sites_actual, 3))
pre_cell = pre_mol.get_cell()
center_frac = pre_mol.get_scaled_positions()[center_index]
for i, site in enumerate(sites):
    site_coords[i] = np.dot(site - center_frac, pre_cell)

A = np.zeros((n_sites_actual, n_surf))
for i in range(n_sites_actual):
    for j in range(n_surf):
        d = np.linalg.norm(site_coords[i] - surf_coords_only[j])
        if d > 1e-10:
            A[i, j] = 1.0 / d
A *= coef_pot

b = exact_pot - finite_pot

# SVD fitting
u, s, vt = np.linalg.svd(A, full_matrices=False)
v = vt.T

sigma_threshold = 1e-8
x = np.zeros(n_surf)
n_sigma = 0
for i in range(n_surf):
    if s[i] > sigma_threshold:
        n_sigma += 1
        x += np.dot(u[:, i], b) / s[i] * v[:, i]
    else:
        break

# Compute residuals
residual = np.dot(A, x) - b
raw_chgs_surf = np.array([formal_charges[s] for s in surf_syms_only])
residual_raw = np.dot(A, raw_chgs_surf) - b

mae_fitted = np.abs(residual).sum() / n_sites_actual
mae_raw = np.abs(residual_raw).sum() / n_sites_actual

print(f"\n  === Fitting Results ===")
print(f"  Number of sites = {n_sites_actual}")
print(f"  Number of surface point charges = {n_surf}")
print(f"  Number of singular values used = {n_sigma}")
print(f"  MAE (raw formal charges) = {mae_raw:.5f} V")
print(f"  MAE (fitted charges)     = {mae_fitted:.5f} V")

np.savetxt('surfaceCharges.dat', x)
print(f"\n  Wrote fitted surfaceCharges.dat: {n_surf} values")

# Write a copy of CONTCAR for reference
import shutil
shutil.copy(post_opt_poscar, 'CONTCAR')

print("\n" + "=" * 60)
print("  DONE! All files written to:")
print(f"  {workdir}")
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
