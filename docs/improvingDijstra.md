# Summary: Robust Track Skeleton Extraction Algorithm

## Purpose

Improve the reliability of track traversal and skeleton extraction for single- and double-electron classification by suppressing diffusion branches and enforcing physics-weighted continuity, while keeping the Li & Zeng endpoint energy method unchanged.

This algorithm addresses failures where geometric Dijkstra traversal:
- jumps onto side branches,
- truncates the physical backbone,
- or misses one of the Bragg peaks.

---

## Processing Pipeline

### 1. Pre-pruning (Recommended)

Remove weak and isolated voxels before graph construction:

- Remove voxel `i` if:
  - `E_i < E_min`
  - AND `deg(i) < 2` after proximity linking

This suppresses:
- delta-electron twigs,
- diffusion noise,
- low-energy dead ends.

---

### 2. Graph Construction with Physics-Weighted Edges

Build a proximity graph:

- Connect voxels `i`, `j` if:
  - `|r_i - r_j| ≤ R_link`

Assign **energy-weighted edge cost**: