# kNN-Based Track Extreme Finding: Implementation Plan

## Problem Statement

### Current Behavior
Extreme finding fails for U-shaped tracks (e.g. xe137 event 12, 2400–2500 keV):
- **Expected:** extremes at track start (low E) and Bragg peak (high E)
- **Observed:** both extremes found near Bragg region

### Root Cause
`build_tracks_kdtree` currently builds a **radius graph**: connect every voxel to **all** neighbors within `max_distance`.

This produces **shortcut edges** across a U-bend when the two arms come within `R_link`, collapsing topological distances and corrupting diameter / farthest-point logic.

---

## Proposed Solution: k-Nearest Neighbor Graph (kNN)

### Why kNN Works
Replace “all neighbors within radius” by “k closest neighbors”:

- Each node has only **local connectivity** (bounded degree)
- Shortcuts across a U-bend are unlikely to appear because they would need to be among the **k closest** connections for both endpoints
- Graph becomes sparse, stable, and much cheaper to process than dense radius graphs

This attacks the root cause directly: **you prevent shortcut edges from being created**, instead of pruning them afterward (MST).

---

## Algorithm Overview

### 1) Build kNN Graph
For each voxel `i`:
- query `k+1` nearest neighbors (includes self)
- add undirected edges `(i, j)` for the `k` nearest `j ≠ i`

Optional safety:
- include a distance gate `d(i,j) ≤ r_max` to avoid forcing edges in very sparse regions

### 2) Connected Components → Tracks
Same as now: connected components define track candidates.

### 3) Extreme Finding
Keep your existing method:
- farthest-from-arbitrary → farthest-from-farthest
- shortest path between extremes as skeleton

kNN will make this robust on U-turns by removing shortcut edges upstream.

---

## Implementation Plan (Julia)

### Dependencies
- Use `NearestNeighbors.jl` (already in your workflow for KDTree / inrange)
- Keep `Graphs.jl` for `SimpleGraph`, components, paths

---

## Phase 1: Replace Radius Graph with kNN Graph in `build_tracks_kdtree`

### Existing (problematic)
- `idxs = inrange(tree, coords[:, i], r, true)` → adds **all** neighbors in radius

### New (kNN)
- `idxs, dists = knn(tree, coords[:, i], k+1, true)` → add only **k** neighbors

#### Recommended parameters
- `k = 8` to `12` (start with 10)
- Optional `r_max = max_distance` (same scale you use today) as a guardrail

---

## Phase 2: Optional “Mutual kNN” (More Conservative)
To suppress asymmetric edges:
- keep edge `(i,j)` only if `i ∈ kNN(j)` AND `j ∈ kNN(i)`

This further reduces spurious cross-arm links in weird geometries.

Use this only if plain kNN still leaves shortcuts.

---

## Phase 3: Optional Energy-Weighted Traversal (Only If Needed)
If you still see mis-traversals (branch capture, endpoint drift):
- keep kNN graph topology
- change shortest-path edge weight to prefer high-support regions:

`w_ij = d_ij / (0.5*(E_i + E_j) + ε)`

This biases paths toward the physical backbone without enforcing a tree.

---

## Integration Options

### Option A (Recommended): New builder `build_tracks_knn`
- Add a new function alongside the current one
- Switch per-run with a keyword flag (safe, no regressions)

### Option B: Replace radius-graph logic in `build_tracks_kdtree`
- Cleaner long-term
- Riskier (affects all downstream behavior immediately)

Start with Option A, validate on known failure events, then decide.

---

## Testing Plan

### Unit Tests
1. **Linear track:** kNN graph stays connected; extremes unchanged
2. **U-shaped synthetic:** kNN graph should not contain cross-arm shortcuts
3. **Dense blob:** kNN graph remains connected but does not explode in edges
4. **Sparse tail:** verify `r_max` (if used) prevents forced long edges

### Integration Tests
1. **xe137 event 12:** extremes should shift to start ↔ Bragg end
2. **bb0nu sample:** extremes should remain Bragg ↔ Bragg with no regression
3. **Performance:** ensure wall-time improves vs radius graph for large tracks

---

## Success Criteria
- Fix U-shaped failures (xe137 event 12) without introducing new ones
- Maintain bb0nu behavior
- Reduce sensitivity to `R_link` tuning
- Keep runtime stable or improved (kNN graph is sparse by construction)

---

## Files to Modify

| File | Changes |
|------|---------|
| `src/tracks.jl` | Add `build_tracks_knn` (or replace radius edge-building block) |
| `src/Petit.jl`  | Export new function (optional) |
| `Project.toml`  | Add `NearestNeighbors` only if not already present |

---

## Complexity Notes

| Approach | Edges | Cost |
|----------|-------|------|
| Radius graph | variable, can be dense | expensive in dense regions; creates shortcuts |
| kNN graph | ~ `N*k/2` | predictable, sparse, stable |
| MST on dense graph | extra `O(E log E)` | can be heavy; not guaranteed to remove shortest shortcuts |

Bottom line: **kNN is the simplest upstream fix**: it prevents shortcut edges rather than trying to delete them later.