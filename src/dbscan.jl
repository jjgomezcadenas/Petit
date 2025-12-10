using Clustering
using LinearAlgebra

"""
    dbscan_ion_cloud(x, y, z, q; eps, minpts)

Run DBSCAN on an ion cloud given by arrays x, y, z (positions) and q (charges).

Arguments
---------
- `x, y, z` : vectors (same length N) with coordinates of ions.
- `q`       : vector of charges (same length N). Not used by DBSCAN itself,
              but kept so you can aggregate charges later.
- `eps`     : neighborhood radius for DBSCAN (e.g. ~ diffusion sigma, in mm).
- `minpts`  : minimum number of points to start a cluster (e.g. 10–50).

Returns
-------
- `labels::Vector{Int}` : cluster label for each point (1,2,...,K), or 0 for noise.
- `K::Int`              : number of clusters found (excluding noise).
"""
function dbscan_ion_cloud(x::AbstractVector,
                          y::AbstractVector,
                          z::AbstractVector,
                          q::AbstractVector;
                          eps::Real,
                          minpts::Integer)

    N = length(x)
    @assert length(y) == N "x,y,z,q must have same length"
    @assert length(z) == N "x,y,z,q must have same length"
    @assert length(q) == N "x,y,z,q must have same length"

    # Build a 3×N matrix of points as required by Clustering.dbscan
    X = hcat(x, y, z)'  # size: 3 × N

    # Run DBSCAN (Euclidean metric by default)
    result = dbscan(X, eps, minpts; dist_metric = Euclidean())

    # result.assignments is a Vector{Int} of length N:
    #   0 = noise, 1..K = cluster id
    labels = result.assignments
    K = maximum(labels)  # number of clusters (noise is 0)

    return labels, K
end