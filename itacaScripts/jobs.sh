#!/bin/bash
# ITACA reconstruction jobs
# Generated: 2026-02-23
# 12 jobs: 3 regions × 2 isotopes × 2 configurations

SCRIPT_DIR="/Users/jjgomezcadenas/Projects/Petit/itacaScripts"
DATA_BASE="/Users/jjgomezcadenas/Data/HD5t/precdr"

# Configuration 1: Electron track (σ=10mm)
ELECTRON_OPTS="--voxel-size=10.0 --sigmat=10.0 --sigmal=1.0 --max-distance=10.0 --energy-threshold=0"

# Configuration 2: Ion track (σ=1mm)
ION_OPTS="--voxel-size=2.0 --sigmat=1.0 --sigmal=1.0 --max-distance=2.0 --energy-threshold=0"

# Common options
COMMON_OPTS="--readdir --batch --ievt=1 --levt=5000"

# =============================================================================
# Bi-214 jobs
# =============================================================================

# bi214 copper_shell - electron track
julia "$SCRIPT_DIR/itaca_simple_reco.jl" "$DATA_BASE/copperbkg/bi214" bi214_copper_shell_filtered_ecut_2400_2500 $ELECTRON_OPTS $COMMON_OPTS

# bi214 copper_shell - ion track
julia "$SCRIPT_DIR/itaca_simple_reco.jl" "$DATA_BASE/copperbkg/bi214" bi214_copper_shell_filtered_ecut_2400_2500 $ION_OPTS $COMMON_OPTS

# bi214 ptfe_barrel - electron track
julia "$SCRIPT_DIR/itaca_simple_reco.jl" "$DATA_BASE/innerbkg/bi214" bi214_ptfe_barrel_filtered_ecut_2400_2500 $ELECTRON_OPTS $COMMON_OPTS

# bi214 ptfe_barrel - ion track
julia "$SCRIPT_DIR/itaca_simple_reco.jl" "$DATA_BASE/innerbkg/bi214" bi214_ptfe_barrel_filtered_ecut_2400_2500 $ION_OPTS $COMMON_OPTS

# bi214 ptfe_endcap - electron track
julia "$SCRIPT_DIR/itaca_simple_reco.jl" "$DATA_BASE/innerbkg/bi214" bi214_ptfe_endcap_filtered_ecut_2400_2500 $ELECTRON_OPTS $COMMON_OPTS

# bi214 ptfe_endcap - ion track
julia "$SCRIPT_DIR/itaca_simple_reco.jl" "$DATA_BASE/innerbkg/bi214" bi214_ptfe_endcap_filtered_ecut_2400_2500 $ION_OPTS $COMMON_OPTS

# =============================================================================
# Tl-208 jobs
# =============================================================================

# tl208 copper_shell - electron track
julia "$SCRIPT_DIR/itaca_simple_reco.jl" "$DATA_BASE/copperbkg/tl208" tl208_copper_shell_filtered_ecut_2400_2500 $ELECTRON_OPTS $COMMON_OPTS

# tl208 copper_shell - ion track
julia "$SCRIPT_DIR/itaca_simple_reco.jl" "$DATA_BASE/copperbkg/tl208" tl208_copper_shell_filtered_ecut_2400_2500 $ION_OPTS $COMMON_OPTS

# tl208 ptfe_barrel - electron track
julia "$SCRIPT_DIR/itaca_simple_reco.jl" "$DATA_BASE/innerbkg/tl208" tl208_ptfe_barrel_filtered_ecut_2400_2500 $ELECTRON_OPTS $COMMON_OPTS

# tl208 ptfe_barrel - ion track
julia "$SCRIPT_DIR/itaca_simple_reco.jl" "$DATA_BASE/innerbkg/tl208" tl208_ptfe_barrel_filtered_ecut_2400_2500 $ION_OPTS $COMMON_OPTS

# tl208 ptfe_endcap - electron track
julia "$SCRIPT_DIR/itaca_simple_reco.jl" "$DATA_BASE/innerbkg/tl208" tl208_ptfe_endcap_filtered_ecut_2400_2500 $ELECTRON_OPTS $COMMON_OPTS

# tl208 ptfe_endcap - ion track
julia "$SCRIPT_DIR/itaca_simple_reco.jl" "$DATA_BASE/innerbkg/tl208" tl208_ptfe_endcap_filtered_ecut_2400_2500 $ION_OPTS $COMMON_OPTS
