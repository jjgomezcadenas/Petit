# Petit HDF5 Batch Processing Script

## Overview
`run_hd5t.jl` is a batch processing script for analyzing HDF5 files using the Petit package. It processes particle physics event data, performs track reconstruction, and generates histograms for different track categories.

## Usage

### Basic Usage
```bash
julia run_hd5t.jl <input_directory>
```

### With Options
```bash
julia run_hd5t.jl <input_directory> [options]
```

### Available Options
- `--input-file=<file>` - Input HDF5 file name (default: `0nubb.next.h5`)
- `--events=<n>` - Number of events to process (default: 100)
- `--voxel-size=<mm>` - Voxel size in millimeters (default: 5)
- `--max-distance=<mm>` - Maximum distance for track building in mm (default: 10)
- `--energy-threshold=<keV>` - Energy threshold in keV (default: 10)
- `--output-dir=<dir>` - Output directory for results (default: `znubb`)

### Examples

Process with default parameters:
```bash
julia run_hd5t.jl /path/to/data/
```

Process 1000 events with custom parameters:
```bash
julia run_hd5t.jl /path/to/data/ --events=1000 --voxel-size=3 --output-dir=results
```

Process a specific file:
```bash
julia run_hd5t.jl /path/to/data/ --input-file=mydata.h5 --events=500
```

## Output Files

The script generates the following output files in the specified output directory:

1. **HSt1.txt** - Histograms for single track events
2. **HSt2p.txt** - Histograms for two-track events (primary track)
3. **HSt2s.txt** - Histograms for two-track events (secondary track)
4. **HSt3p.txt** - Histograms for three+ track events (primary track)
5. **HSt3s.txt** - Histograms for three+ track events (secondary tracks)
6. **analysis_summary.txt** - Summary statistics of the analysis run

Each histogram file contains:
- `hx` - X-coordinate distribution
- `hy` - Y-coordinate distribution
- `hz` - Z-coordinate distribution
- `he` - Energy distribution

## Running as a Batch Job

### On a local machine:
```bash
nohup julia run_hd5t.jl /data/input/ --events=10000 --output-dir=results_batch &
```

### On a cluster with SLURM:
Create a submission script `submit_petit.sh`:
```bash
#!/bin/bash
#SBATCH --job-name=petit_analysis
#SBATCH --output=petit_%j.out
#SBATCH --error=petit_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G

module load julia/1.11  # Adjust to your system
julia run_hd5t.jl /path/to/data/ --events=10000 --output-dir=results_${SLURM_JOB_ID}
```

Submit with:
```bash
sbatch submit_petit.sh
```

## Requirements
- Julia 1.x or higher
- Petit package and its dependencies
- Input HDF5 files in the expected format

## Error Handling
The script includes error handling for:
- Missing input directories or files
- Invalid command line arguments
- Processing errors during analysis

Error messages will indicate the specific issue and suggest solutions.