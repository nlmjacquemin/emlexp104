
# Supplementary Scripts for Reproducible Analysis

Scripts used for the analysis presented in the paper: 

These scripts were executed on the EPFL high-performance computing (HPC) cluster using SLURM (v23.11.10) and Apptainer (v1.2.5). Each node consisted of two Intel(R) Xeon(R) Platinum 8360Y processors running at 2.4 GHz, with 36 cores per processor (72 cores per node), and 3 TB of SSD storage.

## How to Use

### Requirements

Ensure the following are installed and available:
- **SLURM**
- **Apptainer**
- Required container images (**Docker**)

If the dependencies are not available, you can execute each script manually (e.g., `bash run_coverm.sh [arguments]`) with the appropriate container image or install the necessary software on your system.

### General Usage

Scripts can be launched using the provided launcher and Apptainer scripts. Parameters are specified for each script using the `--` separator, as shown below:

```bash
sbatch \
  --cpus-per-task="CPUs (e.g., 72)" \
  --time="HH:MM:SS (e.g., 96:00:00)" \
  launcher.sh \
  -s apptainer \
  --working_directory="path to working directory (e.g., /scratch/username/projects/E80)" \
  -- \
  -i "image name without extension (e.g., coverm)" \
  -s "script name without extension (e.g., run_coverm)" \
  --micromamba_env="conda environment name (e.g., coverm)"
```

#### Adding Arguments for Apptainer

You can specify additional arguments for Apptainer commands as follows:

```bash
sbatch \
  --cpus-per-task="CPUs (e.g., 72)" \
  --time="HH:MM:SS (e.g., 96:00:00)" \
  launcher.sh \
  -s apptainer \
  --working_directory="path to working directory (e.g., /scratch/username/projects/E80)" \
  --bind="bind directory (e.g., /scratch/username/gtdb:/db)" \
  -- \
  -i "image name without extension (e.g., coverm)" \
  -s "script name without extension (e.g., run_coverm)" \
  --micromamba_env="conda environment name (e.g., coverm)"
```

#### Changing Input Parameters for the Script

To modify input parameters for the script, use the following pattern:

```bash
sbatch \
  ... \
  --micromamba_env="conda environment name (e.g., coverm)" \
  -- \
  --output_folder="output directory path (e.g., analysis/coverm_bis)"
```

#### Passing Additional Arguments to the Script Command

You can directly pass additional arguments to the command executed within the script, like this:

```bash
sbatch \
  ... \
  --micromamba_env="conda environment name (e.g., coverm)" \
  -- \
  -- \
  --methods="method(s) for calculating coverage (e.g., mean)"
```

---

### Key Notes:
1. Use `--` to separate script parameters from the SLURM/launcher parameters.
2. Use an additional `--` to separate script-level arguments from the command-line arguments executed within the script.


