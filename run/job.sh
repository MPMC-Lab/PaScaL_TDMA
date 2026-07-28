#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --partition=gpu
#SBATCH --hint=nomultithread
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH -J PaScaLTDMA
#SBATCH -o log_%j.out

if [[ -n "${MODULESHOME:-}" && -f "${MODULESHOME}/init/bash" ]]; then
    source "${MODULESHOME}/init/bash"
fi
if ! command -v module >/dev/null 2>&1; then
    echo "error: environment-modules command is unavailable" >&2
    exit 1
fi
module purge
module load PrgEnv-nvidia
module load craype-accel-nvidia90

export MPICH_GPU_SUPPORT_ENABLED=1

cd "$SLURM_SUBMIT_DIR"

srun -N 1 -n 4 "$SLURM_SUBMIT_DIR/a.out"
