#!/bin/bash
# PaScaL_TDMA build environment setup for KISTI-6 GPU partition
# Usage: source setup_env.sh

module purge
module load PrgEnv-nvidia
module load craype-accel-nvidia90
module load craype-arm-grace

export MPICH_GPU_SUPPORT_ENABLED=1

module li
