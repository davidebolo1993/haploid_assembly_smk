#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate snakemakeenv_latest
snakemake all --profile config/slurm --rerun-triggers mtime

