## Smk pipeline for haploid assembly

```bash
snakemake all --profile config/slurm --rerun-triggers mtime 
```