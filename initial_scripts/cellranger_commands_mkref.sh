#!/bin/bash
#SBATCH --qos=normal            # Quality of Service
#SBATCH --job-name=NICD_FC      # Job Name
#SBATCH --nodes=4               # Number of Nodes
#SBATCH --ntasks-per-node=2     # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=8       # Number of threads per task (OMP threads)

cd /lustre/project/wdeng7/deep/CellRangerData/

cellranger mkref --genome=refgenome4cellranger --fasta=/lustre/project/wdeng7/deep/dm6/genomes/Drosophila_melanogaster.BDGP6.28.dna.toplevel_GAL4_EGFP.fa --genes=/lustre/project/wdeng7/deep/dm6/genes/Drosophila_melanogaster.BDGP6.28.100.chr_filtered.GAL4.EGFP.gtf
