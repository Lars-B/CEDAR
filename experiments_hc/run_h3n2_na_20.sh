#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --mem=8G
#SBATCH --account=def-chauvec
#SBATCH --job-name=h3n2_na_20
#SBATCH --array=1-10
#SBATCH --output=/scratch/chauvec/CEDAR/CEDAR/experiments_hc/log/h3n2_na_20/h3n2_na_20_%a.out
#SBATCH --error=/scratch/chauvec/CEDAR/CEDAR/experiments_hc/log/h3n2_na_20/h3n2_na_20_%a.err

module load StdEnv/2023  gcc/12.3  openmpi/4.1.5 raxml-ng/1.2.0 python/3.11.5

CEDAR_DIR=/scratch/chauvec/CEDAR/CEDAR
EXP_DIR=${CEDAR_DIR}/experiments_hc

RES_OUT=${EXP_DIR}/results/h3n2_na_20/h3n2_na_20_${SLURM_ARRAY_TASK_ID}.tsv
TREE_DIR=${EXP_DIR}/tree_folder/h3n2_na_20/run_${SLURM_ARRAY_TASK_ID}

rm -rf  ${RES_OUT} ${TREE_DIR}
mkdir -p ${TREE_DIR}

python ${CEDAR_DIR}/src/CEDAR.py HOP_hc \
       --fasta_path ${EXP_DIR}/data/h3n2_na_20.fasta \
       --DNA_model "GTR" \
       --tree_folder_path ${TREE_DIR} \
       --out_file_path ${RES_OUT} \
       --seed ${SLURM_ARRAY_TASK_ID}
