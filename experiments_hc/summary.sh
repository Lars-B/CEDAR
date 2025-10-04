#!/bin/bash

DATASET=$1
LOG_ID=$2

module load StdEnv/2023  gcc/12.3  openmpi/4.1.5 raxml-ng/1.2.0 python/3.11.5

RES_DIR=./results/${DATASET}
LOG_DIR=./log/${DATASET}
TREES_FILE_NWK=${RES_DIR}/${DATASET}_best_trees.nwk
rm -f ${TREES_FILE_NWK}
touch ${TREES_FILE_NWK}

for i in {1..10}
do
    RES_FILE=${RES_DIR}/${DATASET}_${i}.csv
    LOG_FILE=${LOG_DIR}/${DATASET}_${LOG_ID}_${i}.err
    if [ ! -s "${LOG_FILE}" ]; then
	tail -1 ${RES_FILE} | cut -f 3 >> ${TREES_FILE_NWK}
    else
	echo "run ${i} terminated with an error"
    fi
done

TREES_FILE_TVC=${RES_DIR}/${DATASET}_best_trees.treevec
python ../src/CEDAR.py fromNewick --input_file ${TREES_FILE_NWK} --output_file ${TREES_FILE_TVC}

DIST_FILE=${RES_DIR}/${DATASET}.dist
python ../src/CEDAR.py HOP_sim --input_file ${TREES_FILE_TVC} --output_file ${DIST_FILE} --mode pairwise
