#!/bin/bash

DATASET=$1

module load StdEnv/2023  gcc/12.3  openmpi/4.1.5 raxml-ng/1.2.0 python/3.11.5

RES_DIR=./results/${DATASET}
LOG_DIR=./log/${DATASET}
TREES_FILE_NWK=${RES_DIR}/${DATASET}_best_trees.nwk
SCORES_FILE=${RES_DIR}/${DATASET}_scores.csv
rm -f ${TREES_FILE_NWK} ${SCORES_FILE}
touch ${TREES_FILE_NWK}

NB_TAXA=`grep ">" -c ./data/${DATASET}.fasta`
echo "Dataset size:" ${NB_TAXA}

echo "Reading best trees and scores"
for i in {1..10}
do
    RES_FILE=${RES_DIR}/${DATASET}_${i}.tsv
    LOG_FILE=${LOG_DIR}/${DATASET}_${i}.err
    if [ ! -s "${LOG_FILE}" ]; then
	SCORE=`tail -1 ${RES_FILE} | cut -f 2`
	echo ${DATASET}","${i}","${SCORE} >> ${SCORES_FILE}
	BEST_TREE_FILE_NWK=${RES_DIR}/${DATASET}_best_tree_${i}.nwk
	tail -1 ${RES_FILE} | cut -f 3 >> ${TREES_FILE_NWK}
	tail -1 ${RES_FILE} | cut -f 3 > ${BEST_TREE_FILE_NWK}
    else
	echo "run ${i} terminated with an error"
    fi
done

echo "HOP distance between best trees"
TREES_FILE_TVC=${RES_DIR}/${DATASET}_best_trees.treevec
python ../src/CEDAR.py fromNewick --input_file ${TREES_FILE_NWK} --output_file ${TREES_FILE_TVC}

HOP_DIST_FILE=${RES_DIR}/${DATASET}_HOP.dist
python ../src/CEDAR.py HOP_sim --input_file ${TREES_FILE_TVC} --output_file ${HOP_DIST_FILE} --mode pairwise

echo "RF distance between best trees"

RF_DIST_FILE=${RES_DIR}/${DATASET}_RF.dist
echo "#tree1,tree2,RFdistance" > ${RF_DIST_FILE}
for i in {1..10}
do
    TREE1=${RES_DIR}/${DATASET}_best_tree_${i}.nwk
    for j in {1..10}
    do
	if [ "${j}" -gt "${i}" ]; then
	    TREE2=${RES_DIR}/${DATASET}_best_tree_${j}.nwk
	    RF_DIST=`rf ${TREE1} ${TREE2}`
	    echo ${i}","${j}","${RF_DIST} >> ${RF_DIST_FILE}
	fi
    done
done
