# CEDAR: hill-climbing tree space exploration

This folder contains the scripts and results of experiments on using the hill-climbing heuristic implemented in CEDAR for tree space exploration.
The experiments repeat the analysis of four datasets used in the paper <a href="https://doi.org/10.1093/sysbio/syae030"> Phylo2Vec: A Vector Representation for Binary Trees</a>.

## Data
Data is taken from the datasets used in the Phylo2Vec paper, as given in the daaets folder of the Phylo2Vec github repo, <a href="https://github.com/sbhattlab/phylo2vec/tree/main/py-phylo2vec/phylo2vec/datasets/descr">Phylo2Vec datasets</a>:
- <a href="https://github.com/neherlab/treetime_examples/blob/master/data/zika/zika.fasta">Zika dataset</a> (86 taxa, 10,807nt sequences)  
- <a href="https://github.com/4ment/phylostan/blob/master/examples/fluA/fluA.fa">fluA dataset</a> (69 taxa, 987nt sequences)  
- <a href="https://github.com/neherlab/treetime_examples/blob/master/data/h3n2_na/h3n2_na_20.fasta">h3n2_na dataset</a> (19 taxa, 1,407nt sequences)  
- <a href="https://github.com/KlausVigo/phangorn/blob/main/data/yeast.RData">yeast dataset</a> (8 taxa, 127,018nt sequences)  

```
> date
Thu Oct  2 18:45:45 PDT 2025
data > cd data
data > wget https://raw.githubusercontent.com/neherlab/treetime_examples/refs/heads/master/data/zika/zika.fasta
data > wget https://raw.githubusercontent.com/4ment/phylostan/refs/heads/master/examples/fluA/fluA.fa
data > mv fluA.fa fluA.fasta
data > wget https://raw.githubusercontent.com/neherlab/treetime_examples/refs/heads/master/data/h3n2_na/h3n2_na_20.fasta
data > R
> library(phangorn)
Loading required package: ape
> data(yeast)
> write.phyDat(yeast, "yeast.fasta", format="fasta")
```

## Method: Hill-climbing heuristic

Starting from a random tree, the heuristic iterates the following steps
  - reorder randomly the leaves of the current tree
  - compute the likelihood of all trees in the HOP neighbourhood of the current tree, using RAxML-NG and the GTR model  
  - select the best likelihood tree
  - if its likelihood is within a given tolerance (`0.001`) of the best tree so far:
    - decrease a patience counter [patience step]
  - otherwise:
    - the best tree becomes the current tree
    - the patience counter is reset to max_patience (`max_patience=5`)
    
until the maximum number of iterations is reached or the patience counter is 0

## Experiments

For each dataset, we run 10 hill-climbing explorations with a random starting tree.
```
sed 's/DATASET/zika/g' run_template.sh > run_zika.sh
sed 's/DATASET/fluA/g' run_template.sh > run_fluA.sh
sed 's/DATASET/h3n2_na_20/g' run_template.sh > run_h3n2_na_20.sh
sed 's/DATASET/yeast/g' run_template.sh > run_yeast.sh
experiments_hc > sbatch run_zika.sh
experiments_hc > sbatch run_fluA.sh
experiments_hc > sbatch run_h3n2_na_20.sh
experiments_hc > sbatch run_yeast.sh
```

To summarize the results, we compute the pairwise HOP similarity between the last tree obtained in each run.
For each dataset `DATASET` we compute several files:
- `results/<DATASET>/<DATASET>_scores.csv`: ML score of the last tree of each run;  
- `results/<DATASET>/<DATASET>_best_tree_[1-10].nwk`: last tree of each run in Newick format;  
- `results/<DATASET>/<DATASET>_best_trees.nwk`: last tree of each run in Newick format;  
- `results/<DATASET>/<DATASET>_best_trees.treevec`: last tree of each run in TreeVec format;  
- `results/<DATASET>/<DATASET>_HOP.dist`: pairwise HOP similarity between trees above.
- `results/<DATASET>/<DATASET>_RF.dist`: pairwise normalized Robinson-Foulds distance between trees above.	
```
experiments_hc > ./summary.sh fluA
experiments_hc > ./summary.sh h3n2_na_20
experiments_hc > ./summary.sh zika
experiments_hc > ./summary.sh yeast
```

Looking at the scores files, all runs end on trees with very similar ML scores.
For the small (8 taxa) yeast dataset, all runs end at the same tree.
These resuls are similar to the ones obained by Phylo2Vec with a similar hill-climbing heuristic.
