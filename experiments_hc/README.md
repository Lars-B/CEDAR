# CEDAR: tree space exploration experiments

## Data
Data is taken from the datasets used in the Phylo2Vec paper, as given in their github repo:
<a href="https://github.com/sbhattlab/phylo2vec/tree/main/py-phylo2vec/phylo2vec/datasets/descr">Phylo2Vec datasets</a>:
- <a href="https://github.com/neherlab/treetime_examples/blob/master/data/zika/zika.fasta">Zika dataset</a>  
- <a href="https://github.com/4ment/phylostan/blob/master/examples/fluA/fluA.fa">fluA dataset</a>
- <a href="https://github.com/neherlab/treetime_examples/blob/master/data/h3n2_na/h3n2_na_20.fasta">h3n2_na dataset</a>
- <a href="">yeast dataset</a>


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

## Experiments

For each dataset, run 10 hill-climbing exploration with a random starting tree.

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

Summarizing: computing the parwise HOP similarity between the last tree obtained in each run.
For each dataset `DATASET` computing three files:  
- `results/<DATASET>/<DATASET>_best_trees.nwk`: last tree of each run in Newick format;  
- `results/<DATASET>/<DATASET>_best_trees.treevec`: last tree of each run in TreeVec format;  
- `results/<DATASET>/<DATASET>.dist`: pairwise HOP similarity between trees above.
```
experiments_hc > ./summary.sh zika 6388415
experiments_hc > ./summary.sh h3n2_na_20 6388418
experiments_hc > ./summary.sh yeast 6388420
experiments_hc > ./summary.sh fluA 6388416
run 2 terminated with an error
run 3 terminated with an error
```
Two runs of the `fluA` dataset did not finish in time, their final trees are not considered.
Overall we observe a wide range of similarity between the final trees of the runs, which is not a good result.