# CEDAR: Gelman-Rubin MCMC convergence

This folder contains the results of experiments on using the Gelman-Rubin statistics described in 
<a href="https://doi.org/10.1109/TCBB.2024.3457875">An Automated Convergence Diagnostic for Phylogenetic MCMC Analyses</a> 
and implemented in CEDAR with the HOP distance instead of the RNNI distance.

## Data and methods
The analyzed data was composed of two <a href="https://www.beast2.org/">BEAST</a> runs, files `ds1-r1.newick.gz` and `ds1-r2.newick.gz`.

We implemented the method described in <a href="https://doi.org/10.1109/TCBB.2024.3457875">An Automated Convergence Diagnostic for Phylogenetic MCMC Analyses</a>
with a single difference: for the computation of the distance between a pair of trees, instead of the ranked NNI (RNNI) distance, we take the minimum, over a prescribed number of random leaves orders, HOP distance;
this approach aims to leverage the fast computation of the HOP distance to obtain an estimate of the SPR distance.
The random leaves orders are detemined at the beginning of the process and are the same for all distance computations.

## Experiments
For the first experiment, we take the last 1000 trees of each BEAST run and consider 5 random leaves orders for the distance computation:
```
gunzip ds1-r1.newick.gz
gunzip ds1-r2.newick.gz
python ../src/CEDAR.py GR \
  --Newick_file_1 ds1-r1.newick --Newick_file_2 ds1-r2.newick \
  --nb_trees 1000 --nb_orders 5 --seed 1 \
  --output_gr_file GR_1000.tsv --output_orders_file orders_1000.tsv
gzip ds1-r1.newick
gzip ds1-r2.newick
```
The output files are `GR_1000.tsv` that records the Gelman-Rubin statistic value for each of the 1000 pair of trees from the two BEAST runs
and `orders_1000.tsv` that recors the 5 random leaves orders (to allow reproducibility of the experiment).

The second experiment is similar to the first one but considers the last 100 trees of each BEAST run and 3 random leaves orders.
The output files are `GR_100.tsv` and `orders_100.tsv`
```
gunzip ds1-r1.newick.gz
gunzip ds1-r2.newick.gz
python ../src/CEDAR.py GR \
  --Newick_file_1 ds1-r1.newick --Newick_file_2 ds1-r2.newick \
  --nb_trees 100 --nb_orders 3 --seed 1 \
  --output_gr_file GR_100.tsv --output_orders_file orders_100.tsv
gzip ds1-r1.newick
gzip ds1-r2.newick
```
