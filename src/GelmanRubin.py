"""
CEDAR: manipulating phylogenetic rooted trees representations as vectors
Convergence of an MCMC
"""

__author__ = "Cedric Chauve"
__credits__ = ["Cedric Chauve", "Louxin Zhang"]
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"

import argparse
import numpy as np

from utils import (
    __read_file
)
from LeavesOrder import (
    order2str
)
from TreeVec import (
    TreeVec,
    get_nb_taxa
)

def __hop_distance_within(in_TreeVec_trees_1):
    """
    Given a list of trees, compute the hop distance between all pairs of trees
    - in_TreeVec_trees_1: list(TreeVec) list of TreeVec objects
    Output:
    - np.array where entry [i][j] is the distance between trees i from list 1 and j from list 1 (0-base index)
    """
    nb_trees = len(in_TreeVec_trees_1)
    nb_taxa = get_nb_taxa(in_TreeVec_trees_1[0])
    distances = np.zeros([nb_trees, nb_trees])
    for i in range(0,nb_trees-1):
        tree1 = in_TreeVec_trees_1[i]
        range_j = [j for j in range(i+1,nb_trees)]
        for j in range_j:
            tree2 = in_TreeVec_trees_1[j]
            similarity = tree1.hop_similarity(tree2, compute_seq=False)
            distances[i][j] = nb_taxa - similarity
            distances[j][i] = distances[i][j]
    return distances

def __hop_distance_between(in_TreeVec_trees_1, in_TreeVec_trees_2):
    """
    Given two list of TreeVec trees, compute the hop distance between trees of list 1 an trees of list 2
    - in_TreeVec_trees_1: list(TreeVec) list of TreeVec objects
    - in_TreeVec_trees_2: list(TreeVec) list of TreeVec objects
    Output:
    - np.array where entry [i][j] is the distance between trees i from list 1 and j from list 2 (0-base index)
    """
    nb_trees = len(in_TreeVec_trees_1)
    nb_taxa = get_nb_taxa(in_TreeVec_trees_1[0])    
    distances = np.zeros([nb_trees, nb_trees])
    for i in range(0,nb_trees-1):
        tree1 = in_TreeVec_trees_1[i]
        for j in range(0,nb_trees):
            tree2 = in_TreeVec_trees_2[j]
            similarity = tree1.hop_similarity(tree2, compute_seq=False)
            distances[i][j] = nb_taxa - similarity
    return distances

def _GelmanRubin(in_TreeVec_trees_1, in_TreeVec_trees_2):
    """
    Given two lists of trees (chains), compute the Gelman Rubin diagnostic value for each tree of each chain
    Assumptions:
    - both chains have the same number of trees
    - all trees are on the same taxon set
    - all TreeVec strings have been obtained with the same leaves order 

    """
    nb_trees = len(in_TreeVec_trees_1)

    squared_distances_within = {}
    squared_distances_between = {}
    # squared distances within chain 1
    squared_distances_within[1] = np.square(__hop_distance_within(in_TreeVec_trees_1))
    # squared distances within chain 2
    squared_distances_within[2] = np.square(__hop_distance_within(in_TreeVec_trees_2))
    # squared distances between chain 1 (index 1) and chain 2 (index 2)
    squared_distances_between[1] = np.square(__hop_distance_between(in_TreeVec_trees_1, in_TreeVec_trees_2))
    # squared distances between chain 2 (index 1) and chain 1 (index 2)
    squared_distances_between[2] = np.transpose(squared_distances_between[1])
    # Remark: waste of space

    # GelmanRubin value for tree of rank i in each chain
    GR =  {c: np.zeros(nb_trees) for c in [1,2]}
    for c1 in [1,2]:
        # Processing trees of chain c1
        for i in range(1,nb_trees):
            # Processing tree i
            for s in range(i+1):
                # Remark: overkill as sums could be updated dynamically
                sum_1 = np.sum(squared_distances_between[c1][s][:i+1])
                sum_2 = np.sum(squared_distances_within[c1][s][:i+1])
                PSRF = np.sqrt(sum_1/sum_2)
                GR[c1][i] += PSRF
            GR[c1][i] /= np.float64(i+1)
    return GR

def GelmanRubin(in_Newick_trees_file_1, in_Newick_trees_file_2, nb_trees, out_file):

    in_Newick_trees_1 = __read_file(in_Newick_trees_file_1)
    in_Newick_trees_2 = __read_file(in_Newick_trees_file_2)
    assert len(in_Newick_trees_1)>=nb_trees, "Chain 1 too short"
    assert len(in_Newick_trees_2)>=nb_trees, "Chain 2 too short"
    
    _tree = TreeVec(newick_str=f"{in_Newick_trees_1[0].rstrip()};")
    leaf2idx,idx2leaf = _tree.extract_leaves_order()    
    leaves_order_str = order2str(idx2leaf)
    # TO DO: random leaves order

    in_trees_1 = [
        TreeVec(newick_str=f"{in_Newick_tree.rstrip()};",leaf2idx=leaf2idx)
        for in_Newick_tree in in_Newick_trees_1[-nb_trees:]
    ]
    in_trees_2 = [
        TreeVec(newick_str=f"{in_Newick_tree.rstrip()};",leaf2idx=leaf2idx)
        for in_Newick_tree in in_Newick_trees_2[-nb_trees:]
    ]

    GR = _GelmanRubin(in_trees_1, in_trees_2)

    with open(out_file, "w") as _out_file:
        _out_file.write(f"#{leaves_order_str}\n")
        for i in range(nb_trees):
            _out_file.write(f"{i},{GR[1][i]},{GR[2][i]}\n")

if __name__ == "__main__":
    description = "CEDAR: manipulating phylogenetic rooted trees representations as vectors; Gelman Rubin diagnostic test for MCMC convergence"

    argparser = argparse.ArgumentParser(prog="CEDAR-GR", description=description)
    subparsers = argparser.add_subparsers(title="commands", help="command help")

    # Creating Gelman Rubin diagnostic values for two chains
    GR = subparsers.add_parser("GR", help="Gelman Rubin diagnosic values")
    GR.set_defaults(cmd="GR")
    GR.add_argument("--Newick_file_1", type=str, help="Input Newick file 1")
    GR.add_argument("--Newick_file_2", type=str, help="Input Newick file 2")
    GR.add_argument("--nb_trees", type=int, help="Number of trees (trail of files to consider")
    GR.add_argument("--output_file", type=str, help="Output CEDAR file")
    
    args = argparser.parse_args()

    GelmanRubin(
        args.Newick_file_1,
        args.Newick_file_2,
        args.nb_trees,
        args.output_file
    )
    
