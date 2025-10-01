"""
CEDAR: manipulating phylogenetic rooted trees representations as vectors
Hill Climbing heuristic to find a tree from squence data
"""

__author__ = "Cedric Chauve"
__credits__ = ["Cedric Chauve", "Louxin Zhang"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"

import argparse
import numpy as np
from Bio import SeqIO
from ete3 import Tree

from _raxml import raxml_loss
from TreeVec import TreeVec

def read_fasta(data_path):
    """
    Read a fasta file

    Input:
    - data_path (str): Path to fasta file
    Output:
    - Generator
        An iterator for all sequences
        Each element has "id" (taxon name) and "seq" (sequence) attributes
    """
    return SeqIO.parse(data_path, "fasta")

def create_random_tree(leaves):
    """
    Create a random tree as a TreeVec object

    Input:
    - leaves (list(str)): List of leaves labels
    Output:
    - TreeVec: TeeVec object
    """    
    tree_aux = Tree()
    nb_leaves = len(leaves)
    tree_aux.populate(nb_leaves, names_library=leaves)
    tree = TreeVec(tree=tree_aux)
    return tree

def optimization_step(current_tree, fasta_path, current_best_score, DNA_model, tree_folder_path):
    """
    Chose the best neighbour

    Input:
    - tree (TreeVec): current tree
    - rng: random numbers generator
    Output:
    - TreeVec: new tree
    """
    tree_ngbs = current_tree.hop_neighbourhood()
    found_better_tree = False
    best_ngb,best_score = None,current_best_score
    for tree in tree_ngbs:
        current_score = raxml_loss(
            fasta_path,
            tree,
            DNA_model,
            tree_folder_path,
            outfile="tmp.tree"
        )
        if current_score > current_best_score:
            found_better_tree = True
            best_ngb = tree
            best_score = current_score

    return best_score,best_ngb
    
def main(
        fasta_path,
        random_seed,
        DNA_model,
        tree_folder_path,
        patience_max
):

    # Read a FASTA file
    sequences = read_fasta(fasta_path)
    leaves = [record.id for record in sequences]
    nb_leaves = len(leaves)

    # Create a random tree
    current_tree = create_random_tree(leaves)
    current_score = raxml_loss(
        fasta_path,
        current_tree,
        DNA_model,
        tree_folder_path,
        outfile="tmp.tree"
    )
    
    # Instantiate a random number generator
    rng = np.random.default_rng(random_seed)
    
    # Loop of optimization
    stop = False
    patience = patience_max
    print(current_tree.treevec2str())
    print(current_score)
    while not stop:
        print("--> new iteration")
        best_score,best_ngb = optimization_step(
            current_tree, fasta_path, current_score, DNA_model, tree_folder_path
        )
        if best_ngb is None:
            print("No better neighbour found")
            _current_tree = current_tree.treevec2tree()
            rng.shuffle(leaves)
            _leaf2idx = {
                leaves[i]: i+1 for i in range(nb_leaves)
            }
            current_tree = TreeVec(tree=_current_tree, leaf2idx=_leaf2idx)
            patience -= 1
        else:
            print("Better neighbour found")
            current_tree = best_ngb
            current_score = best_score
            patience = patience_max
        print(current_tree.treevec2str())
        print(current_score)
        stop = (patience == 0)

main(
    fasta_path="test.fasta",
    random_seed=0,
    DNA_model="JC",
    tree_folder_path="tree_folder",
    patience_max=5
)
