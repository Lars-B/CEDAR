"""
CEDAR-exploration: exploration of the tree space in a maximum likelihood framework
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
from Bio import SeqIO
from ete3 import Tree

from _raxml import raxml_loss
from TreeVec import TreeVec

def _read_fasta(data_path):
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

def _create_random_tree(leaves):
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

def _neighborhood_ML_scores(current_tree, fasta_path, DNA_model, tree_folder_path, iteration_counter):
    """
    Compute the likelihood of all trees in the neighbourhood of a tree using Raxml

    Input:
    - current_tree (TreeVec): current tree
    - fasta_path (str): pah to FASTA file
    - DNA_model (str): DNA evolution model to use in Raxml
    - tree_folder_path (str): pah to folder where Newick trees are written
    - iteration_counter (int): iteration counter in hill-climbing process
    Output:
    - ngb_scores (list((int,int),float)):
      list of the hop (i,j) and likelihood score of all trees in the neighbourhood
    """
    # Set of hops defining the neigbourhood of current tree
    ngb_hops = current_tree.hop_neighbourhood(export_list=True)
    # Output list
    ngb_scores = []
    # Loop on possible hops
    for (i,j) in ngb_hops:
        # New tree obtained by prforming hop (i,j)
        ngb_tree = current_tree.hop(i,j,inplace=False)
        # Liklihood score of new_tree
        ngb_tree_score = raxml_loss(
            fasta_path, ngb_tree, DNA_model,
            tree_folder_path, outfile=f"tmp_{iteration_counter}_{i}_{j}.tree"
        )
        ngb_scores.append(((i,j),ngb_tree_score))
    return ngb_scores

def best_ML_neighbour(current_tree, current_score, fasta_path, DNA_model, tree_folder_path, iteration_counter):
    """
    Chose the next tree as the best neighbouring tree in terms of likelihood score

    Input:
    - current_tree (TreeVec): current tree
    - current_score (float): score of current tree
    - fasta_path (str): pah to FASTA file
    - DNA_model (str): DNA evolution model to use in Raxml
    - tree_folder_path (str): pah to folder where Newick trees are written
    - iteration_counter (int): iteration counter in hill-climbing process
    Output:
    - best_ngb_tree (TreeVec): best neighbour tree
    - best_ngb_score (float): score of best neighbour tree
    """
    # Computing the likelihood of all trees in the neighbourhood
    ngb_scores_and_hops = _neighborhood_ML_scores(
        current_tree, fasta_path, DNA_model, tree_folder_path, iteration_counter
    )
    # Isolating liklihood scores in the same order than in ngb_scores_and_hops
    ngb_scores = [ngb_result[1] for ngb_result in ngb_scores_and_hops]
    # Index of the best likelihood score
    best_score_idx = np.argmax(ngb_scores)
    # HOP defining the best tree
    (i,j) = ngb_scores_and_hops[best_score_idx][0]
    # Creating the best tree by perforing the corresponding hop
    best_ngb_tree = current_tree.hop(i,j,inplace=False)
    best_ngb_score = ngb_scores[best_score_idx]
    return best_ngb_tree,best_ngb_score

def hill_climbing(
        fasta_path, DNA_model, tree_folder_path,
        first_tree, rng,
        tol, max_patience, max_nb_iterations,
        out_file_path
):
    """
    Using a Hill-Climbing heuristic to find a maximum likelihood tree from a set of
    FASTA sequences.

    Algorithm:
    iterate
    - compute the likelihood of all trees in the hop-neighbourhood
    - select the highest likelihood tree
    - if its likelihood is within tol of the best tree so far:
      - reorder leaves of the current tree and decrease the patience counter
    - otherwise
      - the best tree becomes the current tree
      - the patience counter is set to max_patience
    until the max number of iterations is reached or the patence counter is 0

    Input:
    - fasta_path (str): pah to FASTA file
    - DNA_model (str): DNA evolution model to use in Raxml
    - tree_folder_path (str): pah to folder where Newick trees are written
    - first_tree (TreeVec): starting tree
    - rng (Generator): numpy random number generator
    - tol (float): tolerance used to define a stopping criterion
    - max_patience (int): maximum number of iterations of rordering leaves
    - max_nb_iterations (int): maximum number of iterations
    - out_file_path (str): path of the file where all explored trees are recorded
      together with their likelihood score
    Output: Explored trees written in out_file_path
    """
    out_file = open(out_file_path, "w")
    current_tree = first_tree
    current_score = raxml_loss(
        fasta_path, current_tree, DNA_model,
        tree_folder_path, outfile="tmp_0_.tree"
    )
    out_file.write(f"START\t{current_score}\t{current_tree.treevec2newick()}\n")
    stop = False
    patience_counter = max_patience
    iteration_counter = 1
    while not stop:
        best_ngb_tree,best_ngb_score = best_ML_neighbour(
            current_tree, current_score,
            fasta_path, DNA_model, tree_folder_path,
            iteration_counter
        )
        if best_ngb_score - current_score > tol:
            current_tree = best_ngb_tree
            current_score = best_ngb_score
            patience_counter = max_patience
            prefix = "NEIGHBOUR"
        else:
            current_tree = current_tree.reorder_leaves(rng)
            patience_counter -= 1
            prefix = "REORDER"
        stop = (
            (patience_counter==0)
            or
            (
                (max_nb_iterations is not None)
                and
                (iteration_counter>max_nb_iterations)
            )
        )
        iteration_counter += 1
        out_file.write(f"{prefix}\t{current_score}\t{current_tree.treevec2newick()}\n")

def random_walk(
        fasta_path, DNA_model, tree_folder_path,
        first_tree, rng,
        max_nb_iterations,
        out_file_path
):
    """
    Performing a random walk in the tree space for a given number of iterations

    Input:
    - fasta_path (str): pah to FASTA file
    - DNA_model (str): DNA evolution model to use in Raxml
    - tree_folder_path (str): pah to folder where Newick trees are written
    - first_tree (TreeVec): starting tree
    - rng (Generator): numpy random number generator
    - max_nb_iterations (int): maximum number of iterations (steps of the walk)
    - out_file_path (str): path of the file where all explored trees are recorded
      together with their likelihood score
    Output: Explored trees written in out_file_path
    """
    out_file = open(out_file_path, "w")
    current_tree = first_tree
    current_score = None
    current_score = raxml_loss(
        fasta_path, current_tree, DNA_model,
        tree_folder_path, outfile="tmp_0_.tree"
    )
    out_file.write(f"START\t{current_score}\t{current_tree.treevec2newick()}\n")
    for i in range(max_nb_iterations):
        new_tree = current_tree.random_hop(rng, inplace=False)
        current_tree = new_tree
        current_score = raxml_loss(
            fasta_path, current_tree, DNA_model,
            tree_folder_path, outfile=f"tmp_0_{i}.tree"
        )
        out_file.write(f"RANDOM\t{current_score}\t{current_tree.treevec2newick()}\n")

def main(args):
    # Read a FASTA file to record sequences (taxa) names
    leaves = [record.id for record in _read_fasta(args.fasta_path)]

    # Instantiate a random numbers generators
    rng = np.random.default_rng(args.random_seed)

    # Create a random tree to start from
    first_tree = _create_random_tree(leaves)

    if args.cmd == "hc":
        # Hill-Climbing exploration
        hill_climbing(
            args.fasta_path, args.DNA_model, args.tree_folder_path,
            first_tree, rng,
            args.tol, args.max_patience, args.max_nb_iterations,
            args.out_file_path
        )
    elif args.cmd == "rw":
        # Random Walk exploration
        random_walk(
            args.fasta_path, args.DNA_model, args.tree_folder_path,
            first_tree, rng,
            args.max_nb_iterations,
            args.out_file_path
        )

def _parse_arguments():
    description = "CEDAR-exploration: exploring the tree space"

    argparser = argparse.ArgumentParser(prog="CEDAR-exploration", description=description)

    subparsers = argparser.add_subparsers(title="commands", help="command help")

    # Hill-Climbing
    hc = subparsers.add_parser("hc", help="Hill-Climbing")
    hc.set_defaults(cmd="hc")
    hc.add_argument("--fasta_path", type=str, help="Input FASTA file")
    hc.add_argument("--DNA_model", type=str, help="Input DNA model for Raxml")
    hc.add_argument("--tree_folder_path", type=str, help="Path to folder where trees are written for Raxml (not created)")
    hc.add_argument("--tol", type=float, default=0.001, help="Tolerance to dcide if a best tree has been found")
    hc.add_argument("--max_patience", type=int, default=3, help="Max number of leaves reordering steps if no best neighbour is found")
    hc.add_argument("--max_nb_iterations", type=int, default=100000, help="Max number of iterations")
    hc.add_argument("--random_seed", type=int, default=0, help="Random number generator seed")
    hc.add_argument("--out_file_path", type=str, help="Path to output file")

    # Random Walk
    rw = subparsers.add_parser("rw", help="Random Walk")
    rw.set_defaults(cmd="rw")
    rw.add_argument("--fasta_path", type=str, help="Input FASTA file")
    rw.add_argument("--DNA_model", type=str, help="Input DNA model for Raxml")
    rw.add_argument("--tree_folder_path", type=str, help="Path to folder where trees are written for Raxml (not created)")
    rw.add_argument("--max_nb_iterations", type=int, default=100000, help="Max number of iterations")
    rw.add_argument("--random_seed", type=int, default=0, help="Random number generator seed")
    rw.add_argument("--out_file_path", type=str, help="Path to output file")

    return argparser.parse_args()

if __name__ == "__main__":

    args = _parse_arguments()

    main(args)
