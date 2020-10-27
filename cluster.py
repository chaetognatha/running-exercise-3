#!/usr/bin/python3
'''
    Title: cluster.py
    Date: 2020-10-23
    Author: Mattis Knulst
    Description:
        This program will parse a file with a similarity matrix in squareform, then generate a dendrogram from it
        that can be interpreted as a phylogenetic tree.
    List of functions:

    List of non-standard modules:
        Pyplot from matplotlib
        Numpy
        Dendrogram, linkage and pdist from scipy
    Procedure:
        The 4 input files are located in a list that is first parsed to create numpy arrays from the files that can
        then be tested for linkage and used to generate a tree. Thus 4 trees are created.

    Usage:
        ./cluster.py
'''
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

inputs = ["mt_out_identity","mt_out_align","y_out_identity","y_out_align"]


def matrix_parser(fh):
    my_counter = 0
    m = []
    with open(fh) as f:
        for line in f:
            if my_counter == 0:
                names = list(line.strip("\n").split("\t"))
                my_counter += 1
            else:
                row = list(line.strip("\n").split("\t"))
                m.append(row[1:])
    matrix = np.array(m, dtype=float)
    return names, matrix


def make_tree(fh):
    names, matrix = matrix_parser(fh)
    distance_matrix = pdist(matrix)
    print(names)
    print(distance_matrix)

    """
    Hierarchical clustering
    """
    linked = linkage(distance_matrix, 'single')

    labelList = names[1:]

    plt.figure(figsize=(25, 15))
    dendrogram(linked,
               orientation='left',
               labels=labelList,
               distance_sort='descending',
               show_leaf_counts=True)
    #plt.show()
    thing = fh + ".png"
    plt.savefig(thing)


# make_tree("julia_out.txt")
def main():
    counter = 0
    for f in inputs:
        make_tree(f)
main()
