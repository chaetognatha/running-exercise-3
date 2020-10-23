import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

inputs = ["mtDNAtable1","Ychrtable1","mtDNAtable2","Ychrtable2"]


def matrix_parser(fh):
    my_counter = 0
    m = []
    with open(fh) as f:
        for line in f:
            if my_counter == 0:
                names = list(line.strip("\n").split())
                my_counter += 1
            else:
                row = list(line.strip("\n").split())
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

    labelList = names

    plt.figure(figsize=(10, 7))
    dendrogram(linked,
               orientation='top',
               labels=labelList,
               distance_sort='descending',
               show_leaf_counts=True)
    plt.show()


# make_tree("julia_out.txt")
def main():
    for f in inputs:
        make_tree(f)
main()