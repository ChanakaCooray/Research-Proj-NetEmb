import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
import os
from sklearn.preprocessing import normalize
import numpy as np


# matplotlib.use('TkAgg')


def main():
    # analyze_dir = "analyze"

    # for filename in os.listdir(analyze_dir):
    # g = nx.Graph()
    # with open(os.path.join(analyze_dir, filename)) as f:
    #     for line in f:
    #         splitLine = line.split()
    #         edge1 = int(splitLine[0])
    #         edge2 = int(splitLine[1])
    #         weight = int(splitLine[2])
    #         # if edge1 != edge2:
    #         #     g.add_edge(edge1, edge2, weight=weight)
    #         g.add_edge(edge1, edge2, weight=weight)

    # g.add_edge(1,2,weight=3)
    # g.add_edge(2,3,weight=3)
    # g.add_edge(4,4,weight=3)
    #
    #
    # # print("Cell: {}".format(filename))
    # print("# of edges: {}".format(g.number_of_edges()))
    # print("# of nodes: {}".format(g.number_of_nodes()))
    # print("Connected: {}".format(nx.is_connected(g)))
    # print("Connected Components: {}".format(nx.number_connected_components(g)))
    # print("\n")

    # for i in range(0, 5):
    #     print(i)

    # matrix = np.zeros((3, 3), dtype=np.int)

    # for key, val in data.items():
    #     matrix[key[0]][key[1]] = val
    #     matrix[key[1]][key[0]] = val

    # matrix = np.matrix([[8, 4, 0], [4, 2, 0], [0, 0, 9]])

    # matrix = np.arange(0,27,3).reshape(3,3).astype(np.float64)

    # array([[  0.,   3.,   6.],
    #   [  9.,  12.,  15.],
    #   [ 18.,  21.,  24.]])

    # print(matrix)

    # normed_matrix = normalize(matrix, axis=0, norm='l1')

    # print(normed_matrix)

    # a = [[3,4,5],[6,7,8],[9,10,11],[12,13,14]]
    a = np.array([[3,4,5],[6,7,8],[9,10,11],[12,13,14]])
    # print(a[a > 50])

    print(a[:2,:])

    # print(np.linalg.norm(matrix, ord=1))


if __name__ == '__main__':
    main()
