import numpy as np
import os
from scipy import sparse
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm


def main():

    matrix_file = "analyze-1M/sum_matrix_1M_1CDX1.txt"
    # matrix_file = "analyze-1M/sum_matrix_1M_1CDX2.txt"
    # matrix_file = "analyze-1M/sum_matrix_1M_1CDX3.txt"
    # matrix_file = "analyze-1M/sum_matrix_1M_1CDX4.txt"
    # matrix_file = "analyze-500k/sum_matrix_500k_1CDX1.txt"
    # matrix_file = "analyze-500k/sum_matrix_500k_1CDX2.txt"
    # matrix_file = "analyze-500k/sum_matrix_500k_1CDX3.txt"
    # matrix_file = "analyze-500k/sum_matrix_500k_1CDX4.txt"

    data = np.genfromtxt(matrix_file, delimiter=" ")

    ax = sns.heatmap(data, cmap="Reds", robust=True)

    plt.show()

    # print(data)

if __name__ == '__main__':
    main()
