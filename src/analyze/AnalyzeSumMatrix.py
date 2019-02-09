import sys
import numpy as np
import os
from scipy import sparse
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm
import math


def main():
    # matrix_file = "analyze-1M/sum_matrix_1M_1CDX1.txt"
    # matrix_file = "analyze-1M/sum_matrix_1M_1CDX2.txt"
    # matrix_file = "analyze-1M/sum_matrix_1M_1CDX3.txt"
    # matrix_file = "analyze-1M/sum_matrix_1M_1CDX4.txt"
    matrix_file = "analyze-500k/sum_matrix_500k_1CDX1.txt"
    # matrix_file = "analyze-500k/sum_matrix_500k_1CDX2.txt"
    # matrix_file = "analyze-500k/sum_matrix_500k_1CDX3.txt"
    # matrix_file = "analyze-500k/sum_matrix_500k_1CDX4.txt"

    # chrom_bin_range = "metadata/chrom_bins_range_1M.txt"
    chrom_bin_range = "metadata/chrom_bins_range_500k.txt"

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            splitLine = line.split()
            bin_range[splitLine[0]] = (int(splitLine[2]), int(splitLine[3]))

    # print(bin_range)

    data = np.genfromtxt(matrix_file, delimiter=" ")

    rows = data.shape[0]
    cols = data.shape[1]

    # print(rows)
    # print(cols)
    count = 0
    for i in range(0, rows):
        for j in range(0, cols):
            for key, value in bin_range.items():
                if value[0] <= i <= value[1]:
                    if value[0] <= j <= value[1]:
                        data[i][j] = 0
                        count += 1

    M = (rows * cols - count)

    # unique, counts = np.unique(data, return_counts=True)
    # print(unique)
    # print(unique.shape)

    count1 = 0
    p_max = calculate_pmax(M)
    my_list = []
    data = np.triu(data, k=-1)
    output_file = "output/analyze/output_sum_matrix_500k.txt"
    out = open(output_file, "w")
    out.write("{} {} {}\n".format("bin1", "bin2", "count"))
    # print(p_max)
    for i in range(0, rows):
        for j in range(0, cols):
            if data[i][j] != 0:
                value = int(data[i][j])

                if calculate_f(280, value, p_max) <= threshold(M):
                    count1 += 1
                    out.write("{} {} {}\n".format(i, j, value))
                    my_list.append(value)
    out.close()
    my_arr = np.array(my_list)
    # unique1, counts1 = np.unique(my_arr, return_counts=True)
    # print(unique1)
    # print(unique1.shape)
    # print(myset)
    # print(count1)

    # np.savetxt("analyze-1M/sum_matrix_1M_1CDX1_withoutDiag.txt", data, fmt='%d', delimiter=' ', newline='\n')

    # histogram(data)
    # print(data.shape)
    # sys.exit(0)

    # print(count)

    # sum_by_row = np.triu(data, k=-1).sum(axis=0).astype(int)
    # total_sum = sum_by_row.sum(axis=0).astype(int)
    #
    # print(total_sum)

    # A = -np.sort(-sum_by_row)
    # B = np.argsort(-sum_by_row)

    # print(data[data != 0].shape)
    # print(B.shape)

    # print(total_sum)
    # print(total_sum / data[data != 0].ravel().shape[0])

    # ax = sns.heatmap(data, cmap="Reds")
    #
    # plt.show()


def calculate_f(n, t, p_max):
    return nCr(n, t) * (p_max ** t) * ((1 - p_max) ** (n - t))


def calculate_pmax(M):
    sum_max = get_sum_max()
    # M = get_M()

    p_max = sum_max / M
    return p_max


def nCr(n, r):
    f = math.factorial
    return f(n) // f(r) // f(n - r)


def threshold(M):
    # M = get_M()
    return 0.05 / M

def get_sum_max():
    # return 24083    #for 1M cdx1
    return 27539

# def get_M():
#     # return 6730306      #for 1M cdx1
#     return 0

def histogram(matrix):
    b = np.triu(matrix, k=-1)
    b = b.ravel()
    b = b[b >= 20]

    # unique, counts = np.unique(b, return_counts=True)

    # print(b.max().astype(int))
    # print(unique)

    # print(b.shape)

    plt.hist(b, bins=b.max().astype(int) - 1)
    plt.show()


if __name__ == '__main__':
    # calculateProb()
    main()
    # print(threshold())
    # print(calculate_pmax())