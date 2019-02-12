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
    matrix_file = "analyze-1M/sum_matrix_1M_1CDX1.txt"
    # matrix_file = "analyze-1M/sum_matrix_1M_1CDX2.txt"
    # matrix_file = "analyze-1M/sum_matrix_1M_1CDX3.txt"
    # matrix_file = "analyze-1M/sum_matrix_1M_1CDX4.txt"
    # matrix_file = "analyze-500k/sum_matrix_500k_1CDX1.txt"
    # matrix_file = "analyze-500k/sum_matrix_500k_1CDX2.txt"
    # matrix_file = "analyze-500k/sum_matrix_500k_1CDX3.txt"
    # matrix_file = "analyze-500k/sum_matrix_500k_1CDX4.txt"

    chrom_bin_range = "metadata/chrom_bins_range_1M.txt"
    # chrom_bin_range = "metadata/chrom_bins_range_500k.txt"

    number_of_cells = 280

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            splitLine = line.split()
            bin_range[splitLine[0]] = (int(splitLine[2]), int(splitLine[3]))

    data = np.genfromtxt(matrix_file, delimiter=" ")

    rows = data.shape[0]
    cols = data.shape[1]

    zero_bins = 0
    sum_by_row = data.sum(axis=0).astype(int)
    # print(sum_by_row.shape)
    zero_bin_list = {}
    for i, val in enumerate(sum_by_row):
        if val == 0:
            zero_bins += 1
            for key, value in bin_range.items():
                if value[0] <= i <= value[1]:
                    zero_bin_list[i] = key

    # zero_bin_list_out = "output/analyze/zero_bin_list.txt"
    # out_zero_bin = open(zero_bin_list_out, "w")
    # out_zero_bin.write("{} {}\n".format("bin", "chrm"))
    # for key, val in zero_bin_list.items():
    #     out_zero_bin.write("{} {}\n".format(key, val))
    # out_zero_bin.close()

    df_zero_bin = pd.DataFrame(list(zero_bin_list.items()), columns=['bin', 'chrm'])
    df_zero_bin_count = df_zero_bin['chrm'].value_counts().reset_index()
    df_zero_bin_count.columns = ['chrm', 'count']
    df_zero_bin_count.set_index('chrm', inplace=True)

    # M = ( (N-X)^2 - (n1-x1)^2 - (n2-x2)^2 - ... - (n20-x20)^2 - (n21-x21)^2 )/2
    sum_reduction = 0  # stores summation of (n2-x2)^2 + ... + (n20-x20)^2 + (n21-x21)^2
    count_X = 0  # stores X
    count_N = 0  # stores N
    for key, value in bin_range.items():
        n = value[1] - value[0] + 1  # store the number of bins for the chromosome
        zero_bins_chrm = df_zero_bin_count.loc[key, 'count']  # stores zero bins for the chrm
        sum_reduction = sum_reduction + (n - zero_bins_chrm) ** 2  # cal (n2-x2)^2 and add that to sum
        count_X = count_X + zero_bins_chrm
        count_N = count_N + n

    count = 0
    for i in range(0, rows):
        for j in range(0, cols):
            for key, value in bin_range.items():
                if value[0] <= i <= value[1]:
                    if value[0] <= j <= value[1]:
                        data[i][j] = 0
                        count += 1

    # unique, counts = np.unique(data, return_counts=True)

    M = ((count_N - count_X) ** 2 - sum_reduction) / 2

    # print(M)

    count1 = 0
    p_max = calculate_pmax(M)
    # my_list = []
    data = np.triu(data, k=-1)
    output_file = "output/analyze/output_sum_matrix_1M.txt"
    out = open(output_file, "w")
    out.write("{} {} {}\n".format("bin1", "bin2", "count"))

    for i in range(0, rows):
        for j in range(0, cols):
            if data[i][j] != 0:
                value = int(data[i][j])

                if calculate_f(number_of_cells, value, p_max) <= threshold(M):
                    count1 += 1
                    out.write("{} {} {}\n".format(i, j, value))
                    # my_list.append(value)
    out.close()

    # my_arr = np.array(my_list)
    # unique1, counts1 = np.unique(my_arr, return_counts=True)
    # print(unique1)
    # np.savetxt("analyze-1M/sum_matrix_1M_1CDX1_withoutDiag.txt", data, fmt='%d', delimiter=' ', newline='\n')
    # histogram(data)
    # sum_by_row = np.triu(data, k=-1).sum(axis=0).astype(int)
    # total_sum = sum_by_row.sum(axis=0).astype(int)
    # A = -np.sort(-sum_by_row)
    # B = np.argsort(-sum_by_row)
    # print(total_sum / data[data != 0].ravel().shape[0])
    # ax = sns.heatmap(data, cmap="Reds")
    # plt.show()


def calculate_f(n, t, p_max):
    return nCr(n, t) * (p_max ** t) * ((1 - p_max) ** (n - t))


def calculate_pmax(M):
    sum_max = get_sum_max()

    p_max = sum_max / M
    return p_max


def nCr(n, r):
    f = math.factorial
    return f(n) // f(r) // f(n - r)


def threshold(M):
    return 0.05 / M


def get_sum_max():
    return 24083  # for 1M cdx1
    # return 27539 #for 1M cdx1


def histogram(matrix):
    b = np.triu(matrix, k=-1)
    b = b.ravel()
    b = b[b >= 20]

    # unique, counts = np.unique(b, return_counts=True)

    plt.hist(b, bins=b.max().astype(int) - 1)
    plt.show()


if __name__ == '__main__':
    main()
