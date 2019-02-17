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
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def main():
    parser = ArgumentParser("AnalyzeSumMatrix",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument("--cell-type", required=True, help='Cell type, Ex: 1CDX1')
    parser.add_argument("--bin-size", required=True, help='Bin size, Ex: 1M, 500k')

    args = parser.parse_args()

    cell_type = args.cell_type
    bin_size = args.bin_size

    number_of_cells = 0

    # "1CDX1" (280 cells), "1CDX2" (303 cells), "1CDX3" (262 cells) and "1CDX4" (326 cells))
    if cell_type == "1CDX1":
        number_of_cells = 280
    elif cell_type == "1CDX2":
        number_of_cells = 303
    elif cell_type == "1CDX3":
        number_of_cells = 262
    elif cell_type == "1CDX4":
        number_of_cells = 326
    else:
        sys.exit("Cell type is not defined.")

    sum_max = 0
    if bin_size == "1M":
        if cell_type == "1CDX1":
            sum_max = 24083
        elif cell_type == "1CDX2":
            sum_max = 20065
        elif cell_type == "1CDX3":
            sum_max = 29153
        elif cell_type == "1CDX4":
            sum_max = 39006
    elif bin_size == "500k":
        if cell_type == "1CDX1":
            sum_max = 27539
        elif cell_type == "1CDX2":
            sum_max = 25419
        elif cell_type == "1CDX3":
            sum_max = 39575
        elif cell_type == "1CDX4":
            sum_max = 42029
    else:
        sys.exit("Bin size is not defined.")

    matrix_file = "analyze-{}/SumMatrix/sum_matrix_{}_{}.txt".format(bin_size, bin_size, cell_type)
    chrom_bin_range = "metadata/chrom_bins_range_{}.txt".format(bin_size)
    # output_file = "output/analyze/analyze_sum_matrix/{}/output_sum_matrix_{}_{}_v1.txt".format(bin_size, bin_size,
    #                                                                                                       cell_type)
    output_file = "output/analyze/analyze_sum_matrix/{}/threshold-0.1/output_sum_matrix_{}_{}_v1.txt".format(bin_size,
                                                                                                             bin_size,
                                                                                                             cell_type)

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
    # sys.exit(0)

    count1 = 0
    p_max = calculate_pmax(M, sum_max)
    # my_list = []
    data = np.triu(data, k=-1)
    out = open(output_file, "w")
    out.write("{} {} {}\n".format("bin1", "bin2", "count", "p_value"))
    # out.write("{} {} {} {}\n".format("bin1", "bin2", "count", "p_value"))

    for i in range(0, rows):
        for j in range(0, cols):
            if data[i][j] != 0:
                value = int(data[i][j])
                p_value = calculate_f(number_of_cells, value, p_max)
                # if p_value <= threshold(M):
                #     count1 += 1
                #     out.write("{} {} {} {}\n".format(i, j, value, p_value))
                # my_list.append(value)

                if value >= number_of_cells / 10:
                    count1 += 1
                    out.write("{} {} {}\n".format(i, j, value, p_value))

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


def calculate_pmax(M, sum_max):
    p_max = sum_max / M
    return p_max


def nCr(n, r):
    f = math.factorial
    return f(n) // f(r) // f(n - r)


def threshold(M):
    return 0.05 / M


def histogram(matrix):
    b = np.triu(matrix, k=-1)
    b = b.ravel()
    b = b[b >= 20]

    # unique, counts = np.unique(b, return_counts=True)

    plt.hist(b, bins=b.max().astype(int) - 1)
    plt.show()


if __name__ == '__main__':
    main()
