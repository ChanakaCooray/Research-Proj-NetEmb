import numpy as np
import os
from scipy import sparse
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm

# matplotlib.use('TkAgg')


def main():
    # chromosome bins metada
    chromBin_file = "metadata/chrom_bins.txt"
    # chromosome sizes file

    # data file
    # data_file = "analyze/1CDES_p3.H6"
    # data_file = "analyze-1M/1CDU.44"
    # data_file = "analyze-1M/1CDX1.71"
    # data_file = "analyze-1M/1CDX2.372"
    data_file = "analyze-1M/1CDX4.206"

    # output file
    # output_file = "output/matrix/1CDU.1.out"

    # os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # search for last index
    with open(chromBin_file) as f:
        for last in f: pass
        n = int(last.split()[1]) + 15

    # search for last index
    # with open(chromBin_file) as f:
    #     for last in f: pass
    # print(n)

    data = {}
    with open(data_file) as f:
        for line in f:
            splitLine = line.split()
            edge1 = int(splitLine[0])
            edge2 = int(splitLine[1])
            # if edge1 != edge2:
            data[(edge1, edge2)] = int(splitLine[2])

    # df = pd.read_csv(data_file, sep=" ", header=None, names=["bin1", "bin2", "val"])

    matrix = np.zeros((n, n), dtype=np.int)

    for key, val in data.items():
        matrix[key[0]][key[1]] = val
        matrix[key[1]][key[0]] = val


    # print(matrix)

    # print(matrix.sum(axis=1))

    #
    # np.fill_diagonal(matrix, 0)
    my_sum = matrix.sum(axis=1)

    np.savetxt('test.txt', my_sum)

    # fig, axes = plt.subplots(nrows=2, ncols=2)
    # fig.tight_layout()
    #
    # ax1 = axes[0][0]
    #
    # ax1.hist(sum, bins=sum.max() - 1)
    #
    # plt.show()

    # print(matrix[2554][2549])

    # plt.style.use('seaborn-whitegrid')
    # df = df.drop(df[df.bin1 == df.bin2].index)
    # plt.scatter(df["bin1"], df["bin2"], s=df["val"],c=df["val"], cmap="viridis")
    # plt.colorbar()
    # plt.show()

    # np.fill_diagonal(matrix, 0)

    # ax = sns.heatmap(matrix, robust=True)
    # plt.figure(figsize=(100, 100))
    # plt.imshow(matrix, aspect='auto', cmap=plt.get_cmap("Reds"))
    # plt.figure(figsize=(100, 100))
    # fig.set_size_inches(20, 15)
    # plt.savefig("output/analyze/heatmap.png")

    # sparse_mat = sparse.csr_matrix(matrix)
    #
    # print(sparse_mat)

    # np.savetxt(output_file, matrix, fmt="%d")

    # out = open(output_file, "w")
    # for row in matrix:
    #     np.savetxt(output_file, row)
    # out.write(row)

    # out.close()

    # plt.show()


if __name__ == '__main__':
    main()
