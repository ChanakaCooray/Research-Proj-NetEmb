import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    analyze_dir = "analyze-1M"
    # analyze_dir = "test"
    analyze_cat = "1CDX1"
    # size_chromY = 16

    # chromBin_file = "metadata/chrom_bins_1M.txt"
    chrom_bin_range = "metadata/chrom_bins_range_1M.txt"
    # chrom_bin_range = "metadata/chrom_bins_range_500k.txt"

    # search for last index
    # with open(chromBin_file) as f:
    #     for last in f: pass
    #     n = int(last.split()[1]) + size_chromY

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            splitLine = line.split()
            bin_range[splitLine[0]] = (int(splitLine[2]), int(splitLine[3]))

    sum_max = 0
    for filename in os.listdir(analyze_dir):
        if not filename.startswith(analyze_cat):
            continue

        data = {}
        count = 0
        with open(os.path.join(analyze_dir, filename)) as f:
            for line in f:
                splitLine = line.split()
                edge1 = int(splitLine[0])
                edge2 = int(splitLine[1])

                intra_chrm = False
                for key, value in bin_range.items():
                    if value[0] <= edge1 <= value[1]:
                        if value[0] <= edge2 <= value[1]:
                            intra_chrm = True
                            break

                if not intra_chrm:
                    count += 1
                    if edge1 <= edge2:
                        data[(edge1, edge2)] = 1
                    else:
                        data[(edge2, edge1)] = 1

        if count > sum_max:
            sum_max = count

    print(sum_max)
        # print(count)
        #
        # matrix = np.zeros((n, n), dtype=np.int)
        #
        # for key, val in data.items():
        #     matrix[key[0]][key[1]] = val
        #     matrix[key[1]][key[0]] = val
        #
        # print(matrix.sum(axis=0).sum(axis=0))
        #
        # ax = sns.heatmap(matrix, cmap="Reds")
        #
        # plt.show()

        # np.fill_diagonal(matrix, 0)


if __name__ == '__main__':
    main()
