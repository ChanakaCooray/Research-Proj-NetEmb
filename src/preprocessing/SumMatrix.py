import os
import numpy as np
import matplotlib.pyplot as plt


def main():
    analyze_dir = "analyze-1M"
    # analyze_dir = "test"
    analyze_cat = "1CDX1"
    size_chromY = 16

    chromBin_file = "metadata/chrom_bins.txt"

    # search for last index
    with open(chromBin_file) as f:
        for last in f: pass
        n = int(last.split()[1]) + size_chromY

    sum_matrix = np.zeros((n, n), dtype=np.int)
    for filename in os.listdir(analyze_dir):

        if not filename.startswith(analyze_cat):
            continue

        data = {}
        with open(os.path.join(analyze_dir, filename)) as f:
            for line in f:
                splitLine = line.split()
                edge1 = int(splitLine[0])
                edge2 = int(splitLine[1])
                if(edge1<=edge2):
                    data[(edge1, edge2)] = 1
                else:
                    data[(edge2, edge1)] = 1

        matrix = np.zeros((n, n), dtype=np.int)

        for key, val in data.items():
            matrix[key[0]][key[1]] = val
            matrix[key[1]][key[0]] = val

        np.fill_diagonal(matrix, 0)

        sum_matrix = sum_matrix + matrix

    np.savetxt("test.txt", sum_matrix, fmt='%d', delimiter=' ', newline='\n')

if __name__ == '__main__':
    main()
