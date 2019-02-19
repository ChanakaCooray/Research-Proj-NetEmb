import os
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys


def main():
    parser = ArgumentParser("SumMatrix",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument("--analyze-dir", required=True)
    parser.add_argument("--cell-type", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--bin-size", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--shift", default='0')

    args = parser.parse_args()
    analyze_dir = args.analyze_dir
    analyze_cat = args.cell_type
    metadata = args.metadata
    bin_size = args.bin_size
    shift = args.shift
    output_dir = args.output_dir

    # output_file = os.path.join(output_dir, "sum_matrix_{}_{}.txt".format(bin_size, analyze_cat))

    # create output directory if not exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if shift == '0':
        output_file = os.path.join(output_dir, "sum_matrix_{}_{}.txt".format(bin_size, analyze_cat))
        chrom_bin_range = "{}/chrom_bins_range_{}.txt".format(metadata, bin_size)
    else:
        output_file = os.path.join(output_dir, "sum_matrix_{}_{}_shift_{}.txt".format(bin_size, analyze_cat, shift))
        chrom_bin_range = "{}/chrom_bins_range_{}_shift_{}.txt".format(metadata, bin_size, shift)

    # search for last index
    with open(chrom_bin_range) as f:
        for last in f: pass
        n = int(last.split()[3]) + 1

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
                if (edge1 <= edge2):
                    data[(edge1, edge2)] = 1
                else:
                    data[(edge2, edge1)] = 1

        matrix = np.zeros((n, n), dtype=np.int)

        for key, val in data.items():
            matrix[key[0]][key[1]] = val
            matrix[key[1]][key[0]] = val

        # np.fill_diagonal(matrix, 0)

        sum_matrix = sum_matrix + matrix

    np.savetxt(output_file, sum_matrix, fmt='%d', delimiter=' ', newline='\n')


if __name__ == '__main__':
    main()
