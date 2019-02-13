import os
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def main():
    parser = ArgumentParser("SumMatrix",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument("--analyze-dir", required=True)
    parser.add_argument("--cell-type", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--bin-size", required=True)
    parser.add_argument("--size-chromY", required=True)
    parser.add_argument("--output-dir", required=True)

    args = parser.parse_args()
    analyze_dir = args.analyze_dir
    analyze_cat = args.cell_type
    metadata = args.metadata
    bin_size = args.bin_size
    size_chromY = int(args.size_chromY)
    output_dir = args.output_dir

    output_file = os.path.join(output_dir, "sum_matrix_{}_{}.txt".format(bin_size, analyze_cat))

    # size_chromY = 16

    chromBin_file = "{}/chrom_bins_{}.txt".format(metadata, bin_size)

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
