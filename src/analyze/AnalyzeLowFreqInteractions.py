import os
from concurrent.futures import ThreadPoolExecutor
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys
import pandas as pd
import numpy as np
import math
from scipy.stats import binom
import shutil
import tempfile
import matplotlib.pyplot as plt


# determine the bin for the given coordinate and the sliding window
def find_bin(coord, start_index, bin_size, shift):
    bin_num = coord // bin_size
    remainder = coord % bin_size

    bin_num += start_index

    if shift == 0:
        return bin_num

    if remainder >= shift:
        bin_num += 1

    return bin_num


# generate the starting bin indices meta file
def generate_bins(chrom_sizes_file, metadata, bin_size, shift):
    if shift == '0':
        output = os.path.join(metadata, 'chrom_bins_{}.txt'.format(bin_size))
    else:
        output = os.path.join(metadata, 'chrom_bins_{}_shift_{}.txt'.format(bin_size, shift))

    bin_size = convert(bin_size)
    shift = convert(shift)

    df = pd.read_csv(chrom_sizes_file, sep="\t", header=None, names=["chrm", "size"])
    df['size'] = df['size'].apply(lambda x: calc_bin(x, bin_size, shift))

    count = 0

    for i in df.index:
        size = df.at[i, 'size']
        df.at[i, 'size'] = count
        count = count + size + 1

    df.to_csv(output, header=None, index=None, sep=' ', mode='w')


# generate the bin ranges for each chromosome
def generate_ranges(chrom_sizes_file, metadata, bin_size, shift):
    if shift == '0':
        output = os.path.join(metadata, 'chrom_bins_range_{}.txt'.format(bin_size))
    else:
        output = os.path.join(metadata, 'chrom_bins_range_{}_shift_{}.txt'.format(bin_size, shift))

    bin_size = convert(bin_size)
    shift = convert(shift)

    df = pd.read_csv(chrom_sizes_file, sep="\t", header=None, names=["chrm", "size"])
    df['start'] = df['size'].apply(lambda x: calc_bin(x, bin_size, shift))

    sLength = len(df['chrm'])
    df['end'] = np.zeros(sLength).astype(int)

    count = 0

    for i in df.index:
        size = df.at[i, 'start']
        df.at[i, 'start'] = count
        count = count + size + 1
        df.at[i, 'end'] = count - 1

    df.to_csv(output, header=None, index=None, sep=' ', mode='w')


# calculate total number of bins for each chromosome
def calc_bin(x, bin_size, shift):
    bin_num = int(x) // bin_size

    if shift == 0:
        return bin_num

    remainder = int(x) % bin_size

    if remainder > shift:
        bin_num += 1

    return bin_num


# convert the values like 1M, 500k to real integer values
def convert(val):
    if val == '0':
        return 0

    lookup = {'k': 1000, 'M': 1000000, 'B': 1000000000}
    unit = val[-1]
    try:
        number = int(val[:-1])
    except ValueError:
        sys.exit("Value Error.")
    if unit in lookup:
        return lookup[unit] * number
    return int(val)


def main():
    parser = ArgumentParser("Analyze Low Freq Interactions",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument("--data", required=True, help="Path for the root directory of the data files")
    parser.add_argument("--bin-size", required=True, help="Bin size, Eg: 1M, 500k")
    parser.add_argument("--sliding-window", default='0', help="Sliding windows for the bins")
    parser.add_argument("--config-file", required=True, help="Config file specifying the chromosome sizes")

    args = parser.parse_args()
    mat_file = args.data
    config_file = args.config_file
    bin_size = args.bin_size
    shift = args.sliding_window

    temp_dir = tempfile.TemporaryDirectory().name

    metadata = os.path.join(temp_dir, "metadata")

    # create directory if not exists
    if not os.path.exists(metadata):
        os.makedirs(metadata)

    # generate the starting bin indices using the chromosome sizes
    generate_bins(config_file, metadata, bin_size, shift)
    # generate the bin ranges for each chromosome
    generate_ranges(config_file, metadata, bin_size, shift)

    sum_matrix = np.genfromtxt(mat_file, delimiter=" ")

    if shift == '0':
        chrom_bin_range = "{}/chrom_bins_range_{}.txt".format(metadata, bin_size, shift)
    else:
        chrom_bin_range = "{}/chrom_bins_range_{}_shift_{}.txt".format(metadata, bin_size, shift)

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            splitLine = line.split()
            bin_range[splitLine[0]] = (int(splitLine[2]), int(splitLine[3]))

    rows = sum_matrix.shape[0]
    cols = sum_matrix.shape[1]

    value_list = []
    max = 0
    for i in range(0, rows):
        for j in range(0, cols):
            sum_matrix[i][j] = int(sum_matrix[i][j])
            for key, value in bin_range.items():
                if value[0] <= i <= value[1]:
                    if value[0] <= j <= value[1]:
                        sum_matrix[i][j] = 0
            if sum_matrix[i][j] != 0:
                if sum_matrix[i][j] > max:
                    max = sum_matrix[i][j]
                value_list.append(sum_matrix[i][j])

    print(max)
    plt.hist(value_list, bins=int(max))
    plt.show()


if __name__ == '__main__':
    main()
