import pandas as pd
import numpy as np
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def generate_bins(bin_size, shift):
    input = 'metadata/chrom_sizes.txt'
    output = 'metadata/chrom_bins_{}_shift_{}.txt'.format(bin_size, shift)

    bin_size = convert(bin_size)
    shift = convert(shift)

    df = pd.read_csv(input, sep="\t", header=None, names=["chrm", "size"])
    df['size'] = df['size'].apply(lambda x: calc_bin(x, bin_size, shift))

    count = 0

    for i in df.index:
        size = df.at[i, 'size']
        df.at[i, 'size'] = count
        count = count + size + 1

    df.to_csv(output, header=None, index=None, sep=' ', mode='w')


def calc_bin(x, bin_size, shift):
    bin_num = int(x) // bin_size

    if shift == 0:
        return bin_num

    remainder = int(x) % bin_size

    if remainder > shift:
        bin_num += 1

    return bin_num


def generate_ranges(bin_size, shift):
    input = 'metadata/chrom_sizes.txt'
    output = 'metadata/chrom_bins_range_{}_shift_{}.txt'.format(bin_size, shift)

    bin_size = convert(bin_size)
    shift = convert(shift)

    df = pd.read_csv(input, sep="\t", header=None, names=["chrm", "size"])
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


def convert(val):
    lookup = {'k': 1000, 'M': 1000000, 'B': 1000000000}
    unit = val[-1]
    try:
        number = int(val[:-1])
    except ValueError:
        sys.exit("Value Error.")
    if unit in lookup:
        return lookup[unit] * number
    return int(val)


if __name__ == '__main__':
    parser = ArgumentParser("GenerateBins",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument("--bin-size", required=True)
    parser.add_argument("--shift", required=True)

    args = parser.parse_args()
    generate_bins(args.bin_size, args.shift)
    # generate_ranges(args.bin_size, args.shift)
