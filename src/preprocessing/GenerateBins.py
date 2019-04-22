import pandas as pd
import numpy as np
import sys


def generate_bins(binSize_Str):
    input = 'metadata/chrom_sizes.txt'
    output = 'metadata/chrom_bins_{}.txt'.format(binSize_Str)

    binSize = convert(binSize_Str)

    df = pd.read_csv(input, sep="\t", header=None, names=["chrm", "size"])
    df['size'] = df['size'].apply(lambda x: int(x) // binSize)

    count = 0

    for i in df.index:
        size = df.at[i, 'size']
        df.at[i, 'size'] = count
        count = count + size + 1

    df.to_csv(output, header=None, index=None, sep=' ', mode='w')


def generate_ranges(binSize_Str):
    input = 'metadata/chrom_sizes.txt'
    output = 'metadata/chrom_bins_range_{}.txt'.format(binSize_Str)

    binSize = convert(binSize_Str)

    df = pd.read_csv(input, sep="\t", header=None, names=["chrm", "size"])
    df['start'] = df['size'].apply(lambda x: int(x) // binSize)

    sLength = len(df['chrm'])
    df['end'] = np.zeros(sLength).astype(int)

    count = 0

    for i in df.index:
        size = df.at[i, 'start']
        df.at[i, 'start'] = count
        count = count + size + 1
        df.at[i, 'end'] = count - 1

    df.to_csv(output, header=None, index=None, sep=' ', mode='w')


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


if __name__ == '__main__':
    generate_bins("40k")
    generate_ranges("40k")
