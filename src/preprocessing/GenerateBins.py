import pandas as pd
import numpy as np


def generate_bins():
    binSize = 500000
    input = 'metadata/chrom_sizes.txt'
    output = 'metadata/chrom_bins_500k.txt'

    df = pd.read_csv(input, sep="\t", header=None, names=["chrm", "size"])
    df['size'] = df['size'].apply(lambda x: int(x) // binSize)

    count = 0

    for i in df.index:
        size = df.at[i, 'size']
        df.at[i, 'size'] = count
        count = count + size + 1

    df.to_csv(output, header=None, index=None, sep=' ', mode='w')


def generate_ranges():
    binSize = 500000
    input = 'metadata/chrom_sizes.txt'
    output = 'metadata/chrom_bins_range_500k.txt'

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

    # print(df)


if __name__ == '__main__':
    generate_ranges()
