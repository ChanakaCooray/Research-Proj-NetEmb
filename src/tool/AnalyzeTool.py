import os
from concurrent.futures import ThreadPoolExecutor
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys
import pandas as pd
import numpy as np


def process_chrom_file(filename, data_dir, chrom_bin, bin_size, output_dir, shift):
    print("Processing {}...".format(filename))
    with open(os.path.join(data_dir, filename)) as f:
        next(f)
        chr_entry = {}
        for line in f:
            split_line = line.split()

            start_index1 = int(chrom_bin[split_line[0]])
            start_index2 = int(chrom_bin[split_line[2]])

            coord1 = int(split_line[1])
            coord2 = int(split_line[3])

            # figure out the correct bin
            bin1 = find_bin(coord1, start_index1, bin_size, shift)
            bin2 = find_bin(coord2, start_index2, bin_size, shift)

            # check if the entry is already in the map, if it is increase the count number
            if (bin1, bin2) in chr_entry:
                chr_entry[(bin1, bin2)] = chr_entry[(bin1, bin2)] + int(split_line[4])
            elif (bin2, bin1) in chr_entry:
                chr_entry[(bin2, bin1)] = chr_entry[(bin2, bin1)] + int(split_line[4])
            else:
                chr_entry[(bin1, bin2)] = int(split_line[4])

        output_file = os.path.join(output_dir, filename)
        out = open(output_file, "w")
        for key, val in chr_entry.items():
            out.write("{} {} {}\n".format(key[0], key[1], val))
        out.close()


def process_data_files(data_dir, final_output_dir, metadata, shift, bin_size):
    # chromosome bins metadata
    if shift == 0:
        chrom_bin_file = "{}/chrom_bins_{}.txt".format(metadata, bin_size)
    else:
        chrom_bin_file = "{}/chrom_bins_{}_shift_{}.txt".format(metadata, bin_size, shift)

    output_edge_dir = os.path.join(final_output_dir, "temp", "edge_file")

    bin_size = convert(bin_size)

    # create output directory if not exists
    if not os.path.exists(output_edge_dir):
        os.makedirs(output_edge_dir)

    # store bin indexes
    chrom_bin = {}
    with open(chrom_bin_file) as f:
        for line in f:
            split_line = line.split()
            chrom_bin[split_line[0]] = split_line[1]

    print("Processing Files...")
    executor = ThreadPoolExecutor(40)
    for filename in os.listdir(data_dir):
        # process_chrom_file(filename, data_dir, chrom_bin, bin_size, output_edge_dir, shift)
        executor.submit(process_chrom_file, filename, data_dir, chrom_bin, bin_size, output_edge_dir, shift)


def find_bin(coord, start_index, bin_size, shift):
    bin_num = coord // bin_size
    remainder = coord % bin_size

    bin_num += start_index

    if shift == 0:
        return bin_num

    if remainder >= shift:
        bin_num += 1

    return bin_num


def generate_bins(chrom_sizes_file, metadata, bin_size, shift):
    # input = 'metadata/chrom_sizes.txt'
    # output = 'metadata/chrom_bins_{}_shift_{}.txt'.format(bin_size, shift)

    if shift == 0:
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


def generate_ranges(chrom_sizes_file, metadata, bin_size, shift):
    # input = 'metadata/chrom_sizes.txt'

    if shift == 0:
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


def calc_bin(x, bin_size, shift):
    bin_num = int(x) // bin_size

    if shift == 0:
        return bin_num

    remainder = int(x) % bin_size

    if remainder > shift:
        bin_num += 1

    return bin_num


def convert(val):
    if val == '0' or val == 0:
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
    parser = ArgumentParser("GenerateEdges",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument("--data", required=True, help="Path for the root directory of the data files")
    parser.add_argument("--output", required=True, help="output directory")
    parser.add_argument("--bin-size", required=True)
    parser.add_argument("--sliding-window", default='0')
    parser.add_argument("--config-file", required=True)

    args = parser.parse_args()
    data_dir = args.data
    final_output_dir = args.output
    metadata = os.path.join(final_output_dir, "temp", "metadata")
    shift = convert(args.sliding_window)
    config_file = args.config_file
    bin_size = args.bin_size

    # create output directory if not exists
    if not os.path.exists(metadata):
        os.makedirs(metadata)

    generate_bins(config_file, metadata, bin_size, shift)
    generate_ranges(config_file, metadata, bin_size, shift)
    process_data_files(data_dir, final_output_dir, metadata, shift, bin_size)


if __name__ == '__main__':
    main()
