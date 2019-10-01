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


# process the data file to an edge list file of bins
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


# process all the data files concurrently
def generate_data_files(data_dir, output_edge_dir, metadata, shift, bin_size):
    # chromosome bins metadata
    if shift == '0':
        chrom_bin_file = "{}/chrom_bins_{}.txt".format(metadata, bin_size)
    else:
        chrom_bin_file = "{}/chrom_bins_{}_shift_{}.txt".format(metadata, bin_size, shift)

    bin_size = convert(bin_size)
    shift = convert(shift)

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

    executor.shutdown(wait=True)


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


# generate the maximum sum of all the cells
def get_max_sum(data_dir, metadata, bin_size, shift):
    if shift == '0':
        chrom_bin_range = "{}/chrom_bins_range_{}.txt".format(metadata, bin_size)
    else:
        chrom_bin_range = "{}/chrom_bins_range_{}_shift_{}.txt".format(metadata, bin_size, shift)

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            split_line = line.split()
            bin_range[split_line[0]] = (int(split_line[2]), int(split_line[3]))

    max_sum = None
    for filename in os.listdir(data_dir):

        count = 0
        with open(os.path.join(data_dir, filename)) as f:
            for line in f:
                split_line = line.split()
                edge1 = int(split_line[0])
                edge2 = int(split_line[1])

                intra_chrom = False
                for key, value in bin_range.items():
                    if value[0] <= edge1 <= value[1]:
                        if value[0] <= edge2 <= value[1]:
                            intra_chrom = True
                            break

                if not intra_chrom:
                    count += 1

        if not max_sum:
            max_sum = count
        elif count > max_sum:
            max_sum = count

    return max_sum


# generate the minimum sum of all the cells
def get_min_sum(data_dir, metadata, bin_size, shift):
    if shift == '0':
        chrom_bin_range = "{}/chrom_bins_range_{}.txt".format(metadata, bin_size)
    else:
        chrom_bin_range = "{}/chrom_bins_range_{}_shift_{}.txt".format(metadata, bin_size, shift)

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            split_line = line.split()
            bin_range[split_line[0]] = (int(split_line[2]), int(split_line[3]))

    value_list = []
    for filename in os.listdir(data_dir):

        count = 0
        with open(os.path.join(data_dir, filename)) as f:
            for line in f:
                split_line = line.split()
                edge1 = int(split_line[0])
                edge2 = int(split_line[1])

                intra_chrom = False
                for key, value in bin_range.items():
                    if value[0] <= edge1 <= value[1]:
                        if value[0] <= edge2 <= value[1]:
                            intra_chrom = True
                            break

                if not intra_chrom:
                    count += 1

        value_list.append(count)

    np_values = np.array(value_list)
    percent_5 = np.percentile(np_values, 5)

    return percent_5


# generate the mean of sum values of all the cells
def get_mean(data_dir, metadata, bin_size, shift, number_of_cells):
    if shift == '0':
        chrom_bin_range = "{}/chrom_bins_range_{}.txt".format(metadata, bin_size)
    else:
        chrom_bin_range = "{}/chrom_bins_range_{}_shift_{}.txt".format(metadata, bin_size, shift)

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            split_line = line.split()
            bin_range[split_line[0]] = (int(split_line[2]), int(split_line[3]))

    sum = 0
    for filename in os.listdir(data_dir):
        count = 0
        with open(os.path.join(data_dir, filename)) as f:
            for line in f:
                split_line = line.split()
                edge1 = int(split_line[0])
                edge2 = int(split_line[1])

                intra_chrom = False
                for key, value in bin_range.items():
                    if value[0] <= edge1 <= value[1]:
                        if value[0] <= edge2 <= value[1]:
                            intra_chrom = True
                            break

                if not intra_chrom:
                    count += 1
        sum += count

    return sum / number_of_cells


# generate the summation matrix
def generate_sum_matrix(data_dir, metadata, bin_size, shift, output_dir):
    # create output directory if not exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if shift == '0':
        output_file = os.path.join(output_dir, "sum_matrix_{}.txt".format(bin_size))
        chrom_bin_range = "{}/chrom_bins_range_{}.txt".format(metadata, bin_size)
    else:
        output_file = os.path.join(output_dir, "sum_matrix_{}_shift_{}.txt".format(bin_size, shift))
        chrom_bin_range = "{}/chrom_bins_range_{}_shift_{}.txt".format(metadata, bin_size, shift)

    # search for last index
    with open(chrom_bin_range) as f:
        for last in f: pass
        n = int(last.split()[3]) + 1

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            splitLine = line.split()
            bin_range[splitLine[0]] = (int(splitLine[2]), int(splitLine[3]))

    sum_matrix = np.zeros((n, n), dtype=np.int)

    # non_zero_values = []

    for filename in os.listdir(data_dir):

        data = {}
        with open(os.path.join(data_dir, filename)) as f:
            for line in f:
                splitLine = line.split()
                edge1 = int(splitLine[0])
                edge2 = int(splitLine[1])
                val = int(splitLine[2])

                if 2 <= val <= 100:
                    val = 1
                else:
                    val = 0

                for key, value in bin_range.items():
                    if value[0] <= edge1 <= value[1]:
                        if value[0] <= edge2 <= value[1]:
                            val = 1
                            break

                if edge1 <= edge2:
                    data[(edge1, edge2)] = val
                else:
                    data[(edge2, edge1)] = val

        matrix = np.zeros((n, n), dtype=np.int)

        for key, val in data.items():
            matrix[key[0]][key[1]] = val
            matrix[key[1]][key[0]] = val

        # rows = matrix.shape[0]
        # cols = matrix.shape[1]
        #
        # count_non_zero = 0
        # count_total = 0

        # matrix = np.triu(matrix)

        # for i in range(0, rows):
        #     for j in range(0, cols):
        #         for key, value in bin_range.items():
        #             if value[0] <= i <= value[1]:
        #                 if not value[0] <= j <= value[1]:
        #                     count_total += 1
        #                     if matrix[i][j] != 0:
        #                         count_non_zero += 1
        #                     break

        # print(filename)
        # print("total: " + str(count_total))
        # print(str(count_non_zero))

        # non_zero_values.append(count_non_zero)

        # np.fill_diagonal(matrix, 0)

        sum_matrix = sum_matrix + matrix

    # np_non_zero_values = np.array(non_zero_values)
    # percent_5 = np.percentile(np_non_zero_values, 5)
    # percent_25 = np.percentile(np_non_zero_values, 25)
    # percent_50 = np.percentile(np_non_zero_values, 50)
    # percent_75 = np.percentile(np_non_zero_values, 75)
    # percent_95 = np.percentile(np_non_zero_values, 95)
    #
    # print("percentile 5: " + str(percent_5))
    # print("percentile 25: " + str(percent_25))
    # print("percentile 50: " + str(percent_50))
    # print("percentile 75: " + str(percent_75))
    # print("percentile 95: " + str(percent_95))

    # output_file2 = os.path.join(
    #     "/gpfs1/scratch/chanaka.cooray/jobs/work/project/Research-Proj-NetEmb/temp/1CDX4/sum_matrix_{}.txt".format(
    #         bin_size))
    # np.savetxt(output_file2, sum_matrix, fmt='%d', delimiter=' ', newline='\n')

    np.savetxt(output_file, sum_matrix, fmt='%d', delimiter=' ', newline='\n')

    return output_file


def write_zero_bin_output(zero_bin_list, output_file, shift, bin_size, metadata):
    # create directories for output file if not exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    if shift == '0':
        chrom_bin_range = "{}/chrom_bins_range_{}.txt".format(metadata, bin_size)
    else:
        chrom_bin_range = "{}/chrom_bins_range_{}_shift_{}.txt".format(metadata, bin_size, shift)

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            split_line = line.split()
            bin_range[split_line[0]] = (int(split_line[1]), int(split_line[2]))

    bin_size = convert(bin_size)
    shift = convert(shift)

    out = open(output_file, "w")
    out.write("{} {} {} {}\n".format("Chromosome", "Bin", "Start-Index", "End-Index"))
    for key, val in zero_bin_list.items():
        bin_n = key
        chrm_n = val

        chrm_info = bin_range[chrm_n]
        bin_count = bin_n - chrm_info[1]

        if shift == 0:
            start_index = bin_count * bin_size
            end_index = start_index + bin_size - 1
        elif bin_count != 0:
            start_index = (bin_count - 1) * bin_size + shift
            end_index = start_index + bin_size - 1
        else:
            start_index = 0
            end_index = start_index + shift - 1

        if end_index >= chrm_info[0]:
            end_index = chrm_info[0] - 1

        out.write("{} {} {} {}\n".format(val, key, start_index, end_index))
    out.close()


# generate the final analysis
def generate_analyzed_output(sum_matrix, sum_value, metadata, bin_size, shift, number_of_cells, output_file,
                             p_value_user, zero_bin_file):
    if shift == '0':
        chrom_bin_range = "{}/chrom_bins_range_{}.txt".format(metadata, bin_size, shift)
    else:
        chrom_bin_range = "{}/chrom_bins_range_{}_shift_{}.txt".format(metadata, bin_size, shift)

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            splitLine = line.split()
            bin_range[splitLine[0]] = (int(splitLine[2]), int(splitLine[3]))

    data = np.genfromtxt(sum_matrix, delimiter=" ")

    rows = data.shape[0]
    cols = data.shape[1]

    zero_bins = 0
    sum_by_row = data.sum(axis=0).astype(int)

    zero_bin_list = {}
    for i, val in enumerate(sum_by_row):
        if val == 0:
            zero_bins += 1
            for key, value in bin_range.items():
                if value[0] <= i <= value[1]:
                    zero_bin_list[i] = key

    df_zero_bin = pd.DataFrame(list(zero_bin_list.items()), columns=['bin', 'chrm'])

    if not zero_bin_file == 'null':
        write_zero_bin_output(zero_bin_list, zero_bin_file, shift, bin_size, metadata)

    df_zero_bin_count = df_zero_bin['chrm'].value_counts().reset_index()
    df_zero_bin_count.columns = ['chrm', 'count']
    df_zero_bin_count.set_index('chrm', inplace=True)

    # M = ( (N-X)^2 - (n1-x1)^2 - (n2-x2)^2 - ... - (n20-x20)^2 - (n21-x21)^2 )/2
    sum_reduction = 0  # stores summation of (n2-x2)^2 + ... + (n20-x20)^2 + (n21-x21)^2
    count_X = 0  # stores X
    count_N = 0  # stores N
    for key, value in bin_range.items():
        n = value[1] - value[0] + 1  # store the number of bins for the chromosome
        zero_bins_chrm = df_zero_bin_count.loc[key, 'count']  # stores zero bins for the chrm
        sum_reduction = sum_reduction + (n - zero_bins_chrm) ** 2  # cal (n2-x2)^2 and add that to sum
        count_X = count_X + zero_bins_chrm
        count_N = count_N + n

    count = 0
    # count_inter = 0
    for i in range(0, rows):
        for j in range(0, cols):
            for key, value in bin_range.items():
                if value[0] <= i <= value[1]:
                    if value[0] <= j <= value[1]:
                        data[i][j] = 0
                        count += 1
                    # else:
                    #     count_inter += 1

    # print("rows: " + str(rows))
    # print("cols: " + str(cols))
    # print("Intra: " + str(count))
    # print("All: " + str(rows * cols))
    # print("Inter: " + str(count_inter))

    # unique, counts = np.unique(data, return_counts=True)
    M = ((count_N - count_X) ** 2 - sum_reduction) / 2
    count1 = 0
    p_value_thresh = calculate_pval_thresh(M, sum_value)
    data = np.triu(data, k=-1)
    out = open(output_file, "w")
    out.write("{} {} {} {}\n".format("bin1", "bin2", "count", "p_value"))

    for i in range(0, rows):
        for j in range(0, cols):
            if data[i][j] != 0:
                value = int(data[i][j])
                p_value = calculate_f_sum(number_of_cells, value, p_value_thresh)
                if p_value <= threshold(M, p_value_user):
                    count1 += 1
                    out.write("{} {} {} {}\n".format(i, j, value, p_value))

    out.close()


def calculate_f(n, t, p_max):
    return nCr(n, t) * (p_max ** t) * ((1 - p_max) ** (n - t))


# calculate f sum using binomial function
def calculate_f_sum(n, t, p_max):
    sum_f = 1 - binom.cdf(t - 1, n, p_max)
    return sum_f


def calculate_pval_thresh(M, sum_value):
    p_max = sum_value / M
    return p_max


# helper function to calc nCr
def nCr(n, r):
    f = math.factorial
    return f(n) // f(r) // f(n - r)


def threshold(M, p_value_user):
    return p_value_user / M


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
    parser = ArgumentParser("Analyze tool",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument("--data", required=True, help="Path for the root directory of the data files")
    parser.add_argument("--output", required=True, help="Path of the output file")
    parser.add_argument("--bin-size", required=True, help="Bin size, Eg: 1M, 500k")
    parser.add_argument("--sliding-window", default='0', help="Sliding windows for the bins")
    parser.add_argument("--config-file", required=True, help="Config file specifying the chromosome sizes")
    parser.add_argument("--min", action='store_true', default=False)
    parser.add_argument("--mean", action='store_true', default=False)
    parser.add_argument("--max", action='store_true', default=False)
    parser.add_argument("--zero-bin-file", default='null',
                        help="If specified, the program outputs a file that lists the bins with no interactions.")
    parser.add_argument("--p-value", default='0.05', help="P Value")

    args = parser.parse_args()
    data_dir = args.data
    output_file = args.output
    shift = args.sliding_window
    config_file = args.config_file
    bin_size = args.bin_size
    p_value_user = float(args.p_value)
    zero_bin_file = args.zero_bin_file

    if args.min:
        level = "min"
    elif args.mean:
        level = "mean"
    else:
        level = "max"

    temp_dir = tempfile.TemporaryDirectory().name

    # temp directory to store meta data and temporary output files
    output_edge_dir = os.path.join(temp_dir, "edge_files")
    metadata = os.path.join(temp_dir, "metadata")

    # create output directory if not exists
    if not os.path.exists(output_edge_dir):
        os.makedirs(output_edge_dir)

    # create directory if not exists
    if not os.path.exists(metadata):
        os.makedirs(metadata)

    # create directories for output file if not exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # count the number of cells/files in the directory for the analysis calculations
    list_files = os.listdir(data_dir)
    number_of_cells = len(list_files)

    # generate the starting bin indices using the chromosome sizes
    generate_bins(config_file, metadata, bin_size, shift)
    # generate the bin ranges for each chromosome
    generate_ranges(config_file, metadata, bin_size, shift)
    # generate the edge list using bin indices
    generate_data_files(data_dir, output_edge_dir, metadata, shift, bin_size)
    # get the maximum sum using only inter-chromosome interactions

    if level == "max":
        sum_value = get_max_sum(output_edge_dir, metadata, bin_size, shift)
    elif level == "mean":
        sum_value = get_mean(output_edge_dir, metadata, bin_size, shift, number_of_cells)
    else:
        sum_value = get_min_sum(output_edge_dir, metadata, bin_size, shift)

    # generate the summation matrix for all the cells
    sum_matrix = generate_sum_matrix(output_edge_dir, metadata, bin_size, shift,
                                     os.path.join(temp_dir, "sum-matrix"))

    # print("Sum Value: {}".format(sum_value))

    # generate the output file with significant inter-chromosome interactions
    generate_analyzed_output(sum_matrix, sum_value, metadata, bin_size, shift, number_of_cells, output_file,
                             p_value_user, zero_bin_file)

    # clean up temporary files
    try:
        shutil.rmtree(os.path.join(temp_dir))
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))

    print("Output file created.")


if __name__ == '__main__':
    main()
