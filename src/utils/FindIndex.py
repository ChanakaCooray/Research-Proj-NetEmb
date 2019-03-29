import sys
import os


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


def find_index(filename, metadata, bin_size, shift, output_file):
    if shift == '0':
        chrom_bin_range = "{}/chrom_bins_range_{}.txt".format(metadata, bin_size)
    else:
        chrom_bin_range = "{}/chrom_bins_range_{}_shift_{}.txt".format(metadata, bin_size, shift)

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            split_line = line.split()
            bin_range[split_line[0]] = (int(split_line[2]), int(split_line[3]), int(split_line[1]))

    bin_size = convert(bin_size)
    shift = convert(shift)

    out = open(output_file, "w")

    with open(os.path.join(filename)) as f:
        next(f)
        for line in f:
            split_line = line.split()
            bin1 = int(split_line[0])
            bin2 = int(split_line[1])
            count = int(split_line[2])

            chrm1 = ""
            chrm2 = ""
            start_index1 = 0
            start_index2 = 0
            chrm1_size = 0
            chrm2_size = 0
            for key, value in bin_range.items():
                if value[0] <= bin1 <= value[1]:
                    chrm1 = key
                    start_index1 = value[0]
                    chrm1_size = value[2]
                if value[0] <= bin2 <= value[1]:
                    chrm2 = key
                    start_index2 = value[0]
                    chrm2_size = value[2]

            bin1_num = bin1 - start_index1
            bin2_num = bin2 - start_index2

            if shift == 0:
                bin1_start_index = bin1_num * bin_size
                bin1_end_index = bin1_start_index + bin_size - 1
            elif bin1_num != 0:
                bin1_start_index = (bin1_num - 1) * bin_size + shift
                bin1_end_index = bin1_start_index + bin_size - 1
            else:
                bin1_start_index = 0
                bin1_end_index = bin1_start_index + shift - 1

            if bin1_end_index >= chrm1_size:
                bin1_end_index = chrm1_size - 1

            if shift == 0:
                bin2_start_index = bin2_num * bin_size
                bin2_end_index = bin2_start_index + bin_size - 1
            elif bin2_num != 0:
                bin2_start_index = (bin2_num - 1) * bin_size + shift
                bin2_end_index = bin2_start_index + bin_size - 1
            else:
                bin2_start_index = 0
                bin2_end_index = bin2_start_index + shift - 1

            if bin2_end_index >= chrm2_size:
                bin2_end_index = chrm2_size - 1

            # print("{} {} {} {} {} {}".format(bin1, bin2, bin1_start_index, bin1_end_index, bin2_start_index,
            #                                  bin2_end_index))

            out.write("{} {} {} {} {} {} {}\n".format("mm" + chrm1, bin1_start_index, bin1_end_index, "mm" + chrm2,
                                                      bin2_start_index, bin2_end_index, "color=black"))

    out.close()


if __name__ == '__main__':
    for i in range(1, 5):
        find_index("output/tool-output/v2/500k/1CDX{}/output_sum_matrix_500k.txt".format(i), "metadata", "500k", "0",
                   "circos-output/circos-input/input_1CDX{}_500k.txt".format(i))
        find_index("output/tool-output/v2/1M/1CDX{}/output_sum_matrix_1M.txt".format(i), "metadata", "1M", "0",
                   "circos-output/circos-input/input_1CDX{}_1M.txt".format(i))

    for i in range(1, 5):
        for j in range(1, 5):
            find_index("output/tool-output/v2/500k/1CDX{}/output_sum_matrix_500k_shift_{}k.txt".format(i, j * 100),
                       "metadata", "500k",
                       "{}k".format(j * 100),
                       "circos-output/circos-input/input_1CDX{}_500k_shift_{}k.txt".format(i, j * 100))
