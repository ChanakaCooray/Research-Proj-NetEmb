import os
from concurrent.futures import ThreadPoolExecutor
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys
import matplotlib.pyplot as plt
import random
from statistics import mean
from statistics import stdev
from collections import Counter


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


# only used k bins
def main():
    parser = ArgumentParser("EdgeDetector",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    # parser.add_argument("--cell")
    parser.add_argument("--type")
    parser.add_argument("--k")

    args = parser.parse_args()
    # cell = args.cell
    type = args.type
    k = int(args.k)

    chrom_bin_range = "metadata/chrom_bins_range_500k.txt"

    # search for last index
    with open(chrom_bin_range) as f:
        for last in f: pass
        n = int(last.split()[3]) + 1

    # print(n)

    chrom_bin_file = "metadata/chrom_bins_500k.txt"

    # store bin indexes
    chrom_bin = {}
    with open(chrom_bin_file) as f:
        for line in f:
            split_line = line.split()
            chrom_bin["chr" + split_line[0]] = split_line[1]

    mESC_file = "mESC/mESC.{}.txt".format(type)
    # mESC_file = "mESC/mESC.enhancer.txt"
    # mESC_file = "mESC/mESC.h3k4me3.peak.txt"
    # mESC_file = "mESC/mESC.h3k27ac.peak.txt"
    # mESC_file = "mESC/mESC.polII.peak.txt"

    # sig_bin_file = "output/tool-output-v2/{}_500k_max.txt".format(cell)
    # sig_bin_file = "output/tool-output-v2/1CDX2_500k_max.txt"
    # sig_bin_file = "output/tool-output-v2/1CDX3_500k_max.txt"
    # sig_bin_file = "output/tool-output-v2/1CDX4_500k_max.txt"

    for y in range(1, 5):
        sig_bin_file = "output/tool-output-v2/1CDX{}_500k_max.txt".format(y)

        # print(mESC_file)
        # print(sig_bin_file)

        mESC_bin_map = {}

        for i in range(0, n + 1):
            mESC_bin_map[i] = 0

        with open(mESC_file) as f:
            for line in f:
                split_line = line.split()
                bin_n = find_bin(int(split_line[1]), int(chrom_bin[split_line[0]]), 500000, 0)

                mESC_bin_map[bin_n] += 1

        # bin_set = set()
        bin_list = []

        with open(sig_bin_file) as f:
            next(f)
            for line in f:
                split_line = line.split()
                # bin_set.add(int(split_line[0]))
                # bin_set.add(int(split_line[1]))

                bin_list.append(int(split_line[0]))
                bin_list.append(int(split_line[1]))

        bin_map = Counter(bin_list)
        ordered_list_tuples = sorted(bin_map.items(), reverse=True, key=lambda x: x[1])

        # k = 25
        bin_set = []
        frequency_list = []
        cou = 0
        for elem in ordered_list_tuples:
            bin_set.append(elem[0])
            frequency_list.append(elem[1])
            cou += 1

            if cou == k:
                break

        # print(bin_set)
        # print(len(bin_set))
        print(frequency_list)
        # sys.exit(0)

        # sig_bin_values = []

        n_value = 0
        for j in bin_set:
            # sig_bin_values.append(mESC_bin_map[j])
            n_value += mESC_bin_map[j]

        # n_value = sum(sig_bin_values)

        # print("Unique: " + str(len(bin_set)))
        # print(n_value)

        all_bin_list = list(mESC_bin_map.keys())

        random_value_list = []
        for z in range(0, 10000):
            random_list = random.sample(all_bin_list, len(bin_set))

            random_value = 0
            for j in random_list:
                random_value += mESC_bin_map[j]

            random_value_list.append(random_value)

        # print(mean(random_value_list))
        # print(stdev(random_value_list))

        z_score = (n_value - mean(random_value_list)) / stdev(random_value_list)

        print(z_score)

        # all_bin_values = list(mESC_bin_map.values())

        # fig = plt.figure(figsize=(8,6))
        # plt.boxplot([x for x in [sig_bin_values,all_bin_values]], 0, 'rs', 1)
        # plt.xticks([y+1 for y in range(len([sig_bin_values,all_bin_values]))], ['Significant Bins','All Bins'])
        # plt.xlabel('measurement x')
        # t = plt.title('Box plot')
        # plt.show()


# old method using all of the bins
def main2():
    parser = ArgumentParser("EdgeDetector",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument("--cell")
    parser.add_argument("--type")

    args = parser.parse_args()
    cell = args.cell
    type = args.type

    chrom_bin_range = "metadata/chrom_bins_range_500k.txt"

    # search for last index
    with open(chrom_bin_range) as f:
        for last in f: pass
        n = int(last.split()[3]) + 1

    # print(n)

    chrom_bin_file = "metadata/chrom_bins_500k.txt"

    # store bin indexes
    chrom_bin = {}
    with open(chrom_bin_file) as f:
        for line in f:
            split_line = line.split()
            chrom_bin["chr" + split_line[0]] = split_line[1]

    mESC_file = "mESC/mESC.{}.txt".format(type)
    # mESC_file = "mESC/mESC.enhancer.txt"
    # mESC_file = "mESC/mESC.h3k4me3.peak.txt"
    # mESC_file = "mESC/mESC.h3k27ac.peak.txt"
    # mESC_file = "mESC/mESC.polII.peak.txt"

    sig_bin_file = "output/tool-output-v2/{}_500k_max.txt".format(cell)
    # sig_bin_file = "output/tool-output-v2/1CDX2_500k_max.txt"
    # sig_bin_file = "output/tool-output-v2/1CDX3_500k_max.txt"
    # sig_bin_file = "output/tool-output-v2/1CDX4_500k_max.txt"

    print(mESC_file)
    print(sig_bin_file)

    mESC_bin_map = {}

    for i in range(0, n + 1):
        mESC_bin_map[i] = 0

    with open(mESC_file) as f:
        for line in f:
            split_line = line.split()
            bin_n = find_bin(int(split_line[1]), int(chrom_bin[split_line[0]]), 500000, 0)

            mESC_bin_map[bin_n] += 1

    bin_set = set()

    with open(sig_bin_file) as f:
        next(f)
        for line in f:
            split_line = line.split()
            bin_set.add(int(split_line[0]))
            bin_set.add(int(split_line[1]))

    print(len(bin_set))

    # sys.exit(0)

    # sig_bin_values = []

    n_value = 0
    for j in bin_set:
        # sig_bin_values.append(mESC_bin_map[j])
        n_value += mESC_bin_map[j]

    # n_value = sum(sig_bin_values)

    print("Unique: " + str(len(bin_set)))
    print(n_value)

    all_bin_list = list(mESC_bin_map.keys())

    random_value_list = []
    for z in range(0, 10000):
        random_list = random.sample(all_bin_list, len(bin_set))

        random_value = 0
        for j in random_list:
            random_value += mESC_bin_map[j]

        random_value_list.append(random_value)

    print(mean(random_value_list))
    print(stdev(random_value_list))

    z_score = (n_value - mean(random_value_list)) / stdev(random_value_list)

    print(z_score)

    # all_bin_values = list(mESC_bin_map.values())

    # fig = plt.figure(figsize=(8,6))
    # plt.boxplot([x for x in [sig_bin_values,all_bin_values]], 0, 'rs', 1)
    # plt.xticks([y+1 for y in range(len([sig_bin_values,all_bin_values]))], ['Significant Bins','All Bins'])
    # plt.xlabel('measurement x')
    # t = plt.title('Box plot')
    # plt.show()


# only used k bins
def main3():
    parser = ArgumentParser("EdgeDetector",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    # parser.add_argument("--cell")
    parser.add_argument("--type")
    # parser.add_argument("--k")

    args = parser.parse_args()
    # cell = args.cell
    type = args.type
    # k = int(args.k)

    chrom_bin_range = "metadata/chrom_bins_range_500k.txt"

    # search for last index
    with open(chrom_bin_range) as f:
        for last in f: pass
        n = int(last.split()[3]) + 1

    # print(n)

    chrom_bin_file = "metadata/chrom_bins_500k.txt"

    # store bin indexes
    chrom_bin = {}
    with open(chrom_bin_file) as f:
        for line in f:
            split_line = line.split()
            chrom_bin["chr" + split_line[0]] = split_line[1]

    mESC_file = "mESC/mESC.{}.txt".format(type)
    # mESC_file = "mESC/mESC.enhancer.txt"
    # mESC_file = "mESC/mESC.h3k4me3.peak.txt"
    # mESC_file = "mESC/mESC.h3k27ac.peak.txt"
    # mESC_file = "mESC/mESC.polII.peak.txt"

    # sig_bin_file = "output/tool-output-v2/{}_500k_max.txt".format(cell)
    # sig_bin_file = "output/tool-output-v2/1CDX2_500k_max.txt"
    # sig_bin_file = "output/tool-output-v2/1CDX3_500k_max.txt"
    # sig_bin_file = "output/tool-output-v2/1CDX4_500k_max.txt"

    for y in range(1, 5):
        # print("1CDX{}".format(y))
        sig_bin_file = "output/tool-output-v2/1CDX{}_500k_max.txt".format(y)

        # print(mESC_file)
        # print(sig_bin_file)

        mESC_bin_map = {}

        for i in range(0, n + 1):
            mESC_bin_map[i] = 0

        with open(mESC_file) as f:
            for line in f:
                split_line = line.split()
                bin_n = find_bin(int(split_line[1]), int(chrom_bin[split_line[0]]), 500000, 0)

                mESC_bin_map[bin_n] += 1

        # bin_set = set()
        bin_list = []

        with open(sig_bin_file) as f:
            next(f)
            for line in f:
                split_line = line.split()
                # bin_set.add(int(split_line[0]))
                # bin_set.add(int(split_line[1]))

                bin_list.append(int(split_line[0]))
                bin_list.append(int(split_line[1]))

        bin_map = Counter(bin_list)
        ordered_list_tuples = sorted(bin_map.items(), reverse=True, key=lambda x: x[1])

        # k = 25
        bin_set = []
        frequency_list = []
        # cou = 0
        for elem in ordered_list_tuples:
            # print(str(elem[0]) + "a" + str(elem[1]))
            if elem[1] < 3:
                break
            bin_set.append(elem[0])
            frequency_list.append(elem[1])
            # cou += 1

            # if cou == k:
            #     break

        # print(bin_set)
        # print(len(bin_set))
        # print(frequency_list)
        # print(len(ordered_list_tuples))
        # sys.exit(0)

        # sig_bin_values = []

        n_value = 0
        for j in bin_set:
            # sig_bin_values.append(mESC_bin_map[j])
            n_value += mESC_bin_map[j]

        # n_value = sum(sig_bin_values)

        # print("Unique: " + str(len(bin_set)))
        # print(n_value)

        all_bin_list = list(mESC_bin_map.keys())

        random_value_list = []
        for z in range(0, 50000):
            random_list = random.sample(all_bin_list, len(bin_set))

            random_value = 0
            for j in random_list:
                random_value += mESC_bin_map[j]

            random_value_list.append(random_value)

        # print(mean(random_value_list))
        # print(stdev(random_value_list))

        z_score = (n_value - mean(random_value_list)) / stdev(random_value_list)

        print(z_score)

        # all_bin_values = list(mESC_bin_map.values())

        # fig = plt.figure(figsize=(8,6))
        # plt.boxplot([x for x in [sig_bin_values,all_bin_values]], 0, 'rs', 1)
        # plt.xticks([y+1 for y in range(len([sig_bin_values,all_bin_values]))], ['Significant Bins','All Bins'])
        # plt.xlabel('measurement x')
        # t = plt.title('Box plot')
        # plt.show()


if __name__ == '__main__':
    main3()
