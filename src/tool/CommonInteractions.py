import numpy as np


def main():
    matrix1 = np.genfromtxt("output/1CDX-SumMatrix/1CDX3/sum_matrix_500k.txt", delimiter=" ")
    matrix2 = np.genfromtxt("output/1CDX-SumMatrix/1CDX4/sum_matrix_500k.txt", delimiter=" ")

    rows = matrix1.shape[0]
    cols = matrix1.shape[1]

    chrom_bin_range = "metadata/chrom_bins_range_500k.txt"

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            splitLine = line.split()
            bin_range[splitLine[0]] = (int(splitLine[2]), int(splitLine[3]))

    count_inter = 0
    common_inter = 0
    matrix1_inter = 0
    matrix2_inter = 0
    both_zero = 0

    for i in range(0, rows):
        for j in range(0, cols):
            for key, value in bin_range.items():
                if value[0] <= i <= value[1]:
                    if not value[0] <= j <= value[1]:
                        count_inter += 1

                        if matrix1[i][j] == 0 and matrix2[i][j] == 0:
                            both_zero += 1

                        if matrix1[i][j] != 0:
                            matrix1_inter += 1
                            if matrix2[i][j] != 0:
                                common_inter += 1

                        if matrix2[i][j] != 0:
                            matrix2_inter += 1

    print("Total: " + str(count_inter))
    print("Common: " + str(common_inter))
    print("1CDX3: " + str(matrix1_inter))
    print("1CDX4: " + str(matrix2_inter))
    print("Both Zero: " + str(both_zero))


def main2():
    file1 = "output/tool-output-v2/pronucleus_female_500k_min.txt"
    file2 = "output/tool-output-v2/pronucleus_male_500k_min.txt"

    file1_entry = []
    with open(file1) as f:
        next(f)
        for line in f:
            split_line = line.split()
            file1_entry.append((int(split_line[0]), int(split_line[1])))

    file2_entry = []
    with open(file2) as f:
        next(f)
        for line in f:
            split_line = line.split()
            file2_entry.append((int(split_line[0]), int(split_line[1])))

    count_common = 0
    count_unique = 0

    for entry1 in file1_entry:
        found = False
        for entry2 in file2_entry:
            if entry1 == entry2:
                count_common += 1
                found = True
        if not found:
            count_unique += 1


    # for entry1 in file1_entry:
    #     found = False
    #     for entry2 in file2_entry:
    #         if (entry1[0] == entry2[0] and entry1[1] == entry2[1]) or (entry1[1] == entry2[0] and entry1[0] == entry2[1]):
    #             count_common += 1
    #             found = True
    #     if not found:
    #         count_unique += 1

    print("Total 1: "+str(len(file1_entry)))
    print("Total 2: "+str(len(file2_entry)))
    print("Common: "+ str(count_common))
    print("Unique 1: "+str(count_unique))

    print("{},{},{}".format(count_common,len(file1_entry)-count_common,len(file2_entry)-count_common))


if __name__ == '__main__':
    main2()

