import os


def main():
    analyze_dir = "analyze-1M"
    analyze_cat = "1CDX1"

    # chromBin_file = "metadata/chrom_bins_1M.txt"
    chrom_bin_range = "metadata/chrom_bins_range_1M.txt"
    # chrom_bin_range = "metadata/chrom_bins_range_500k.txt"

    bin_range = {}
    with open(chrom_bin_range) as f:
        for line in f:
            splitLine = line.split()
            bin_range[splitLine[0]] = (int(splitLine[2]), int(splitLine[3]))

    sum_max = 0
    for filename in os.listdir(analyze_dir):
        if not filename.startswith(analyze_cat):
            continue

        count = 0
        with open(os.path.join(analyze_dir, filename)) as f:
            for line in f:
                splitLine = line.split()
                edge1 = int(splitLine[0])
                edge2 = int(splitLine[1])

                intra_chrm = False
                for key, value in bin_range.items():
                    if value[0] <= edge1 <= value[1]:
                        if value[0] <= edge2 <= value[1]:
                            intra_chrm = True
                            break

                if not intra_chrm:
                    count += 1
        if count > sum_max:
            sum_max = count

    print(sum_max)


if __name__ == '__main__':
    main()
