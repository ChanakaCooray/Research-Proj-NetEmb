import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def main():
    parser = ArgumentParser("GetMaxSumMatrix",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument("--analyze-dir", required=True)
    parser.add_argument("--cell-type", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--bin-size", required=True)

    args = parser.parse_args()
    analyze_dir = args.analyze_dir
    analyze_cat = args.cell_type
    metadata = args.metadata
    bin_size = args.bin_size

    # chromBin_file = "metadata/chrom_bins_1M.txt"
    chrom_bin_range = "{}/chrom_bins_range_{}.txt".format(metadata, bin_size)
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
