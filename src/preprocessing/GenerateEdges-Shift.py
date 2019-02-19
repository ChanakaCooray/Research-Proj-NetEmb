import os
from concurrent.futures import ThreadPoolExecutor
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys


def processFile(filename, dataDir, chromBin, binSize, outputDir, shift):
    print("Processing {}...".format(filename))
    with open(os.path.join(dataDir, filename)) as f:
        next(f)
        chrEntry = {}
        for line in f:
            splitLine = line.split()

            startIndex1 = int(chromBin[splitLine[0]])
            startIndex2 = int(chromBin[splitLine[2]])

            coord1 = int(splitLine[1])
            coord2 = int(splitLine[3])

            # figure out the correct bin
            bin1 = find_bin(coord1, startIndex1, binSize, shift)
            bin2 = find_bin(coord2, startIndex2, binSize, shift)
            # bin1 = int(splitLine[1]) // binSize + startIndex1
            # bin2 = int(splitLine[3]) // binSize + startIndex2

            # check if the entry is already in the map, if it is increase the count number
            if (bin1, bin2) in chrEntry:
                chrEntry[(bin1, bin2)] = chrEntry[(bin1, bin2)] + int(splitLine[4])
            elif (bin2, bin1) in chrEntry:
                chrEntry[(bin2, bin1)] = chrEntry[(bin2, bin1)] + int(splitLine[4])
            else:
                chrEntry[(bin1, bin2)] = int(splitLine[4])

        outputFile = os.path.join(outputDir, filename)
        out = open(outputFile, "w")
        for key, val in chrEntry.items():
            out.write("{} {} {}\n".format(key[0], key[1], val))
        out.close()


def main():
    parser = ArgumentParser("GenerateEdges",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument("--input", required=True, help="Path for the root directory of the data files")
    parser.add_argument("--output", required=True, help="output directory")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--bin-size", required=True)
    parser.add_argument("--shift", required=True)

    args = parser.parse_args()
    dataDir = args.input
    outputDir = args.output
    metadata = args.metadata
    shift = convert(args.shift)

    binSize = convert(args.bin_size)
    # if args.bin_size == "1M":
    #     binSize = 1000000
    # elif args.bin_size == "500k":
    #     binSize = 500000
    # elif args.bin_size == "100k":
    #     binSize = 100000
    # else:
    #     sys.exit("Bin size is not defined.")

    # chromosome bins metada
    chromBin_file = "{}/chrom_bins_{}_shift_{}.txt".format(metadata, args.bin_size, args.shift)

    # create output directory if not exists
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # store bin indexes
    chromBin = {}
    with open(chromBin_file) as f:
        for line in f:
            splitLine = line.split()
            chromBin[splitLine[0]] = splitLine[1]

    print("Processing Files...")
    executor = ThreadPoolExecutor(40)
    for filename in os.listdir(dataDir):
        executor.submit(processFile, filename, dataDir, chromBin, binSize, outputDir, shift)


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


def find_bin(coord, start_index, bin_size, shift):
    bin_num = coord // bin_size
    remainder = coord % bin_size

    bin_num += start_index

    if shift == 0:
        return bin_num

    if remainder >= shift:
        bin_num += 1

    return bin_num


if __name__ == '__main__':
    main()
