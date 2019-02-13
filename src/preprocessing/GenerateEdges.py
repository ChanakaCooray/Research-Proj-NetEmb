import os
from concurrent.futures import ThreadPoolExecutor
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def processFile(filename, dataDir, chromBin, binSize, outputDir):
    print("Processing {}...".format(filename))
    with open(os.path.join(dataDir, filename)) as f:
        next(f)
        chrEntry = {}
        for line in f:
            splitLine = line.split()

            startIndex1 = int(chromBin[splitLine[0]])
            startIndex2 = int(chromBin[splitLine[2]])

            # figure out the correct bin
            bin1 = int(splitLine[1]) // binSize + startIndex1
            bin2 = int(splitLine[3]) // binSize + startIndex2

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

    args = parser.parse_args()
    dataDir = args.input
    outputDir = args.output
    metadata = args.metadata
    binSize = 0
    if args.bin_size=="1M":
        binSize = 1000000
    elif args.bin_size=="500k":
        binSize = 500000
    elif args.bin_size=="100k":
        binSize = 100000

    # chromosome bins metada
    chromBin_file = "{}/chrom_bins.txt".format(metadata)

    # create output directory if not exists
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # store bin indexes
    chromBin = {}
    with open(chromBin_file) as f:
        for line in f:
            splitLine = line.split()
            chromBin[splitLine[0]] = splitLine[1]

    # print(chromBin["M"])

    print("Processing Files...")
    executor = ThreadPoolExecutor(40)
    for filename in os.listdir(dataDir):
        executor.submit(processFile, filename, dataDir, chromBin, binSize, outputDir)


if __name__ == '__main__':
    main()
