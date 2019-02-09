import os
from concurrent.futures import ThreadPoolExecutor

# input data file
dataDir = os.path.join("output", "chromosomeMap")
# chromosome bins metada
chromBin_file = "metadata/chrom_bins.txt"
# bin size
binSize = 100000
# output directory
outputDir = os.path.join("output", "EdgesBin")

# create output directory if not exists
if not os.path.exists(outputDir):
    os.makedirs(outputDir)

# store bin indexes
chromBin = {}
with open(chromBin_file) as f:
    for line in f:
        splitLine = line.split()
        chromBin[splitLine[0]] = splitLine[1]


def processFile(filename):
    print("Processing {}...".format(filename))
    with open(os.path.join(dataDir, filename)) as f:
        next(f)
        chrEntry = {}
        for line in f:
            splitLine = line.split()

            startIndex = int(chromBin[splitLine[0]])

            # figure out the correct bin
            bin1 = int(splitLine[1]) // binSize + startIndex
            bin2 = int(splitLine[3]) // binSize + startIndex

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


print("Processing Files...")
executor = ThreadPoolExecutor(40)
for filename in os.listdir(dataDir):
    executor.submit(processFile, filename)
