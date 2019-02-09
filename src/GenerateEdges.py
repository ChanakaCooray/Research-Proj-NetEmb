import os
import pandas as pd

# input data file
dataDir = os.path.join("output", "chromosomeMap")
# chromosome bins metada
chrom_bin="metadata/chrom_bins.txt"
# bin size
binSize = 100000
# output directory
outputDir = os.path.join("output", "EdgesBin")

# create output directory if not exists
if not os.path.exists(outputDir):
    os.makedirs(outputDir)

df_chrombin = pd.read_csv(chrom_bin, sep=" ", header=None, names=["chrm", "startIndex"], index_col=0)

for filename in os.listdir(dataDir):
    with open(os.path.join(dataDir, filename)) as f:
        next(f)
        chrEntry = {}
        for line in f:
            splitLine = line.split()

            startIndex = df_chrombin.loc[splitLine[0]]['startIndex']

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

        for key, val in chrEntry.items():
            outputFile = os.path.join(outputDir, filename)
            out = open(outputFile, "a")
            out.write("{} {} {}\n".format(key[0], key[1], val))
            out.close()
