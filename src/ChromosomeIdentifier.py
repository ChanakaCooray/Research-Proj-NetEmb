import os
import pandas as pd
# import threading
# from concurrent.futures import ThreadPoolExecutor
# from time import sleep

# Path for GATC.fends file
fendsFile = "metadata/GATC.fends"
# Path for the root directory of the data files
dataDir = "../data"
# path to metadata of Diploids 2i
metaDip2iFile = "metadata/Diploids_2i.txt"
# path to metadata of Diploids serum
metaDipSerumFile = "metadata/Diploids_serum.txt"
# path to metadata of Haploids
metaHaploidsFile = "metadata/Haploids.txt"

# output directory
outputDir = os.path.join("../output", "chromosomeMap")

# create output directories if not exist
if not os.path.exists(outputDir):
    os.makedirs(outputDir)


# store metadata of Diploids abd Haploids
print("Processing Metadata {}".format(metaDip2iFile))
df_Dip2i = pd.read_csv(metaDip2iFile, sep="\t", header=0, index_col=0)
df_Dip2i.index = df_Dip2i.index.map(lambda x: x.replace("_", "."))

print("Processing Metadata {}".format(metaDipSerumFile))
df_DipSerum = pd.read_csv(metaDipSerumFile, sep="\t", header=0, index_col=0)
df_DipSerum.index = df_DipSerum.index.map(lambda x: x.replace("_", "."))

print("Processing Metadata {}".format(metaHaploidsFile))
df_Hap = pd.read_csv(metaHaploidsFile, sep="\t", header=0, index_col=1)
df_Hap.index = df_Hap.index.map(lambda x: x.replace("_", "."))

# create a map using the GATC.fends file
print("Processing GATC.fends file...")
fends = {}
with open(fendsFile) as f:
    next(f)
    for line in f:
        splitLine = line.split()
        fends[splitLine[0]] = (splitLine[1], splitLine[2])

print("Processing data...")
for subdir, dirs, files in (os.walk(dataDir)):
    dirName = subdir.split("/")[-1]

    if dirName not in df_Dip2i.index and dirName not in df_DipSerum.index and dirName not in df_Hap.index:
        continue

    if dirName in df_Dip2i.index:
        if df_Dip2i.loc[dirName]['passed_qc'] == 0:
            continue

    if dirName in df_DipSerum.index:
        if df_DipSerum.loc[dirName]['passed_qc'] == 0:
            continue

    if dirName in df_Hap.index:
        if df_Hap.loc[dirName]['passed_qc'] == 0:
            continue

    out = open(os.path.join(outputDir, dirName), "w")
    out.write("chr1 coord1 chr2 coord2 count\n")
    with open(os.path.join(subdir, "adj")) as f:
        "Processing " + os.path.join(subdir, "adj") + " ..."
        next(f)
        for line in f:
            splitLine = line.split()
            out.write(fends[splitLine[0]][0] + " " + fends[splitLine[0]][1] + " " + fends[splitLine[1]][0] + " " +
                      fends[splitLine[1]][1] + " 1\n")
