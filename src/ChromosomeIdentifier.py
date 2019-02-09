import os
import pandas as pd

# Path for GATC.fends file
fendsFile = "metadata/GATC.fends"
# Path for the root directory of the data files
dataDir = "data.dir"
# path to metadata of Diploids 2i
metaDip2iFile = "metadata/Diploids_2i.txt"
# path to metadata of Diploids serum
metaDipSerumFile = "metadata/Diploids_serum.txt"
# path to metadata of Haploids
metaHaploidsFile = "metadata/Haploids.txt"

# output directory
outputDir = os.path.join("output", "chromosomeMap")

# create output directories if not exist
if not os.path.exists(outputDir):
    os.makedirs(outputDir)


# helper function to generate expected output
def generate_data(row, fends):
    fend1 = fends.loc[row['fend1']]
    fend2 = fends.loc[row['fend2']]

    return pd.Series({'chr1': fend1['chr'], 'coord1': fend1['coord'], 'chr2': fend2['chr'], 'coord2': fend2['coord'],
                      'count': 1})


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
df_fends = pd.read_csv(fendsFile, sep="\t", header=0, index_col=0, dtype={'fend': 'int', 'chr': 'str', 'coord': 'int'})

for subdir, dirs, files in os.walk(dataDir):
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

    # for fileName in files:
    #     print fileName
    df_adj = pd.read_csv(os.path.join(subdir, "adj"), sep="\t", header=0,
                         dtype={'fend1': 'int', 'fend2': 'int', 'count': 'int'})
    df_output = df_adj.apply(lambda row: generate_data(row, df_fends), axis=1)
    df_output.to_csv(os.path.join(outputDir, dirName), header=True, index=None, sep=' ',
                     columns=['chr1', 'coord1', 'chr2', 'coord2', 'count'])
