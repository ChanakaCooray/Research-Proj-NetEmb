import pandas as pd

binSize = 100000
input = 'metadata/chrom_sizes.txt'
output = 'output/chrom_bins.txt'

df = pd.read_csv(input, sep="\t", header=None, names=["chrm", "size"])
df['size'] = df['size'].apply(lambda x: int(x) // binSize)

count = 0

for i in df.index:
    size = df.at[i, 'size']
    df.at[i, 'size'] = count
    count = count + size + 1

df.to_csv(output, header=None, index=None, sep=' ', mode='w')
