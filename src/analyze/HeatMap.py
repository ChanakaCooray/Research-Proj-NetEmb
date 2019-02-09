# import numpy as np
# import os
# from scipy import sparse
# import matplotlib.pyplot as plt
# import matplotlib
# import pandas as pd
# import seaborn as sns
# import matplotlib.cm as cm
#
# # matplotlib.use('TkAgg')
#
#
# def main():
#     # chromosome bins metada
#     # chromBin_file = "metadata/chrom_bins.txt"
#
#     # data file
#     # data_file = "analyze/1CDES_p3.H6"
#     # data_file = "analyze/1CDU.44"
#     # data_file = "analyze/1CDX1.71"
#     # data_file = "analyze/1CDX2.372"
#     data_file = "analyze/1CDX4.206"
#
#     df = pd.read_csv(data_file, sep=" ", header=None, names=["bin1", "bin2", "val"])
#
#     df = df.drop(df[df.bin1 == df.bin2].index)
#     plt.scatter(df["bin1"], df["bin2"], s=df["val"], c=df["val"], cmap="viridis")
#     plt.colorbar()
#
#     plt.show()
#
#
# if __name__ == '__main__':
#     main()

import numpy as np
import os
from scipy import sparse
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm

# matplotlib.use('TkAgg')


def main():
    # chromosome bins metada
    chromBin_file = "metadata/chrom_bins.txt"
    # chromosome sizes file

    # data file
    # data_file = "analyze/1CDES_p3.H6"
    # data_file = "analyze/1CDU.44"
    # data_file = "analyze/1CDX1.71"
    # data_file = "analyze/1CDX2.372"
    # data_file = "analyze/1CDX4.206"
    analyze_dir = "analyze"
    # outputFile = "output/analyze/components2.txt"

    # search for last index
    with open(chromBin_file) as f:
        for last in f: pass
        n = int(last.split()[1]) + 15

    for filename in os.listdir(analyze_dir):

        # data = {}
        # with open(os.path.join(analyze_dir,filename)) as f:
        #     for line in f:
        #         splitLine = line.split()
        #         edge1 = int(splitLine[0])
        #         edge2 = int(splitLine[1])
        #         if edge1 != edge2:
        #             data[(edge1, edge2)] = int(splitLine[2])

        df = pd.read_csv(os.path.join(analyze_dir,filename), sep=" ", header=None, names=["bin1", "bin2", "val"])

        # matrix = np.zeros((n, n), dtype=np.int)
        #
        # for key, val in data.items():
        #     matrix[key[0]][key[1]] = val
        #     matrix[key[1]][key[0]] = val
        #
        # np.fill_diagonal(matrix, 0)

        # plt.figure(figsize=(20, 20))
        # ax = sns.heatmap(matrix, robust=True)

        sns.jointplot(x=df["bin1"], y=df["bin2"], kind='kde')

        plt.savefig("output/analyze/kde-{}.png".format(filename))

        # plt.show()


if __name__ == '__main__':
    main()
