# import pandas as pd
#
# metaDip2iFile = "metadata/Diploids_2i.txt"
#
# # store metadata of Diploids 2i
# metaDip2i = {}
# print("Processing Metadata {}".format(metaDip2iFile))
# with open(metaDip2iFile) as f:
#     next(f)
#     for line in f:
#         splitLine = line.split()
#         if splitLine[1] == "1":
#             metaDip2i[splitLine[0].replace("_", ".")] = (splitLine[1], splitLine[2])
#
#
# print (len(metaDip2i))
#
# # store metadata of Diploids abd Haploids
# print("Processing Metadata {}".format(metaDip2iFile))
# df_Dip2i = pd.read_csv(metaDip2iFile, sep="\t", header=0, index_col=0)
#
# print(len(df_Dip2i))


# class Point:
#     def __init__(self, x, y):
#         self._x = x
#         self._y = y
#
#     def func(self, x):
#         print(x + str(self._x))
#
# def main2():
#     print ("aaa")
#
# # def main():
# #     print ("Aaa")
#
# # p = Point(1, 2)
# # p.func("aaa")
#
# # point = Point(444, 2)
#
# # print(type(point._x))
#
# if __name__ == "__main__":
#     main2()

import numpy as np
from scipy.io import savemat
from scipy import sparse
import pandas as pd

# data = [
#     (1, 2, 3),
#     (2, 3, 2),
#     (2, 2, 4)
# ]
#
# matrix = np.zeros((4, 4), dtype=np.int)
# for user, item, rating in data:
#     matrix[user][item] = rating
#     matrix[item][user] = rating
#
# np.fill_diagonal(matrix, 0)
#
# Asp = sparse.csr_matrix(matrix)
#
# savemat('output/temp', {'Asp':Asp})

# np.savetxt('output/outfile.txt', matrix, fmt="%d")

# with open('analyze/1CDES_p3.H6') as f:
#     next(f)
#     for line in f:
#         splitLine = line.split()
#
#         if int(splitLine[0]) >26398 or int(splitLine[1]) >26398:
#             print(line)

# with open('metadata/GATC.fends') as f:
#     next(f)
#     for line in f:
#         splitLine = line.split()
#
#         if int(splitLine[2]) == 183083638:
#             print(line)

# df = pd.read_csv("output/EdgesBin/1CDU.1", sep=" ", header=None, names=["bin1", "bin2", "val"])
# df2 = df.drop(df[df.bin1 == df.bin2].index)
#
# print(df2)

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns

matplotlib.use('TkAgg')

a = np.random.random((16, 16))
print(a)
# ax = sns.heatmap(a, linewidth=0.5)
plt.imshow(a)
plt.show()