# import pandas as pd
import math
import time
from scipy.stats import binom
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

# import matplotlib.pyplot as plt
# import numpy as np
# import matplotlib
# import seaborn as sns
#
# bin_size = "500k"
# cell_type = "1CDX1"
# shift = "100k"
#
# # matrix_file = "analyze-{}/SumMatrix/sum_matrix_{}_{}.txt".format(bin_size, bin_size, cell_type)
# # matrix_file = "analyze-{}/SumMatrix/test/sum_matrix_{}_{}test.txt".format(bin_size, bin_size, cell_type)
# matrix_file = "analyze-{}/SumMatrix/test/sum_matrix_{}_{}_shift_{}.txt".format(bin_size, bin_size, cell_type, shift)
# data = np.genfromtxt(matrix_file, delimiter=" ")
#
# np.fill_diagonal(data, 0)
#
# sum_by_row = data.sum(axis=0).astype(int)
# print(sum_by_row.shape)
# print(sum_by_row[6])
# total_sum = sum_by_row.sum(axis=0).astype(int)
#
# print(total_sum)

# matplotlib.use('TkAgg')
#
# a = np.random.random((16, 16))
# print(a)
# # ax = sns.heatmap(a, linewidth=0.5)
# plt.imshow(a)
# plt.show()

from concurrent.futures import ThreadPoolExecutor


def a(i):
    # print("1")
    # print("2")
    # print("3")
    l = [i + 1, i - 1]
    return l


def calculate_f_sum(n, t, p_max):
    sum_f = 0

    for i in range(t, n + 1):
        # print(i)
        sum_f += nCr(n, i) * (p_max ** i) * ((1 - p_max) ** (n - i))

    return sum_f


def calculate_pmax(M, sum_max):
    p_max = sum_max / M
    return p_max


def calculate_f_sum2(n, t, p_max):
    # sum_f = 0
    #
    # for i in range(t, n + 1):
    #     # print(i)
    #     p_value = 1-binom.cdf(i-1, n, p_max)
    #     sum_f += p_value
    #     # sum_f += nCr(n, i) * (p_max ** i) * ((1 - p_max) ** (n - i))

    sum_f = 1 - binom.cdf(t-1, n, p_max)

    return sum_f


def calculate_f(n, t, p_max):
    return nCr(n, t) * (p_max ** t) * ((1 - p_max) ** (n - t))


# helper function to calc nCr
def nCr(n, r):
    f = math.factorial
    return f(n) // f(r) // f(n - r)


def main():
    matrix = np.array([ [1,2,3,4],
                        [5,6,7,8],
                        [9,10,11,12],
                        [13,14,15,16]])

    rows = matrix.shape[0]
    cols = matrix.shape[1]

    for i in range(0, rows):
        for j in range(i+1, cols):
            print(matrix[i][j])


if __name__ == '__main__':
    main()
