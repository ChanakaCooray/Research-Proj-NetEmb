import numpy as np
import pandas as pd


def compare_range(file1, file2):
    df_file1 = pd.read_csv(file1, sep=" ", header=None,
                           names=["chrm1", "chrm1_start_index", "chrm1_end_index", "chrm2", "chrm2_start_index",
                                  "chrm2_end_index"])

    df_file2 = pd.read_csv(file2, sep=" ", header=None,
                           names=["chrm1", "chrm1_start_index", "chrm1_end_index", "chrm2", "chrm2_start_index",
                                  "chrm2_end_index"])

    # print(df_file1.dtypes)
    # print(df_file2.dtypes)

    total = 0
    match = 0
    for index1, row1 in df_file1.iterrows():
        # print(row["chrm1"], row["chrm1_start_index"], row["chrm1_end_index"], row["chrm2"], row["chrm2_start_index"],
        #       row["chrm2_end_index"])
        # print(row1["chrm1_start_index"]-row1["chrm1_end_index"])
        total += 1
        match_found = False
        for index2, row2 in df_file2.iterrows():
            if row1["chrm1"] == row2["chrm1"]:
                if row1["chrm2"] == row2["chrm2"]:
                    if row1["chrm1_start_index"] <= row2["chrm1_end_index"] and row2["chrm1_start_index"] <= row1[
                        "chrm1_end_index"] and row1["chrm2_start_index"] <= row2["chrm2_end_index"] and row2[
                        "chrm2_start_index"] <= row1["chrm2_end_index"]:
                        match += 1
                        match_found = True
                        break

            if row1["chrm1"] == row2["chrm2"]:
                if row1["chrm2"] == row2["chrm1"]:
                    if row1["chrm1_start_index"] <= row2["chrm2_end_index"] and row2["chrm2_start_index"] <= row1[
                        "chrm1_end_index"] and row1["chrm2_start_index"] <= row2["chrm1_end_index"] and row2[
                        "chrm1_start_index"] <= row1["chrm2_end_index"]:
                        match += 1
                        match_found = True
                        break

        # if not match_found:
        #     print(row1["chrm1"], row1["chrm1_start_index"], row1["chrm1_end_index"], row1["chrm2"],
        #           row1["chrm2_start_index"], row1["chrm2_end_index"])

    print("{} {}".format(match, total))


if __name__ == '__main__':
    file1 = "output/range_output/output_1CDX2_500k.txt"
    file2 = "output/range_output/output_1CDX2_500k_shift_100k.txt"

    for i in range(1, 5):
        for j in range(1, 5):
            file1 = "output/range_output/output_1CDX{}_500k.txt".format(i)
            file2 = "output/range_output/output_1CDX{}_500k_shift_{}k.txt".format(i, j * 100)
            print("1CDX{} shift {}k".format(i, j * 100))
            compare_range(file1, file2)
