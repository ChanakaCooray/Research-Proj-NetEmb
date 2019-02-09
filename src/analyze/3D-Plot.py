import os
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits import mplot3d

def color_1M(x):
    if x < 198:
        return "#e6194B"
    elif 198 <= x < 380:
        return "#3cb44b"
    elif 380 <= x < 540:
        return "#ffe119"
    elif 540 <= x < 696:
        return "#4363d8"
    elif 696 <= x < 849:
        return "#f58231"
    elif 849 <= x < 999:
        return "#911eb4"
    elif 999 <= x < 1152:
        return "#42d4f4"
    elif 1152 <= x < 1284:
        return "#f032e6"
    elif 1284 <= x < 1409:
        return "#bfef45"
    elif 1409 <= x < 1539:
        return "#fabebe"
    elif 1539 <= x < 1661:
        return "#469990"
    elif 1661 <= x < 1783:
        return "#e6beff"
    elif 1783 <= x < 1904:
        return "#9A6324"
    elif 1904 <= x < 2030:
        return "#fffac8"
    elif 2030 <= x < 2134:
        return "#800000"
    elif 2134 <= x < 2233:
        return "#aaffc3"
    elif 2233 <= x < 2329:
        return "#808000"
    elif 2329 <= x < 2420:
        return "#ffd8b1"
    elif 2420 <= x < 2482:
        return "#000075"
    elif 2482 <= x < 2649:
        return "#a9a9a9"
    elif x >= 2649:
        return "#000000"


def color_100k(x):
    if x < 1972:
        return "#e6194B"
    elif 1972 <= x < 3790:
        return "#3cb44b"
    elif 3790 <= x < 5386:
        return "#ffe119"
    elif 5386 <= x < 6943:
        return "#4363d8"
    elif 6943 <= x < 8469:
        return "#f58231"
    elif 8469 <= x < 9965:
        return "#911eb4"
    elif 9965 <= x < 11491:
        return "#42d4f4"
    elif 11491 <= x < 12809:
        return "#f032e6"
    elif 12809 <= x < 14050:
        return "#bfef45"
    elif 14050 <= x < 15350:
        return "#fabebe"
    elif 15350 <= x < 16569:
        return "#469990"
    elif 16569 <= x < 17782:
        return "#e6beff"
    elif 17782 <= x < 18985:
        return "#9A6324"
    elif 18985 <= x < 20237:
        return "#fffac8"
    elif 20237 <= x < 21272:
        return "#800000"
    elif 21272 <= x < 22256:
        return "#aaffc3"
    elif 22256 <= x < 23209:
        return "#808000"
    elif 23209 <= x < 24117:
        return "#ffd8b1"
    elif 24117 <= x < 24731:
        return "#000075"
    elif 24731 <= x < 26398:
        return "#a9a9a9"
    elif x >= 26398:
        return "#000000"


def main():
    analyze_dir = "NET-Output/100k/Node2vec-3"
    output_dir = "output/graphs/Node2vec3-100k"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(analyze_dir):
        df = pd.read_csv(os.path.join(analyze_dir, filename), delim_whitespace=True, header=None, skiprows=[0],
                         names=["node", "x", "y", "z"])

        # update here 1M 100k
        df2 = df["node"].apply(lambda x: color_100k(x))

        fig = plt.figure()
        ax = plt.axes(projection='3d')

        ax.scatter3D(df["x"], df["y"], df["z"], c=df2)

        # print(df2)
        # plt.figure(figsize=(20, 20))
        # plt.scatter(df["x"], df["y"], marker='o', c=df2)
        plt.savefig("{}/{}.png".format(output_dir, filename))
        # plt.show()

        # break


if __name__ == '__main__':
    main()
