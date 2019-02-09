import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import os

matplotlib.use('TkAgg')


def main():
    analyze_dir = "analyze"
    output_file = "output/analyze/histogram.png"

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    fig, axes = plt.subplots(nrows=5, ncols=2)
    fig.tight_layout()

    counter = 0

    for filename in os.listdir(analyze_dir):
        file = os.path.join(analyze_dir, filename)
        df = pd.read_csv(file, sep=" ", header=None, names=["bin1", "bin2", "val"])

        ax1 = axes[counter][0]
        ax2 = axes[counter][1]

        ax1.hist(df["val"], bins=df['val'].max() - 1)
        ax1.set_xlabel('Values')
        ax1.set_title("{} with Diagonals".format(filename))

        df2 = df.drop(df[df.bin1 == df.bin2].index)

        ax2.hist(df2["val"], bins=df2['val'].max() - 1)
        ax2.set_xlabel('Values')
        ax2.set_title("{} without Diagonals".format(filename))

        counter += 1

    fig.set_size_inches(20, 15)
    plt.savefig(output_file, dpi=100)
    # plt.show()


if __name__ == '__main__':
    main()
