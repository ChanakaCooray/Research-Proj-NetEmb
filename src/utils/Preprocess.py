import sys
import os
import csv


def preprocess():
    data_dir = "new-data/data/pronucleus_male"
    output_dir = "new-data/preprocessed-data/pronucleus_male"

    # create output directory if not exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(data_dir):
        if not filename.endswith(".csv"):
            continue
        with open(os.path.join(data_dir, filename)) as f:
            csv_reader = csv.reader(f, delimiter=',')
            next(f)
            outputFile = os.path.join(output_dir, filename)
            out = open(outputFile, "w")
            out.write("chr1 coord1 chr2 coord2 count\n")
            for row in csv_reader:
                out.write("{} {} {} {} 1\n".format(row[0], row[2], row[1], row[3]))
            out.close()


if __name__ == '__main__':
    preprocess()
