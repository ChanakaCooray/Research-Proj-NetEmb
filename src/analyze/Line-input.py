import os


def main():
    analyze_dir = "analyze-1M"
    output_dir = "output/LINE-1M"

    # create output directory if not exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(analyze_dir):
        outputFile = os.path.join(output_dir, "LINE-{}".format(filename))
        out = open(outputFile, "w")
        with open(os.path.join(analyze_dir, filename)) as f:
            for line in f:
                splitLine = line.split()
                edge1 = int(splitLine[0])
                edge2 = int(splitLine[1])
                if edge1 != edge2:
                    out.write("{} {} {}\n".format(edge1, edge2, splitLine[2]))
                    out.write("{} {} {}\n".format(edge2, edge1, splitLine[2]))
        out.close()


if __name__ == '__main__':
    main()
