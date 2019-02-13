import os
from concurrent.futures import ThreadPoolExecutor
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def processFile(subdir, dirName, fends, outputDir):
    out = open(os.path.join(outputDir, dirName), "w")
    out.write("chr1 coord1 chr2 coord2 count\n")
    with open(os.path.join(subdir, "adj")) as f:
        print("Processing {}...".format(os.path.join(subdir, "adj")))
        next(f)
        for line in f:
            splitLine = line.split()
            out.write(fends[splitLine[0]][0] + " " + fends[splitLine[0]][1] + " " + fends[splitLine[1]][0] + " " +
                      fends[splitLine[1]][1] + " 1\n")
    out.close()


def main():
    parser = ArgumentParser("ChromosomeIdentifier",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument("--input", required=True, help="Path for the root directory of the data files")
    parser.add_argument("--output", required=True, help="output directory")
    parser.add_argument("--metadata", required=True)

    args = parser.parse_args()
    dataDir = args.input
    outputDir = args.output
    metadata = args.metadata

    # Path for GATC.fends file
    fendsFile = "{}/GATC.fends".format(metadata)
    # path to metadata of Diploids 2i
    metaDip2iFile = "{}/Diploids_2i.txt".format(metadata)
    # path to metadata of Diploids serum
    metaDipSerumFile = "{}/Diploids_serum.txt".format(metadata)
    # path to metadata of Haploids
    metaHaploidsFile = "{}/Haploids.txt".format(metadata)

    # create output directories if not exist
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # store metadata of Diploids 2i
    metaDip2i = {}
    print("Processing Metadata {}".format(metaDip2iFile))
    with open(metaDip2iFile) as f:
        next(f)
        for line in f:
            splitLine = line.split()
            if splitLine[1] == "1":
                metaDip2i[splitLine[0]] = (splitLine[1], splitLine[2])

    # store metadata of Diploids Serum
    metaDipSerum = {}
    print("Processing Metadata {}".format(metaDipSerumFile))
    with open(metaDipSerumFile) as f:
        next(f)
        for line in f:
            splitLine = line.split()
            if splitLine[1] == "1":
                metaDipSerum[splitLine[0]] = (splitLine[1], splitLine[2])

    # store metadata of Haploids
    metaHap = {}
    print("Processing Metadata {}".format(metaHaploidsFile))
    with open(metaHaploidsFile) as f:
        next(f)
        for line in f:
            splitLine = line.split()
            if splitLine[2] == "1":
                metaHap[splitLine[1]] = (splitLine[0], splitLine[2], splitLine[3])

    # create a map using the GATC.fends file
    print("Processing GATC.fends file...")
    fends = {}
    with open(fendsFile) as f:
        next(f)
        for line in f:
            splitLine = line.split()
            fends[splitLine[0]] = (splitLine[1], splitLine[2])
        print(len(fends))

    print("Processing data...")
    executor = ThreadPoolExecutor(40)
    for subdir, dirs, files in (os.walk(dataDir)):
        dirName = subdir.split("/")[-1]
        dirNameTem = dirName.replace(".", "_")

        if dirNameTem not in metaDip2i and dirNameTem not in metaDipSerum and dirNameTem not in metaHap:
            continue

        executor.submit(processFile, subdir, dirName, fends, outputDir)


if __name__ == '__main__':
    main()
