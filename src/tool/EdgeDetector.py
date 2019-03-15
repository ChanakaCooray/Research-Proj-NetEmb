import os
from concurrent.futures import ThreadPoolExecutor
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys


# process the data file to an edge list file of bins
def process_chrom_file(filename, data_dir, chrm_num, start_index, end_index, output_dir):
    print("Processing {}...".format(filename))
    with open(os.path.join(data_dir, filename)) as f:
        next(f)
        entries = []
        for line in f:
            split_line = line.split()

            chrm1 = split_line[0]
            chrm2 = split_line[2]

            coord1 = int(split_line[1])
            coord2 = int(split_line[3])

            if chrm1 == chrm_num:
                if start_index <= coord1 <= end_index:
                    entries.append(line)
            elif chrm2 == chrm_num:
                if start_index <= coord2 <= end_index:
                    entries.append(line)

        if len(entries) > 0:
            output_file = os.path.join(output_dir, filename)
            out = open(output_file, "w")
            for val in entries:
                out.write(val)
            out.close()


# process all the data files concurrently
def generate_data_files(data_dir, output_dir, chrm_num, start_index, end_index):
    print("Processing Files...")
    executor = ThreadPoolExecutor(40)
    for filename in os.listdir(data_dir):
        executor.submit(process_chrom_file, filename, data_dir, chrm_num, start_index, end_index, output_dir)
    executor.shutdown(wait=True)


# convert the values like 1M, 500k to real integer values
def convert(val):
    if val.isdigit():
        return int(val)

    lookup = {'k': 1000, 'M': 1000000, 'B': 1000000000}
    unit = val[-1]
    try:
        number = int(val[:-1])
    except ValueError:
        sys.exit("Value Error.")
    if unit in lookup:
        return lookup[unit] * number
    return int(val)


def parse_file(file, data_dir, output_dir):
    with open(file) as f:
        for line in f:
            split_line = line.split()

            chrm_num = split_line[1]
            start_index = split_line[2]
            end_index = split_line[3]

            output_dir_line = os.path.join(output_dir, "chrm{}_{}_{}".format(chrm_num, start_index, end_index))

            # create output directory if not exists
            if not os.path.exists(output_dir_line):
                os.makedirs(output_dir_line)

            generate_data_files(data_dir, output_dir_line, chrm_num, start_index, end_index)


def main():
    parser = ArgumentParser("EdgeDetector",
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            conflict_handler='resolve')
    parser.add_argument("--file", default='None', help="A file containing all the ranges to check")
    parser.add_argument("--data", required=True, help="Path for the root directory of the data files")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--chrm", required=False, help="Chromosome Number")
    parser.add_argument("--start-index", required=False, help="Start index of the range, Ex: 500k, 1.5M")
    parser.add_argument("--end-index", required=False, help="End index of the range, Ex: 500k, 1.5M")

    args = parser.parse_args()
    data_dir = args.data
    output_dir = args.output
    file = args.file

    if file != 'None':
        parse_file(file, data_dir, output_dir)
        return

    chrm_num = args.chrm
    start_index = convert(args.start_index)
    end_index = convert(args.end_index)

    # create output directory if not exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    generate_data_files(data_dir, output_dir, chrm_num, start_index, end_index)


if __name__ == '__main__':
    main()
