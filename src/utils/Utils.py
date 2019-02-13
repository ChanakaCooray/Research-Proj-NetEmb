
def compare():
    with open('analyze-1M/1CDES_p3.H6', 'r') as file1:
        with open('temp/out/1CDES_p3.H6', 'r') as file2:
            same = set(file1).intersection(file2)

    same.discard('\n')

    with open('temp/out/out.txt', 'w') as file_out:
        for line in same:
            file_out.write(line)


if __name__ == '__main__':
    compare()
