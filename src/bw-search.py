import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='search subsequence in a genome following the burrows wheeler '
                                                    'transform')
    parser.add_argument('infile', action='store', help='path to the index file generated by bw-build')
    parser.add_argument('q', action='store', help='sequence to search')
    parser.add_argument('--count-only', help='show only the correspondance number', action='store_true')
    return parser.parse_args()


def read_infile(infile):
    with open(infile, 'r') as input_file:
        c, n, p, f = input_file.readline().split()
        input_file.readline()
        bw = input_file.readline()
    return bw


def inverse_burrow_wheeler_transform(bw):
    table = [''] * len(bw)
    for _ in range(len(bw)):
        table = sorted(bw[i] + table[i] for i in range(len(bw)))
    sequence = [row for row in table if row.endswith('$')][0]
    return sequence.rstrip('$')


def bw_search():
    args = parse_arguments()
    bw = read_infile(args.infile)
    print(inverse_burrow_wheeler_transform(bw))


if __name__ == "__main__":
    bw_search()