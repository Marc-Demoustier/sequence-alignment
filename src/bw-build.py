import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description='Apply the Burrow-Wheeler transform on the infile and output the '
                                               'transformation in the outfile')
    parser.add_argument('infile', action='store', help='path to the fasta input file')
    parser.add_argument('outfile', action='store', help='path to the output file')
    parser.add_argument('f', action='store', help='define the index creation frequency', type=int)
    parser.add_argument('--compress', help='save the compressed version of the Burrow-Wheeler transform',
                        action='store_true')
    return parser.parse_args()


def read_infile(infile):
    sequence = ''

    with open(infile, 'r') as input_file:
        # process the file
        line = input_file.readline()
        while (line):
            line = input_file.readline()[:-1] # remove the \n at the end of each line
            sequence += line
    return sequence


def burrow_wheeler_transform(sequence):
    sequence += '$'

    table = sorted(sequence[i:] + sequence[:i] for i in range(len(sequence)))  # Table of rotations of string
    last_column = [row[-1:] for row in table]  # Last characters of each row
    return "".join(last_column)


def generate_outfile(args, bw):
    c = 1 if args.compress else 0
    n = bw.index('$') if args.compress else 0
    p = 0 # FIXME: A adds at the end of the transform with compression for being a multiple of 4
    f = args.f

    with open(args.outfile, 'w') as output_file:
        output_file.write(f'{c} {n} {p} {f}\n')
        output_file.write('0\n')
        output_file.write(burrow_wheeler_transform(bw))


def bw_build():
    args = parse_arguments()

    sequence = read_infile(args.infile)
    bw = burrow_wheeler_transform(sequence)
    generate_outfile(args, bw)


if __name__ == "__main__":
    bw_build()