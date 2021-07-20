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
            line = input_file.readline()[:-1]
            sequence += line
    return sequence


def generate_outfile(args, sequence):
    c = 1 if args.compress else 0
    n = 0 # FIXME: $ position
    p = 0 # FIXME: A adds at the end of the transform with compression for being a multiple of 4
    f = args.f

    print(c, n, p, f)
    print(sequence)

def bw_build():
    args = parse_arguments()

    sequence = read_infile(args.infile)
    generate_outfile(args, sequence)

bw_build()