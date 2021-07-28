import argparse
import random

def parse_arguments():
    parser = argparse.ArgumentParser(description='Apply the Burrow-Wheeler transform on the infile and output the '
                                               'transformation in the outfile')
    parser.add_argument('infile', action='store', help='path to the fasta input file')
    parser.add_argument('outfile', action='store', help='path to the output file')
    parser.add_argument('f', action='store', help='define the index creation frequency', type=int)
    parser.add_argument('--compress', help='save the compressed version of the Burrow-Wheeler transform',
                        action='store_true')
    parser.add_argument('--progressive', help='create the index with a progressive version using parameter k', type=int,
                        action='store', metavar='k')
    return parser.parse_args()


def read_infile(infile):
    sequence = ''

    with open(infile, 'r') as input_file:
        line = input_file.readline()
        while (line):
            line = input_file.readline()[:-1] # remove the \n at the end of each line
            sequence += line
    return sequence


def burrow_wheeler_transform(sequence, frequency, k):
    sequence += '$'

    if k:
        bwt = ''
        frequency_counter = 0
        positions_index = []

        # create splitters array
        splitters = sorted(sequence[i:] + sequence[:i] for i in random.sample(range(len(sequence)), k))
        splitters.insert(0, chr(0))
        splitters.append(chr(127))

        for i in range(len(splitters) - 1):
            bucket = []
            for j in range(len(sequence)):
                rotation = sequence[j:] + sequence[:j]
                if splitters[i] <= rotation < splitters[i + 1]:
                    bucket.append(rotation)

            bucket = sorted(bucket)
            last_column = [row[-1:] for row in bucket]  # Last characters of each row

            bwt += ''.join(last_column)

            for j in range(len(bucket)):
                if frequency_counter == 0:
                    positions_index.append(str(len(sequence) - bucket[j].index('$') - 1))

                frequency_counter = (frequency_counter + 1) % frequency

    else:
        table = sorted(sequence[i:] + sequence[:i] for i in range(len(sequence)))  # Table of rotations of string
        last_column = [row[-1:] for row in table]  # Last characters of each row
        bwt = ''.join(last_column)

        # create list of indexes (position of $)
        positions_index = [str(len(table) - table[i].index('$') - 1) for i in range(0, len(table), frequency)]

    return bwt, positions_index


def compress_sequence(bw, p):
    compress_dict = {
        'A': 0b00,
        'C': 0b01,
        'G': 0b10,
        'T': 0b11
    }

    bw = bw.replace('$', '') # remove the $
    bw += 'A' * p

    bw_compressed = []
    for i in range(0, len(bw), 4):
        byte_char = compress_dict[bw[i + 3]]
        for j in range(2, -1, -1):
            byte_char = byte_char << 2
            byte_char += compress_dict[bw[i + j]]
        bw_compressed.append(byte_char)

    return bw_compressed


def generate_outfile(args, sequence):
    bw, positions_index = burrow_wheeler_transform(sequence, args.f, args.progressive)

    c = 1 if args.compress else 0
    n = bw.index('$') if args.compress else 0
    p = 4 - (len(bw) - 1) % 4 if args.compress and (len(bw) - 1) % 4 != 0 else 0
    f = args.f

    positions_index_str = ','.join(positions_index)

    with open(args.outfile, 'w') as output_file:
        output_file.write(f'{c} {n} {p} {f}\n')
        output_file.write(f'{positions_index_str}\n')

    if args.compress:
        bw_compressed = compress_sequence(bw, p)

        with open(args.outfile, 'ab') as output_file:
            output_file.write(bytearray(bw_compressed))
    else:
        with open(args.outfile, 'a') as output_file:
            output_file.write(bw)


def bw_build():
    args = parse_arguments()
    sequence = read_infile(args.infile)
    generate_outfile(args, sequence)


if __name__ == "__main__":
    bw_build()