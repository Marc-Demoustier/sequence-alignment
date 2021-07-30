#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import random

def parse_arguments():
    """
    Parse all the arguments we need for this function

        Returns:
            args: a dict like with all the arguments parsed
    """
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
    """
    Read the input file and get the sequence from it

        Parameters:
            infile (string): path to the input file

        Returns:
            sequence (string): the sequence from the file
    """
    sequence = ''

    with open(infile, 'r') as input_file:
        line = input_file.readline()
        while (line):
            line = input_file.readline()[:-1] # remove the \n at the end of each line
            sequence += line
    return sequence


def burrow_wheeler_transform(sequence, frequency, k):
    """
    Apply the Burrow Wheeler Transform algorithm on a sequence and return it

        Parameters:
            sequence (string): the sequence to transform
            frequency (int): the frequency to use for the position indexes
            k (int): if it is not None then we apply the progressive version of the algorithm by cutting the suffix
                     array in k parts

        Returns:
            bwt (string): the Burrow Wheeler Transform of the given sequence
            position_index (int array): array containing all the position indexes (created with the frequency)
    """
    sequence += '$'

    # progressive case
    if k:
        bwt = ''
        frequency_counter = 0
        positions_index = []

        # create splitters array for a sample sort
        splitters = sorted(sequence[i:] + sequence[:i] for i in random.sample(range(len(sequence)), k))
        splitters.insert(0, chr(0))
        splitters.append(chr(127))

        for i in range(len(splitters) - 1):
            bucket = []
            for j in range(len(sequence)):
                rotation = sequence[j:] + sequence[:j]
                # add elements in the right bucket
                if splitters[i] <= rotation < splitters[i + 1]:
                    bucket.append(rotation)

            # sort the bucket and get the corresponding Burrow Wheeler Transform
            bucket = sorted(bucket)
            last_column = [row[-1:] for row in bucket]  # Last characters of each row

            # append the Burrow Wheeler Transform with all the buckets
            bwt += ''.join(last_column)

            # append the positions index with a step of frequency
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
    """
    Compress a Burrow Wheeler Transform sequence and returns it

        Parameters:
            bw (string): the Burrow Wheeler Transform sequence
            p (int): The number of A to add at the end to the sequence for having a length divisible by 4

        Returns:
            bw_compressed (array of bytes): the Burrow Wheeler Transform compressed sequence
    """
    compress_dict = {
        'A': 0b00,
        'C': 0b01,
        'G': 0b10,
        'T': 0b11
    }

    bw = bw.replace('$', '') # remove the $
    bw += 'A' * p

    bw_compressed = []
    # compress 4 character in 1 byte
    for i in range(0, len(bw), 4):
        byte_char = compress_dict[bw[i + 3]]
        for j in range(2, -1, -1):
            byte_char = byte_char << 2
            byte_char += compress_dict[bw[i + j]]
        bw_compressed.append(byte_char)

    return bw_compressed


def generate_outfile(args, sequence):
    """
    Apply the Burrow Wheeler Transform to the sequence and generate the output file

        Parameters:
            args (dict like): the parsed arguments
            sequence (string): the sequence to transform
    """
    bw, positions_index = burrow_wheeler_transform(sequence, args.f, args.progressive)

    c = 1 if args.compress else 0
    n = bw.index('$') if args.compress else 0
    p = 4 - (len(bw) - 1) % 4 if args.compress and (len(bw) - 1) % 4 != 0 else 0
    f = args.f

    positions_index_str = ','.join(positions_index)

    # write the output file
    with open(args.outfile, 'w') as output_file:
        output_file.write(f'{c} {n} {p} {f}\n')
        output_file.write(f'{positions_index_str}\n')

    if args.compress:
        # complete the output file with the compressed sequence
        bw_compressed = compress_sequence(bw, p)

        with open(args.outfile, 'ab') as output_file:
            output_file.write(bytearray(bw_compressed))
    else:
        with open(args.outfile, 'a') as output_file:
            output_file.write(bw)


def bw_build():
    """
    Parse the arguments, read the input file and generate the output file with the Burrow Wheeler Transform
    """
    args = parse_arguments()
    sequence = read_infile(args.infile)
    generate_outfile(args, sequence)


if __name__ == "__main__":
    bw_build()