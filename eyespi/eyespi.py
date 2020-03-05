#!/usr/bin/python3

import sys
import getopt
import subprocess
from Bio import SeqIO


def main(argv):
    """"This script removes a target sequence from an input fasta based on a blastn search"""
    inputfile = ''  # type: str
    matchfile = ''  # type: str
    identity = ''   # type: int
    match_length = ''  # type: int
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:c:x:m:o:", ["ifile=", "cfile=", "ival1=", "ival2=", "ofile="])
    except getopt.GetoptError:
        print ('presence-cut.py -i <inputfile> -c <matchfile> -x <identity> -m <match_length> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('presence.py -i <inputfile> -c <matchfile> -x <identity> -m <match_length>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-c", "--cfile"):
            matchfile = arg
        elif opt in ("-x", "--ival1"):
            identity = arg
        elif opt in ("-m", "--ival2"):
            match_length = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    print ('Input file is', inputfile)
    print ('Sequence to search for is', matchfile)
    print ('blast identity is', identity, '%')

    command = ['blastn', '-query', matchfile,'-subject', inputfile,'-perc_identity', identity,
               '-outfmt', '6', '-max_target_seqs', '1']

    output = subprocess.Popen(command, stdout=subprocess.PIPE)
    output = output.communicate()[0]
    output = output.decode()
    firstline = output.split("\n")[0]
    length_match = firstline.split("\t")[3]
    start_match = firstline.split("\t")[8]
    end_match = firstline.split("\t")[9]

    if length_match > match_length:
        print('SGI4 Present')
        print('start =', start_match)
        print('end =', end_match)
        print('length =', length_match)
        if end_match > start_match:
            record = SeqIO.read(inputfile, "fasta")
            with open(outputfile, "w") as out:
                SeqIO.write(record[:int(start_match) + 1] + record[int(end_match) + 1:], out, "fasta")
        elif start_match > end_match:
            record = SeqIO.read(inputfile, "fasta")
            with open(outputfile, "w") as out:
               SeqIO.write(record[:int(end_match) + 1] + record[int(start_match) + 1:], out, "fasta")
    elif length_match < match_length:
        print('SGI4 absent')

    print(firstline)
    print('length =', length_match)
    print('start =', start_match)
    print('end =', end_match)


if __name__ == "__main__":
    main(sys.argv[1:])
