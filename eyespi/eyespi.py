#!/usr/bin/python3

import sys
import getopt
import subprocess
from Bio import SeqIO
import argparse
import logging

log = logging.getLogger()

def run_blast(matchfile, inputfile, identity='90' ):
    command = ['blastn', '-query', matchfile,'-subject', inputfile,'-perc_identity', identity,
               '-outfmt', '6', '-max_target_seqs', '1']
    log.debug(f"Running {' '.join(command)}")
    output = subprocess.Popen(command, stdout=subprocess.PIPE)
    output = output.communicate()[0]
    output = output.decode()
    return output

def check_startup():
    return True

def main(args):
    """"This script removes a target sequence from an input fasta based on a blastn search"""
    inputfile = args.ifile 
    matchfile = args.cfile
    identity = args.ival1
    match_length = args.ival2
    outputfile = args.output

    # Testing 
    inputfile = 'spi_seq/SPI-1.gbk'
    matchfile = 'test-data/A130.fasta'
    identity = '90'
    # This should be a proportion of the ref sequence 
    match_length = 90
    outputfile = 'test.fna'
    log.setLevel(logging.DEBUG)
    ## 
    if inputfile.endswith('.gbk'):
        clean_input = inputfile + '.fna'
        SeqIO.convert(inputfile, 'genbank', inputfile + '.fna', 'fasta')
    else:
        clean_input = inputfile
    output = run_blast(matchfile, clean_input, identity)
    firstline = output.split("\n")[0]
    length_match = firstline.split("\t")[3]
    start_match = firstline.split("\t")[8]
    end_match = firstline.split("\t")[9]

    if int(length_match) > int(match_length):
        print('start =', start_match)
        print('end =', end_match)
        print('length =', length_match)
        if end_match > start_match:
            record = SeqIO.read(clean_input, "fasta")
            with open(outputfile, "w") as out:
                SeqIO.write(record[:int(start_match) + 1] + record[int(end_match) + 1:], out, "fasta")
        elif start_match > end_match:
            record = SeqIO.read(clean_input, "fasta")
            with open(outputfile, "w") as out:
               SeqIO.write(record[:int(end_match) + 1] + record[int(start_match) + 1:], out, "fasta")
    elif length_match < match_length:
        print('SGI4 absent')

    print(firstline)
    print('length =', length_match)
    print('start =', start_match)
    print('end =', end_match)


if __name__ == "__main__":
    if check_startup():
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--ifile', help='input file')
        parser.add_argument('-c', '--cfile', help='match file')
        parser.add_argument('-x', '--ival1', help='identity')        
        parser.add_argument('-m', '--ival2', help='match length')  
        parser.add_argument('-o', '--output', help='output')  
        parser.add_argument('-v', '--verbose', help='verbose') 
        args = parser.parse_args()
        main(args)
