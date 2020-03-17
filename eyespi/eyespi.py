#!/usr/bin/python3

import sys
import subprocess
from Bio import SeqIO
import argparse
import logging
import os 
import csv 
import shutil
import re
import textwrap

def run_blast(matchfile, inputfile, output_dir, identity='90' ):
    log = logging.getLogger('run_blast')

    if not os.path.exists(output_dir):
        log.debug(f'Creating output dir {output_dir}')
        os.mkdir(output_dir)
    if not os.path.exists(os.path.join(output_dir, 'scratch')):
        log.debug(f'Creating output SCRATCH dir {output_dir}')
        os.mkdir(os.path.join(output_dir, 'scratch'))

    out_file = os.path.join(output_dir, 'scratch', f'{os.path.basename(inputfile)}vs{os.path.basename(matchfile)}.tab')
    if os.path.exists(out_file):
        log.debug(f'Using cached blast results: {out_file}')
        return open(out_file).read()
    else:
        command = ['blastn', '-query', matchfile,'-subject', inputfile,'-perc_identity', identity,
                '-outfmt', '6']
        log.debug(f"Running {' '.join(command)}")
        output = subprocess.Popen(command, stdout=subprocess.PIPE)
        output = output.communicate()[0]
        output = output.decode()
        with open(out_file, 'w') as f:
            f.write(output)
        return output

def check_startup():
    return True

def eyespi(args):
    """"This script removes a target sequence from an input fasta based on a blastn search"""
    log = logging.getLogger('eyespi')

    inputfolder = args.ifile 
    matchfile = args.cfile
    identity = args.ival1
    output_dir = args.output_dir
    # Testing 
    inputfolder = 'spi_seq/'
    matchfile_dir = 'test-data/'
    ref_length_percent = 0.5 
    log.setLevel(logging.DEBUG)
    ## 
    # Create output dir. 
    # Construct list of all regions from reference 
    # Check reference files are ok (and format if needed)
    all_values = []
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    outfile = os.path.join(output_dir, args.output_prefix + '.table.csv')
    headers = [] 
    count = 0 
    easy_file = os.path.join(output_dir, args.output_prefix + '.easyfig.sh')
    easy_file_handle = open(easy_file, 'w')
    island_names =  [ ]
    for query in os.listdir(matchfile_dir):
        count += 1
#        if count > 30:
#            continue
        if query.endswith('.fasta'):#  or query.endswith('.fa'):
            values = dict(name=query.split('.')[0])
            group = re.match('ragout_(\d+)_all', values['name'])
            if group:
                values['name'] = f'mysnps_{group.group(1)}'

            matchfile = os.path.join(matchfile_dir, query)
            
            for seq in sorted(os.listdir(inputfolder)):
                inputfile = os.path.join(inputfolder,seq)
                clean_input = None
                if inputfile.endswith('.gbk'):
                    clean_input = inputfile + '.fna'
                    if not os.path.exists(os.path.exists(inputfile + '.fna')):
                        SeqIO.convert(inputfile, 'genbank', inputfile + '.fna', 'fasta')
                elif inputfile.endswith('.fasta'):
                    clean_input = inputfile
                if clean_input:
                    island_name = os.path.basename(clean_input).split('.')[0]
                    island_names.append(island_name)
                    score = [0.0]
                    ref_length = sum([len(x) for x in SeqIO.parse(clean_input, 'fasta')])
                    output = run_blast(matchfile, clean_input, output_dir, identity)
                    
                    blast_results  = output.split("\n")
                    easy_dir = os.path.join(output_dir, island_name)     
                    if not os.path.exists(easy_dir):
                        os.mkdir(easy_dir)
                    easy_dir_ref = os.path.join(output_dir, island_name, island_name + '.gbk')
                    if not os.path.exists(easy_dir_ref):
                        shutil.copy(inputfile, easy_dir_ref)
                    easy_seq_file = os.path.join(output_dir, island_name, values['name'] + '.fasta')
                    loci_coverage = [0] * ref_length
                    clips = [] 
                    
                    if len(blast_results) > 1:
                        log.debug(f'Multple BLAST Result: {len(blast_results)}')
                    log.debug(f'Running SPI: {island_name} {ref_length}')
                    for result in blast_results:
                        log.debug(f'BLAST Result: {result}')
                        if len(result.split('\t')) > 9:
                            match_name = result.split("\t")[0]
                            length_match = int(result.split("\t")[3])
                            start_match = int(result.split("\t")[8])
                            end_match = int(result.split("\t")[9])
                            if start_match > end_match:
                                t_end_match = end_match
                                end_match = start_match
                                start_match = t_end_match
                            query_start_match = result.split("\t")[6]
                            query_end_match = result.split("\t")[7]
                            if query_start_match > query_end_match:
                                t_query_end_match = query_end_match
                                query_end_match = query_start_match
                                query_start_match = t_query_end_match                                                    
                            if sum(loci_coverage[start_match : end_match]) < 10:

                                this_score = length_match / float(ref_length)
                                if this_score > 1.0: 
                                    this_score = 1.0
                                if this_score > ref_length_percent:
                                    records = [x for x in SeqIO.parse(matchfile, "fasta") if x.name == match_name]
                                    score.append(this_score)
                                    loci_coverage[start_match:end_match]  = [1] * (end_match - start_match)

                                    if records:
                                        record = records[0]
                                        clip = record[int(query_start_match): int(query_end_match)]
                                        clip.name += f':{query_start_match}-{query_end_match}'
                                            #if query_start_match > query_end_match:
                                            #    clip = record[int(query_end_match) + 1]:+ record[int(query_start_match) + 1:]
                                            #else:
                                            #    clip = record[:int(query_start_match) + 1] + record[int(query_end_match) + 1:]
                                        clips.append(clip)

                                else:
                                    score.append(0.0)
                        else: 
                            score.append(0.0)
                    values[island_name] = round((sum(score) / float(len(score))) * 100, 2)
                    if len(clips) > 1 : 
                        print('multiple alignment')
                    with open(easy_seq_file, "w") as out_seq:
                        final_string = ''
                        for clip in clips:
                            final_string += f'>{clip.name}\n{textwrap.fill(str(clip.seq), 90)}\n'
                        out_seq.write(final_string)
                        #SeqIO.write(clips, out, "fasta")
            ref_file_name = os.path.basename(easy_dir_ref)

            headers = list(set(headers + list(values.keys())))
            all_values.append(values)
    for island in list(set(island_names)):
        easy_file_handle.write(f"python2 ../Easyfig.py  -o {island}.bmp  -f1 T -legend top -leg_name gene  {island}/{island}.gbk {island}/A130.fasta {island}/S09304.fasta {island}/SO7676.fasta {island}/SO9207.fasta {island_name}/L01157.fasta {island}/O1960.fasta\n")            
    table = csv.DictWriter(open(outfile, 'w'),  fieldnames = sorted(headers))
    table.writeheader()
    table.writerows(all_values)


if __name__ == "__main__":
    logging.basicConfig()
    if check_startup():
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--ifile', help='input file')
        parser.add_argument('-c', '--cfile', help='match file')
        parser.add_argument('-x', '--ival1', help='identity', default = '90' )        
        parser.add_argument('-o', '--output_dir', help='output directory', default='eyespy_out')  
        parser.add_argument('-p', '--output_prefix', help='output prefix', default='myspi')  
        parser.add_argument('-v', '--verbose', help='verbose') 
        args = parser.parse_args()
        eyespi(args)
