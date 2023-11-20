# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:00:37 2023

@author: Genglin Guo
@e-mail: 2019207025.njau.edu.cn
"""

import argparse
import pathlib
import time
import subprocess
import multiprocessing
from Bio import SeqIO

def get_argument():
    # Parsers
    parser = argparse.ArgumentParser(description = 'CPtoxin', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_group_1 = parser.add_argument_group('Input and Output')
    parser_group_2 = parser.add_argument_group('Parameters')

    # Input and output
    parser_group_1.add_argument('-i', '--input', required = True, nargs = '+', type = str, 
                                help = 'Input FASTA file')
    parser_group_1.add_argument('-r', '--reference', required = False, type = str, default = 'CPtoxin',
                                help = 'Reference CPtoxin sequence file')
    parser_group_1.add_argument('-o', '--output', required = False, type = str, default = 'CPtoxin_output.txt',
                              help = 'Output file')
    # Parameters
    parser_group_2.add_argument('-t', '--threads', required = False, type = int, default = min(multiprocessing.cpu_count(), 4), 
                        help = 'Threads to use for BLAST searches')
    parser_group_2.add_argument('--min_cov', required = False, type = float, default = 80.0, 
                               help = 'Minimum percentage coverage to consider a single gene complete. [default: 80.0%]')
    parser_group_2.add_argument('--min_id', required = False, type = float, default = 70.0, 
                               help = 'Minimum percentage identity to consider a single gene complete. [default: 70.0%]')
    return parser
    
def parse_refs(reference):
    ref_path = pathlib.Path(reference).resolve()
    refs = {}
    for contig in SeqIO.parse(ref_path, 'fasta'):
        refs[contig.name] = contig.seq
    return refs

def generate_temp_file(ref_name, ref_seq):
    with open('temp_blast_file.txt', 'wt') as file:
        file.write('>')
        file.write(ref_name)
        file.write('\n')
        file.write(str(ref_seq))
    return pathlib.Path('temp_blast_file.txt').resolve()
        
def run_blast(inputfile, refseq):
    input_path = pathlib.Path(inputfile).resolve()
    # Do blast, iterator the result to a list
    blast_hits = []
    command = ['blastn', '-query', refseq, '-subject', input_path, '-outfmt', 
               '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq']
    process = subprocess.run(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = process.stdout.decode()
    pathlib.Path(refseq).unlink()
    for line in line_iterator(out):
        blast_hits.append(BlastResult(line))
    return blast_hits

def line_iterator(line_breaks):
    # Handle the BLAST output and remove the line breaks 
    line = -1
    while True:
        nextline = line_breaks.find('\n', line + 1)
        if nextline < 0:
            break
        yield line_breaks[line + 1:nextline]
        line = nextline

class BlastResult(object):
    # Handle the BLAST output
    def __init__(self, hit_string):
        parts = hit_string.split('\t')
        self.length = int(parts[8])
        self.pident = float(parts[9])

def pending_result(blast_hits, min_cov, min_id, refseq):
    # Find the best serotype of inputfile by sequence alignment
    best_coverage = 0.0
    best_identity = 0.0
    for i in blast_hits:
        coverage = 100.0 * i.length / len(str(refseq))
        identity = i.pident
        if coverage > best_coverage:
            best_coverage = coverage
            best_identity = identity
        elif coverage == best_coverage and identity > best_identity:
            best_identity = identity
    if best_coverage >= min_cov and best_identity >= min_id:
        return '1'
    else:
        return '0'
    
def typing(toxincode):
    types = {
        '1000000' : 'A',
        '1010001' : 'B',
        '1000001' : 'C',
        '1100001' : 'C',
        '1010000' : 'D',
        '1110000' : 'D',
        '1001100' : 'E',
        '1101100' : 'E',
        '1001000' : 'E',
        '1000100' : 'E',
        '1101000' : 'E',
        '1100100' : 'E',
        '1100000' : 'F',
        '1000010' : 'G'
        }
    if toxincode in types:
        toxintype = types[toxincode]
    else:
        toxintype = 'NA'
    return toxintype

def generate_output(output):
    # Generate a blank output table file
    if pathlib.Path(output).is_file():
        return
    headers = ['Isolate', 'Toxintype', 'plc', 'cpe', 'etx', 'iap', 'ibp', 'netB', 'cpb']
    with open(output, 'wt') as file:
        file.write('\t'.join(headers))
        file.write('\n')

def output(output, inputfile, toxintype, toxincode):
    # Generate output
    inputfile = strip_suffix(inputfile)
    simple_output = inputfile + ' : '  + ' ' + 'toxintype' + ' ' + toxintype
    line = [inputfile, toxintype]
    for i in toxincode:
        if i == '1':
            line.append('+')
        else:
            line.append('-')
    print(simple_output)
    with open(output, 'at') as file:
        file.write('\t'.join(line))
        file.write('\n')

def strip_suffix(inputfile):
    # Strip the suffix of inputfile
    if inputfile.lower().endswith('.fa'):
        inputfile = inputfile[:-3]
    elif inputfile.lower().endswith('.fna'):
        inputfile = inputfile[:-4]
    elif inputfile.lower().endswith('.fas'):
        inputfile = inputfile[:-4]
    elif inputfile.lower().endswith('.fasta'):
        inputfile = inputfile[:-6]
    return inputfile

def main():
    starttime = time.perf_counter()
    # Initialize
    args = get_argument().parse_args()
    reference_dict = parse_refs(args.reference)
    for inputfile in args.input:
        toxincode = ''
        for ref_name, ref_seq in reference_dict.items():
            temp_file = generate_temp_file(ref_name, ref_seq)
            blast_hits = run_blast(inputfile, temp_file)
            presence = pending_result(blast_hits, args.min_cov, args.min_id, ref_seq)
            toxincode += presence
        toxintype = typing(toxincode)
    # Generate output
        generate_output(args.output)
        output(args.output, inputfile, toxintype, toxincode)
    endtime = time.perf_counter() - starttime
    per_genome_time = endtime / len(args.input)
    print('{:.1f}h{:.1f}m{:.1f}s for one subject'.format(per_genome_time // 3600, per_genome_time % 3600 // 60, per_genome_time % 60))
    print('Total time consumed : {:.1f}h{:.1f}m{:.1f}s'.format(endtime // 3600, endtime % 3600 // 60, endtime % 60))
   
main()