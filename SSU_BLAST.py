import argparse
import sys

#function to parse command line arguments #
def check_arg(args=None): 
	parser = argparse.ArgumentParser(description= 'Produces a brief analysis of RNA-seq data')
	parser.add_argument('-i', '--input',
		help='path to input file',
		required = 'True')
	parser.add_argument('-o', '--output',
		help = 'output file name',
		required = 'True')
	return parser.parse_args(args)

#retrieve command line arguments and assign to variables
args = check_arg(sys.argv[1:])
infile = args.input
outfile = args.output

con = open(infile)
use = con.read().split()

ids = list()
for i in range(2, len(use), 2):
	ids.append(use[i])

fin = list()
for i in range(0, len(ids)):
	fin.append(ids[i].replace('"', ''))

import Bio

from Bio import Entrez
from Bio import SeqIO

Entrez.email = 'cdotson@luc.edu'

hold1 = Entrez.efetch(db = 'nucleotide', id = fin, rettype = 'gb', retmode = 'text')
records = SeqIO.parse(hold1, "genbank") 
all_na = list()
for taxa in records:
	all_na.append(taxa.annotations['taxonomy'])

import csv

with open(outfile, 'w') as f:
	write = csv.writer(f)
	write.writerows(all_na)