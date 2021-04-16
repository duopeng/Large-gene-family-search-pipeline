#!/usr/bin/python3
from collections import defaultdict
import argparse
import sys
import linecache
import re
from Bio import SeqIO
import collections

parser = argparse.ArgumentParser(description='This script prints out the length of all seqs in a fasta file.  please specify all optional parameters. ')
parser.add_argument('--file', help='fasta file input')

args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
input_file=args_dict['file'] #extract the file option value from dict

if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)
#print (input_file)
def main():
	try:
		with open(input_file, "r") as handle, open("{}.seqlengths".format(input_file), "w") as writehandle, open("stats.final.txt", "a") as writehandle2: # rU means open for reading using universal readline mode this means you dont have to worry if the file uses Unix, Mac or DOS Windows style newline characters The with- statement makes sure that the file is properly closed after reading it
			fasta_sequences = SeqIO.parse(handle,'fasta')
			writehandle2.write("\n\n\n\n\n\n\n\n\n\n\n\n\n")
			for entry in fasta_sequences:
				writehandle.write(entry.id)
				writehandle.write("\t")
				writehandle.write(str(len(str(entry.seq))))
				writehandle.write("\n")
				writehandle2.write(str(len(str(entry.seq))))
				writehandle2.write("\n")
			writehandle.close()
	except Exception  as e:
		print("Unexpected error:", str(sys.exc_info()))
		print("additional information:", e)
		PrintException()

##########################
## function definitions ##
##########################
def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

if __name__ == "__main__": main()
