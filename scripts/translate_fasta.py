#!/usr/bin/python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import argparse
import sys
from collections import defaultdict
import re
import linecache



##############
## arguments##
##############

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
parser= MyParser(description='This script translates a fasta file')
#parser = argparse.ArgumentParser()
parser.add_argument('--fastafile', help='fasta file input')
parser.add_argument('--lenCutoff', help='minium peptide length, peptides below this cutoff will be discarded')
args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
input_fastafile=args_dict['fastafile'] #extract the file option value from dict
lengthcutoff=args_dict['lenCutoff'] #extract the file option value from dict
if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)
	
if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)
#print (input_file)


startcodonpattern=re.compile("^ATG", flags=re.IGNORECASE)
	
		
def main():

	try:
	#start filtering fasta file
		with open(input_fastafile, "rU") as handle, open("{}.translation.{}lenCutoff.fasta".format(input_fastafile,lengthcutoff),'w') as writehandle: # rU means open for reading using universal readline mode this means you dont have to worry if the file uses Unix, Mac or DOS Windows style newline characters The with- statement makes sure that the file is properly closed after reading it
			fasta_sequences = SeqIO.parse(handle,'fasta')
			for entry in fasta_sequences:
				name,seq, desc = entry.id, str(entry.seq), entry.description
				coding_dna=Seq(seq,generic_dna)
				if (not (re.search(startcodonpattern,str(coding_dna)) is None)):
					#print(coding_dna)
					seq_translation=str(coding_dna.translate(to_stop=True))
					if len(seq_translation)>=int(lengthcutoff):
						#print(name)
						writehandle.write(">{}\n{}\n".format(name,seq_translation))
		
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