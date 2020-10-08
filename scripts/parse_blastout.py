#!/usr/bin/python3
from Bio import SeqIO
from Bio.Blast import NCBIXML
import argparse
import sys
import re

parser = argparse.ArgumentParser(description='please specify input file.')
parser.add_argument('--blastout', help='blastout xml file')
parser.add_argument('--fasta', help='blast input fasta file, for retrieving best match information')


args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
input_file=args_dict['blastout'] #extract the file option value from dict
input_fastafile=args_dict['fasta'] #extract the file option value from dict
#print (input_file)
genesPassFilter=dict()
gene_anno_pattern=re.compile('FAMILYNAMEPLACEHOLDER')
max_aln_num=2


if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)

def main():

	try:
	
		with open(input_file, "rU") as filehandle, open("{}.parsedNfilterout".format(input_file),'w') as writefilehandle: # rU means open for reading using universal readline mode this means you dont have to worry if the file uses Unix, Mac or DOS Windows style newline characters The with- statement makes sure that the file is properly closed after reading it=			
			
			writelog=open("log.parseblast",'w')
			blast_records = NCBIXML.parse(filehandle) # we have a pair of input functions, read and parse, where read is for when you have exactly one object, and parse is an iterator for when you can have lots of objects, but instead of getting SeqRecord or MultipleSeqAlignment objects, we get BLAST record objects.
			
			for blast_record in blast_records:	
				break_flag=0
				for alignment in blast_record.alignments:
					for hsp in alignment.hsps:
						#coverage=hsp.align_length/blast_record.query_length
						perc_iden=hsp.identities/hsp.align_length
						aligned_gene_name=alignment.title
						#print (aligned_gene_name)
						if not (re.search(gene_anno_pattern,aligned_gene_name) is None): # match gene_anno_pattern
							genesPassFilter[blast_record.query]=1						
							break
						else:
							writefilehandle.write(blast_record.query)
							writefilehandle.write("\n")
							writefilehandle.write(hsp.query)
							writefilehandle.write("\n")
							writefilehandle.write(aligned_gene_name)			
							writefilehandle.write("\n")
							break
					break_flag+=1
					if (break_flag>=max_aln_num):
						break
					
							

		
	except Exception  as e:
		print("Unexpected error:", str(sys.exc_info()))
		print("additional information:", e)

	try:
		with open(input_fastafile, "rU") as handle, open("{}.parsedNfiltered".format(input_file),'w') as writefilehandle2: # rU means open for reading using universal readline mode this means you dont have to worry if the file uses Unix, Mac or DOS Windows style newline characters The with- statement makes sure that the file is properly closed after reading it
			fasta_sequences = SeqIO.parse(handle,'fasta')
			for entry in fasta_sequences:
				name,seq,desc = entry.id, str(entry.seq), str(entry.description)
				#print (name)
				if name in genesPassFilter.keys(): # ortholog does not exists, print seq											
					writefilehandle2.write(">")
					writefilehandle2.write(name)
					writefilehandle2.write("\n")
					writefilehandle2.write(seq)
					writefilehandle2.write("\n")

				else:
					writelog.write("{} have orthologs, excluded from output\n".format(name))
					

		
	except Exception  as e:
		print("Unexpected error:", str(sys.exc_info()))
		print("additional information:", e)
		
if __name__ == "__main__": main()
