#!/usr/bin/python3
from collections import defaultdict
import argparse
import sys
import linecache
import re
from Bio import SeqIO
import collections
import subprocess
##############
## arguments##
##############
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
parser= MyParser(description='This script does xxx')
#parser = argparse.ArgumentParser()
parser.add_argument('--tab', help='table of parsed blast out')
parser.add_argument('--fasta', help='gene seqeuences in fasta format,')
parser.add_argument('--dist', help='will not update ts boundary if match region is beyond [dist],')
args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
tabfile=args_dict['tab'] #extract the file option value from dict
fastafile=args_dict['fasta'] #extract the file option value from dict
dist=int(args_dict['dist'])
if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)


#####################
## initialization  ##
#####################	

## dict intialization
gene_boundary=defaultdict(dict) ## initialized a multi-dimension dict,
gene_boundary_thres=defaultdict(dict) 
gene_length=defaultdict(dict) ## initialized a multi-dimension dict
gene_bounds_reached_flags=defaultdict(dict)
gene_seqs=defaultdict(dict)

#####################
##      main       ##
#####################	
def main():

## go through a text file
	try: 
		with open(tabfile, "rU") as handle: 
				next(handle) #skip header
				for line_raw in handle:
					lineSplitArray=re.split("\t",line_raw.strip())
					q_id=lineSplitArray[0]
					q_len=int(lineSplitArray[1])
					q_match_start=int(lineSplitArray[2])
					q_match_end=int(lineSplitArray[3])
					#print(q_id)
					if not q_id in gene_length:
						gene_length[q_id]=q_len

##update boundary not thresholded
					#find the earliest start and latest end match point
					if q_id in gene_boundary:
						if q_match_start < gene_boundary[q_id]['start']:
							gene_boundary[q_id]['start']=q_match_start

						if q_match_end > gene_boundary[q_id]['end']:
							#print ("{} > {}".format(q_match_end,gene_boundary[q_id]['end']))
							gene_boundary[q_id]['end']=q_match_end

					else:
						gene_boundary[q_id]['start']=q_match_start
						gene_boundary[q_id]['end']=q_match_end

##update boundary_thresholded
					for q_id in gene_boundary.keys():
						q_len=gene_length[q_id]
						q_match_start_earliest=int(gene_boundary[q_id]['start'])
						q_match_end_latest=int(gene_boundary[q_id]['end'])
						
						if q_match_start_earliest<=dist:
							gene_boundary_thres[q_id]['start']=q_match_start_earliest
						else:
							gene_boundary_thres[q_id]['start']=0
							
						if (q_len - q_match_end_latest)<=dist:
							gene_boundary_thres[q_id]['end']=q_match_end_latest
						else:
							gene_boundary_thres[q_id]['end']=q_len

					#print ("{},{},{},{},{},{},{},{},{},{},{}\n".format(q_id,q_len,q_match_start,q_match_end,gene_boundary[q_id]['start'],gene_boundary[q_id]['end'],gene_boundary_thres[q_id]['start'],gene_boundary_thres[q_id]['end']))
	except Exception  as e:
		print("Unexpected error:", str(sys.exc_info()))
		print("additional information:", e)
		PrintException()

## print out list of original coordinates and modified coordinates
	try:
		with open("{}.oldNnew.coords".format(fastafile.split("/")[-1]), "w") as whandle:
			for id in gene_length.keys():
				current_gene_length=gene_length[id]
				refined_start=gene_boundary[id]['start']
				refined_end=gene_boundary[id]['end']
				refined_start_thres=gene_boundary_thres[id]['start']
				refined_end_thres=gene_boundary_thres[id]['end']
				whandle.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(id,current_gene_length,refined_start,refined_end,refined_start_thres,refined_end_thres))
				

	except Exception  as e:
		print("Unexpected error:", str(sys.exc_info()))
		print("additional information:", e)
		PrintException()

## print out sequence
	#load TcTS sequences
	with open(fastafile, "rU") as handle:
		fasta_sequences = SeqIO.parse(handle,'fasta')
		for entry in fasta_sequences:
			name, desc,seq = entry.id, entry.description, str(entry.seq)
			gene_seqs[name]=seq

	##export TcTS sequences##
	try:
		with open(fastafile, "rU") as handle, open("{}.shrk".format(fastafile.split("/")[-1]),'w') as writehandle:
			fasta_sequences = SeqIO.parse(handle,'fasta')
			for entry in fasta_sequences:
				genename, desc = entry.id, entry.description
				print(genename)
				if genename in gene_boundary_thres:
					new_start=int(gene_boundary_thres[genename]['start'])
					new_end=int(gene_boundary_thres[genename]['end'])
				else:
					new_start=0
					new_end=len(gene_seqs[genename])
				seq=gene_seqs[genename][(new_start):(new_end)]
				writehandle.write(">{}_{}_{}\n{}\n".format(genename,new_start,new_end,seq))
				
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

