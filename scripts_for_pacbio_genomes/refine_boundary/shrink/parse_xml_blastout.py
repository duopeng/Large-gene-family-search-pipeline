#!/usr/bin/python3
from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast import NCBIXML
import argparse
import sys
import re
import linecache

parser = argparse.ArgumentParser(description='please specify input file.')
parser.add_argument('--blastout', help='blastout xml file')

args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
input_file=args_dict['blastout'] #extract the file option value from dict

if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)

def main():

	try:
		with open("{}.parsed.tab".format(input_file),'w') as writefilehandle: 			
			writefilehandle.write("query_id\tquery_length\tquery_start\tquery_end\tquery_frame\tquery_coverage\thit_acc\thit_desc\thit_id\thit_length\thit_start\thit_end\thit_frame\thit_coverage\tE_val\tpercent_identity\tgap_num\n")#write header to parsed output file
			blast_Qresults = SearchIO.parse(input_file, 'blast-xml') #open blast output file , format is 'blast-xml'
			for QueryResults in blast_Qresults:	
				query_id=QueryResults.id
				query_len=QueryResults.seq_len
				for hit in QueryResults:
					hit_acc=hit.accession
					hit_desc=hit.description
					hit_id=hit.id
					hit_len=hit.seq_len
					for hsp in hit:
						hsp_eval=hsp.evalue
						hsp_identicalnum=hsp.ident_num
						hsp_gapnum=hsp.gap_num
						hsp_hitstart=hsp.hit_start
						hsp_hitend=hsp.hit_end
						hsp_querystart=hsp.query_start
						hsp_queryend=hsp.query_end
						hsp_hitframe=hsp.hit_frame
						hsp_queryframe=hsp.query_frame
						hsp_hitseq=hsp.hit
						hsp_alnlen=hsp.aln_span
						percent_ident=100*hsp_identicalnum/hsp_alnlen
						query_coverage=(int(hsp_queryend)-int(hsp_querystart))/int(query_len)
						hit_coverage=(int(hsp_hitend)-int(hsp_hitstart))/int(hit_len)
						writefilehandle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(query_id,query_len,hsp_querystart,hsp_queryend,hsp_queryframe,str(query_coverage),hit_acc,hit_desc,hit_id,hit_len,hsp_hitstart,hsp_hitend,hsp_hitframe,str(hit_coverage),hsp_eval,str(percent_ident),hsp_gapnum))#write current hsp to output file
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
