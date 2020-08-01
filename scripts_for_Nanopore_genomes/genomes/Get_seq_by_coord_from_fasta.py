#!/usr/bin/python3
from collections import defaultdict
import argparse
import sys
import linecache
import re
from Bio import SeqIO
from Bio import SeqFeature

##############
## arguments##
##############
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
parser= MyParser(description='This script fetches seqs by coordinates from a fasta file.  please specify all optional parameters. ')
#parser = argparse.ArgumentParser()
parser.add_argument('--fafile', help='fasta file')
parser.add_argument('--coordfile', help='coordinate file, first 4 columns must be chr, start, end, strand, [gene ID], [gene name]')
args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
fafile=args_dict['fafile'] #extract the file option value from dict
coordfile=args_dict['coordfile'] #extract the file option value from dict
if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)


#####################
## initialization  ##
#####################	

## dict intialization
coord_dict=defaultdict(dict) ## initialized a multi-dimension dict

#####################
##      main       ##
#####################	
def main():
	try: 
		## go through a text file
		with open(coordfile, "rU") as handle: # rU means open for reading using universal readline mode this means you dont have to worry if the file uses Unix, Mac or DOS Windows style newline characters The with- statement makes sure that the file is properly closed after reading it
			i=1 #counter initialization
			next(handle)#skip header
			for line_raw in handle:
				fields=line_raw.strip().split("\t")
				chr=fields[0]
				start=fields[1]
				end=fields[2]
				strand=fields[3]
				geneID=""
				geneName=""
				if 4 < len(fields): geneID=fields[4] 
				if 5 < len(fields): geneName=fields[5] 
				coord_dict[i]['chr']=chr
				coord_dict[i]['start']=start
				coord_dict[i]['end']=end
				coord_dict[i]['strand']=strand
				coord_dict[i]['geneID']=geneID
				coord_dict[i]['geneName']=geneName
				#print("{}\t{}\t{}\t{}\t{}\t{}".format(chr,start,end,strand,geneID,geneName))
				i+=1
		## go through fasta file
		with open(fafile, "rU") as handle, open("{}.{}.fa".format(fafile,coordfile),"w") as wh:
			fasta_sequences = SeqIO.parse(handle,'fasta')
			for entry in fasta_sequences:
				name, desc, seq = entry.id, entry.description, str(entry.seq)
				#print(name)
				#print(desc)
				for i in coord_dict.keys():
					chr = coord_dict[i]['chr']
					if name == chr:
						#print("start to extract genes")
						start=coord_dict[i]['start']
						end=coord_dict[i]['end']
						strand=coord_dict[i]['strand']
						geneID=coord_dict[i]['geneID']
						geneName=coord_dict[i]['geneName']						
						start_pos = SeqFeature.ExactPosition(start)
						end_pos = SeqFeature.ExactPosition(end)
						if strand == "+": # .extract method only takes 1 or -1 strands
							strand=1
						if strand == "-":
							strand = -1
							
						start_pos = SeqFeature.ExactPosition(int(start)-1) #later extract method will use slice[start:end], which result in start-1(?? now sure how but emprically determined to be so) 
						gene_location=SeqFeature.FeatureLocation(start_pos, end_pos, strand=strand)
						geneseq=gene_location.extract(seq)
						#print("{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,start,end,strand,geneID,geneName))
						#print(geneseq)
						wh.write(">{}\t{}\t{}\t{}\t{}\t{}\n{}\n".format(geneID,geneName,chr,start,end,strand,geneseq))
				#flag=0
				#for ID in ID_dict:
				#	regex_pattern = r"\b(?=\w)" + re.escape(ID) + r"\b(?!\w)" #build regex pattern
				#	if not (re.search(regex_pattern,name) is None):
				#		flag=1
				#if flag==1:
				#	wh.write(">{}\t{}\n{}\n".format(name,desc,seq))
	
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







