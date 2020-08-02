#!/usr/bin/python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import argparse
import sys
from collections import defaultdict
import re
import linecache
import os
import shutil
import subprocess

##############
## arguments##
##############

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
parser= MyParser(description='This script ')
#parser = argparse.ArgumentParser()
parser.add_argument('--listfile', help='3 columns, gene family, genome, full/all')
args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
listfile=args_dict['listfile'] #extract the file option value from dict
if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)
	
if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)
#print (input_file)

genes=dict()
genomes=dict()
all_fulls=dict()

ID_pattern = re.compile('ID=(.+?)(?:\.|:|;|$)')
ID_pattern = re.compile('ID=(.+?)(?:\.|:|;|$)')
ID_pattern = re.compile('ID=(.+?)(?:\.|:|;|$)')
		
def main():

	try:
	#start filtering fasta file
		with open(listfile, "rU") as handle: # rU means open for reading using universal readline mode this means you dont have to worry if the file uses Unix, Mac or DOS Windows style newline characters The with- statement makes sure that the file is properly closed after reading it
			
			##read in gene name, genomes and create working dir and copy commands, modify command
			for line in handle:
				fields=line.strip().split("\t")
				try:
					if fields[0] != '':
						genes[fields[0]]=fields[0]
				except:
					pass
				try:
					if fields[1] != '':
						genomes[fields[1]]=fields[1]
				except:
					pass
				try:
					if fields[2] != '':
						all_fulls[fields[2]]=fields[2]
				except:
					pass
			for gene in genes.keys():
				for genome in genomes.keys():
					for all_full in all_fulls.keys():
						newpath="family_search_{}_{}_{}".format(gene,genome,all_full)
						if (os.path.isdir(newpath) and os.path.exists(newpath)): #delete existing dir
							shutil.rmtree(newpath)
						os.mkdir(newpath) #create new dir
						##copy commands
						subprocess.call("cp scripts/* {}/".format(newpath),shell=True )
						##change blastcmd
						with open("{}/1.blast.cmd.sh".format(newpath), "rt") as fin:
							with open("{}/1.blast.cmd.edited.sh".format(newpath), "wt") as fout:
								for line in fin:
									line=line.replace('familynameplaceholder', gene)
									line=line.replace('genomenameplaceholder', genome)
									line=line.replace('all_or_full_placeholder', all_full)
									fout.write(line)

						#change family name in reverse parse blastout
						with open("{}/parse_blastout.py".format(newpath), "rt") as fin:
							with open("{}/parse_blastout.edited.py".format(newpath), "wt") as fout:
								for line in fin:
									if gene == "beta-gal_GT":
										line=line.replace('FAMILYNAMEPLACEHOLDER', "beta galactofuranosyl glycosyltransferase")
									elif gene == "mucins":
										line=line.replace('FAMILYNAMEPLACEHOLDER', "mucin ")
									elif gene == "MASPs":
										line=line.replace('FAMILYNAMEPLACEHOLDER', "MASP")										
									else:
										line=line.replace('FAMILYNAMEPLACEHOLDER', gene)
									fout.write(line)
						
						#start family search
						subprocess.call("nohup bash {}/1.blast.cmd.edited.sh >{}/nohup.log &".format(newpath,newpath),shell=True )
			
			
			
			
			
		
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