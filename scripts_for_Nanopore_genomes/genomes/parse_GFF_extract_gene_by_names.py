#!/usr/bin/python3
#this script extract GO terms associated with each AGAP identifier
from collections import defaultdict
import argparse
import sys
import linecache
import re
from Bio import SeqIO

##############
## arguments##
##############
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
parser= MyParser(description='this script extract genes with names that includes the keyword provided')
#parser = argparse.ArgumentParser()
parser.add_argument('-g', help='GFF3 file')
parser.add_argument('-k', help='keyword')
parser.add_argument('-f', help='feature type, e.g. gene')

args = parser.parse_args()
args_dict=vars(args) # convert namespace(args) to dict style
gff=args_dict['g'] #extract the file option value from dict
keyword=args_dict['k'] #extract the file option value from dict
keyword_nospace=re.sub("\s|\t","_",keyword)

feature_type=args_dict['f']
if len(sys.argv)==1: # print help message if arguments are not valid
    parser.print_help()
    sys.exit(1)


#####################
## initialization  ##
#####################	

## dict intialization

## defined pattern that needs to be searched in the main program
keyword_pattern=re.compile(keyword,re.IGNORECASE)
feature_type_pattern = re.compile(feature_type)
ID_pattern = re.compile('ID=(.+?)(?:\.|:|;|$)')
NAME_pattern = re.compile('Name=(.+?)(;|$)',re.IGNORECASE)

#####################
##      main       ##
#####################
def main():
	try: 
		## go through gff file
		with open(gff, "rU") as handle, open("{}.{}.{}.tab".format(gff,feature_type,keyword_nospace),"w") as writehandle, open("{}.{}.{}.gff".format(gff,feature_type,keyword_nospace),"w") as writegffhandle : 
			writehandle.write("chr\tstart\tend\tstrand\tgeneID\tngeneName\n")
			for line_raw in handle:
				if (not (re.search(feature_type_pattern,line_raw) is None)):
					if ( not (re.search(keyword_pattern,line_raw) is None)):
						#print (line_raw)
						line_raw=line_raw.strip()
						m=re.search(ID_pattern,line_raw)
						m2=re.search(NAME_pattern, line_raw)					
						gene_ID='N/A'
						NAME = 'N/A'
						if m and m2: 
							gene_ID=m.group(1)								
							NAME=m2.group(1)
							#print(line_raw)
							#print(gene_ID)
							#print(NAME)
							fields=line_raw.split("\t")
							chr=fields[0]
							start=fields[3]
							end=fields[4]
							strand=fields[6]
							#print("{} {} {} {}\n".format(chr,start,end,strand))
							writehandle.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,start,end,strand,gene_ID,NAME))
							writegffhandle.write("{}\n".format(line_raw.strip()))
						else:
							print("error getting ID or name from line:{}".format(line_raw))
						#gene_ID=re.sub("-..$",'',gene_ID)

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

