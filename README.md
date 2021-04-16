# Pipeline that search for members of large gene families in Pacbio/Nanopore genomes

## Prerequisites
Linux-like environment

Two executable programs "blastn" and "makeblastdb" from the ncbi-blast suite (tested with version 2.8.1+)

Python3 (tested with version 3.8.5)

PERL (tested with v5.18.2)

Bioperl (tested with v1.006924)


## Getting Started
The scripts folder contains scripts and input files for a working example, in which the the T cruzi genome is searched for mucin gene family members


To start the example search, you have to first configure the following two parameter in the shell script file "start_search.sh"

BLASTN_path=""      #path to blastn program

makeblastdb_path="" #path to makeblastdb program

then execute the shell script file "start_search.sh". For example "bash start_search.sh"


To start your custom search, please configure input file parameters (lines 8-11) in the shell script file "start_analysis.sh" (using a text editor).

Next, open script file "parse_blastout.py" and edit line 19. Change the gene name that is between the two quotation marks.

Then you can either execute the shell script file "start_search.sh" or execute line-by-line in your shell command line

After running all steps, you should see a file named "result.fasta" which contains gene family members found.

##### Folder Structure:
```
.
├── scripts             # scripts used to do family search
│     │
│     └──input_files    #input files for a working example
│
└── README.md
