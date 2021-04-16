# Pipeline that search for members of large gene families in Pacbio/Nanopore genomes

## Prerequisites
Linux-like environment

Two executable programs "blastn" and "makeblastdb" from the ncbi-blast suite (tested with version 2.8.1+)

Python3 (tested with version 3.7.0)

PERL (tested with v5.18.2)

Bioperl (tested with v1.006924)


## Getting Started
The scripts folder contains scripts and input files for a working example, in which the the T cruzi genome is searched for mucin gene family members

To start the analysis, first configure parameters in the shell script file "start_analysis.sh". Then you can either call your shell to run the file or execute line-by-line in your shell command line

##### Folder Structure:
```
.
├── scripts             # scripts used to do family search
│     │
│     └──input_files    #input files for a working example
│
└── README.md
