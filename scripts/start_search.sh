#############################################
#Plese check these parameter

# files and inputs:
BLASTN_path="/home/duopeng/bin/x86_64/ncbi-blast-2.8.1+/bin/blastn" #path to blastn program
makeblastdb_path="/home/duopeng/bin/x86_64/ncbi-blast-2.8.1+/bin/makeblastdb" #path to blastn program

genomefasta="input_files/Tcruzi_Y.fasta" # genome to perform search in (fasta format)
genefamilyfasta='input_files/mucins_DNA_sequence.fasta' # a collection of family member genes
all_transcripts="input_files/TriTrypDB-34_TcruziCLBrener_AnnotatedTranscripts.fasta" #All transcritps in the genome
genefamilygff="input_files/empty.gff" # coordinates of known family member genes, use empty.gff if you don't want to provide this
#
#Cut-offs:
perc_identity_cutoff=85
mapHspGap=0
minTcTSLength=150
#############################################

#make blast database
printf "\n\nmaking blast database from fasta file...\n"
${makeblastdb_path} -dbtype nucl -in ${genomefasta}

#blast input genes to the genome
printf "\n\nBLASTing input genes to the genome..."
${BLASTN_path} -db ${genomefasta} -query ${genefamilyfasta} -out percid_${perc_identity_cutoff}.blastout -num_threads 1 -num_alignments 100 -max_hsps 100 -perc_identity ${perc_identity_cutoff}

#parse blast, get gene candidates
printf "\n\nparsing BLAST results, obtaining gene candidates..."
perl mark_homology_pieces.pl -i percid_${perc_identity_cutoff}.blastout -g ${genomefasta} -maxhspgap ${mapHspGap} -minTcTSLength ${minTcTSLength} -gff ${genefamilygff} -familyfasta ${genefamilyfasta}  > log_of_gene_num.txt
rm percid_${perc_identity_cutoff}.blastout

#blast back to all transcripts
printf "\n\nmaking blast database from fasta file...\n"
${makeblastdb_path} -dbtype nucl -in ${all_transcripts}
printf "\n\nBLASTing gene candidates to all transcripts in the genome..."
${BLASTN_path} -outfmt 5 -db ${all_transcripts} -query percid_${perc_identity_cutoff}.blastout.mintctslengthcutoff${minTcTSLength}.maxgaplenintcts${mapHspGap}.fasta -out ${perc_identity_cutoff}perc_${minTcTSLength}minlen_${mapHspGap}maxgap.blast2transcripts.blastoutxml -num_threads 1 -num_alignments 50 -max_hsps 50

#parse blast
printf "\n\nparsing BLAST results, filtering gene candidates..."
python3 parse_blastout.py --blastout ${perc_identity_cutoff}perc_${minTcTSLength}minlen_${mapHspGap}maxgap.blast2transcripts.blastoutxml --fasta percid_${perc_identity_cutoff}.blastout.mintctslengthcutoff${minTcTSLength}.maxgaplenintcts${mapHspGap}.fasta
rm ${perc_identity_cutoff}perc_${minTcTSLength}minlen_${mapHspGap}maxgap.blast2transcripts.blastoutxml

#adjust boundary
printf "\n\nadjusting gene boundaries..."
python3 adjust_start_stop_codon.py --fastafile 85perc_150minlen_0maxgap.blast2transcripts.blastoutxml.parsedNfiltered

printf '\n\ndone\n'
