#############################################
#change reverse blast gene name filter!!
#############################################

genomedir=./genomes
perc_identity_cutoff=85
mapHspGap=0
minTcTSLength=150

geneFamilyInput=./all_and_full
genefamilygff=dummy.gff

#parameter setting [case to case bases]
familyname=familynameplaceholder
genome_name=genomenameplaceholder #Brazil or Y
all_or_full=all_or_full_placeholder #full or all
#end of parameter setting


genefamilyfasta=''
if [ ${geneFamilyInput} == "full" ]; then
	genefamilyfasta=${geneFamilyInput}/${familyname}/Full_ClBrener_${familyname}_DNA_sequence.fasta
fi
if [ ${geneFamilyInput} == "all" ]; then
	genefamilyfasta=${geneFamilyInput}/${familyname}/All_ClBrener_${familyname}_DNA_sequence.fasta	
fi	
if [ ${geneFamilyInput} == "WW" ]; then
	genefamilyfasta=${geneFamilyInput}/${familyname}/WW_${genome_name}_${familyname}_DNA_sequence.fasta	
fi	
if [ ${geneFamilyInput} == "40ts" ]; then
	genefamilyfasta=${geneFamilyInput}/${familyname}/40ts_${familyname}_DNA_sequence.fasta	
fi	

genomefasta=${genome_name}.fasta  


genefamilydir=./family_search_${familyname}_${genome_name}_${geneFamilyInput}
codedir=${genefamilydir}
cd ${genefamilydir}

#blast
/home/duopeng/bin/x86_64/ncbi-blast-2.8.1+/bin/blastn \
-db ${genomedir}/${genomefasta} \
-query ${genefamilyfasta} \
-out ${perc_identity_cutoff}percid.2.8.1.blastout \
-num_threads 1 \
-num_alignments 100 \
-max_hsps 100 \
-perc_identity ${perc_identity_cutoff}

#parse blast, get gene candidates
perl ${codedir}/mark_homology_pieces.pl \
-i ${perc_identity_cutoff}percid.2.8.1.blastout \
-g ${genomedir}/${genomefasta} -maxhspgap \
${mapHspGap} -minTcTSLength ${minTcTSLength} \
-gff ${genefamilygff} -familyfasta ${genefamilyfasta} \
> log_of_gene_num.txt

rm ${perc_identity_cutoff}percid.2.8.1.blastout


#blast back to all transcripts
/home/duopeng/bin/x86_64/ncbi-blast-2.8.1+/bin/blastn \
-outfmt 5 \
-db /home/duopeng/TcTS_transcripts/TriTrypDB-34_TcruziCLBrener_AnnotatedTranscripts.fasta \
-query ${genefamilydir}/${perc_identity_cutoff}percid.2.8.1.blastout.${genomefasta}.mintctslengthcutoff${minTcTSLength}.maxgaplenintcts${mapHspGap}.fasta \
-out ${perc_identity_cutoff}perc_${minTcTSLength}minlen_${mapHspGap}maxgap.blast2transcripts.blastoutxml -num_threads 1 \
-num_alignments 50 \
-max_hsps 50

#parse blast
python3 parse_blastout.edited.py \
--blastout ${perc_identity_cutoff}perc_${minTcTSLength}minlen_${mapHspGap}maxgap.blast2transcripts.blastoutxml \
--fasta ${genefamilydir}/${perc_identity_cutoff}percid.2.8.1.blastout.${genomefasta}.mintctslengthcutoff${minTcTSLength}.maxgaplenintcts${mapHspGap}.fasta \

rm ${perc_identity_cutoff}perc_${minTcTSLength}minlen_${mapHspGap}maxgap.blast2transcripts.blastoutxml

#adjust boundary
python3 adjust_start_stop_codon.py \
--fastafile 85perc_150minlen_0maxgap.blast2transcripts.blastoutxml.parsedNfiltered

#get length distribution
python3 Get_seq_lengths_from_fasta.v1.1.py \
--file 85perc_150minlen_0maxgap.blast2transcripts.blastoutxml.parsedNfiltered.adjusted.fasta

