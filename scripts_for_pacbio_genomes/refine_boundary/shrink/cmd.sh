perc_identity_cutoff=85
family_name=Mucins
blast_db_1=fasta/${family_name}_models.fasta # the model ts genes
blast_query_1=/home/duopeng/gene_member_search/refine_boundary_${family_name}/refine_boundary_extend/Pred.${family_name}.Brazil.fasta.ext 
blast_db_2=fasta/${family_name}_models.fasta # genome
blast_query_2=/home/duopeng/gene_member_search/refine_boundary_${family_name}/refine_boundary_extend/Pred.${family_name}.Y.fasta.ext # genome

/home/duopeng/bin/x86_64/ncbi-blast-2.8.1+/bin/makeblastdb -dbtype nucl -in  ${blast_db_1}
/home/duopeng/bin/x86_64/ncbi-blast-2.8.1+/bin/makeblastdb -dbtype nucl -in ${blast_db_2}

echo "blasting..."
/home/duopeng/bin/x86_64/ncbi-blast-2.8.1+/bin/blastn -outfmt 5 -db ${blast_db_1} -query ${blast_query_1} -out ${blast_query_1##*\/}-${blast_db_1##*\/}.blastout.xml -num_threads 1 -num_alignments 2000 -max_hsps 2000 -perc_identity ${perc_identity_cutoff}
echo "done"
echo "blasting..."
/home/duopeng/bin/x86_64/ncbi-blast-2.8.1+/bin/blastn -outfmt 5 -db ${blast_db_2} -query ${blast_query_2} -out ${blast_query_2##*\/}-${blast_db_2##*\/}.blastout.xml -num_threads 1 -num_alignments 2000 -max_hsps 2000 -perc_identity ${perc_identity_cutoff}
echo "done"

#${blast_db_1  <-- from variable blast_db_1
#  ##   <-- greedy front trim
#  *    <-- matches anything
#  :    <-- until the last ':'
 #}
 
#parse blast
python3 parse_xml_blastout.py --blastout ${blast_query_2##*\/}-${blast_db_2##*\/}.blastout.xml
python3 parse_xml_blastout.py --blastout ${blast_query_1##*\/}-${blast_db_1##*\/}.blastout.xml

#refine boundary
python3 refine_boundary_based_on_Blast.py --tab ${blast_query_2##*\/}-${blast_db_2##*\/}.blastout.xml.parsed.tab --fasta ${blast_query_2} --dist 100
python3 refine_boundary_based_on_Blast.py --tab ${blast_query_1##*\/}-${blast_db_1##*\/}.blastout.xml.parsed.tab --fasta ${blast_query_1}  --dist 100