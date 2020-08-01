genomefile=Brazil.fasta
gfffile=Brazil_0201.gff3
keyword="galactofuranosyl glycosyltransferase"
keyword_noblank=${keyword// /_}
featuretype="gene"
python3 parse_GFF_extract_gene_by_names.py -g ${gfffile}  -k "${keyword}" -f ${featuretype}
python3 Get_seq_by_coord_from_fasta.py --fafile ${genomefile} --coordfile ${gfffile}.${featuretype}.${keyword_noblank}.tab 



genomefile=Y.fasta
gfffile=Y_all_curated_annotated.gff
keyword="galactofuranosyl glycosyltransferase"
keyword_noblank=${keyword// /_}
featuretype="gene"
python3 parse_GFF_extract_gene_by_names.py -g ${gfffile}  -k "${keyword}" -f ${featuretype}
python3 Get_seq_by_coord_from_fasta.py --fafile ${genomefile} --coordfile ${gfffile}.${featuretype}.${keyword_noblank}.tab 
