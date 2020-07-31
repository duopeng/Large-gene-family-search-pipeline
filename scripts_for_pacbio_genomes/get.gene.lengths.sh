declare -a file_array=( 
Pred.DGF-1.Brazil.boundaryRefined.fasta
Pred.DGF-1.Y.boundaryRefined.fasta
Pred.GP63.Brazil.boundaryRefined.fasta
Pred.GP63.Y.boundaryRefined.fasta
Pred.MASPs.Brazil.boundaryRefined.fasta
Pred.MASPs.Y.boundaryRefined.fasta
Pred.Mucins.Brazil.boundaryRefined.fasta
Pred.Mucins.Y.boundaryRefined.fasta
Pred.RHS.Brazil.boundaryRefined.fasta
Pred.RHS.Y.boundaryRefined.fasta
Pred.ts.Brazil.boundaryRefined.fasta
Pred.ts.Y.boundaryRefined.fasta
)

for i in "${file_array[@]}"
do
   
   echo "python3 Get_seq_lengths_from_fasta.v1.1.py --file ${i}"
   python3 Get_seq_lengths_from_fasta.v1.1.py --file $i
done
