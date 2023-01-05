# we do here the multiple sequence alignment with muscle

cat ../3.concat_genes/*.fna > all_seq.fna
muscle -in all_seq.fna -out all_seq.ali
# we change the type of alignment file
esl-reformat phylip all_seq.ali > all_seq.phylip
