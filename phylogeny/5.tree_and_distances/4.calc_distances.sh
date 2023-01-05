# we can calculate the all vs all distances using vsearch
vsearch --allpairs_global /Users/milanese/Dropbox/PhD/other/help_elisabetta_paper/4.do_MSA/all_seq.fna --id 0 --blast6out all_dist.tsv

# the result is in:
# 4.dist_all.tsv
