# we used ../4.do_MSA/all_seq.ali in
# http://www.atgc-montpellier.fr/job_status/index.php?jobid=2640537&path=20210303-182710_Bv13/cmd.txt

# (PhyML 3.0)

#Note that I used default parameters except for doing 100 bootstraps and do not
# use the fast alignment methods (so that it should be more precise)

# The result is in the dir:
#  all_seq_phylip_phyml


#-------------------------------------------------------------------------------
# Second tree:
# We also try to create a tree with fasttree.
# But the one created above should be of better quality!

# run in dir 4:
fasttree -nt all_seq.ali > ../5.tree_and_distances/fasttree.tree
