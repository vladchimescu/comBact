# since for all genomes we identified all the 40 MGs,
# here we do a concatenation of all genes.
# we will use this to calculate distances and create the tree.

echo ">s_typhimurium_LT2" > s_typhimurium_LT2.fna
grep -v ">" ../2.progenomes/s_typhimurium_LT2.fna >> s_typhimurium_LT2.fna

echo ">s_typhimurium_14028" > s_typhimurium_14028.fna
grep -v ">" ../2.progenomes/s_typhimurium_14028.fna >> s_typhimurium_14028.fna

echo ">s_pneumoniae_D39V" > s_pneumoniae_D39V.fna
grep -v ">" ../2.progenomes/s_pneumoniae_D39V.fna >> s_pneumoniae_D39V.fna

echo ">s_aureus_newman" > s_aureus_newman.fna
grep -v ">" ../2.progenomes/s_aureus_newman.fna >> s_aureus_newman.fna

echo ">s_aureus_DSM_20231" > s_aureus_DSM_20231.fna
grep -v ">" ../2.progenomes/s_aureus_DSM_20231.fna >> s_aureus_DSM_20231.fna

echo ">P_aeruginosa_PAO1" > P_aeruginosa_PAO1.fna
grep -v ">" ../2.progenomes/P_aeruginosa_PAO1.fna >> P_aeruginosa_PAO1.fna

echo ">p_aeruginosa_PA14" > p_aeruginosa_PA14.fna
grep -v ">" ../2.progenomes/p_aeruginosa_PA14.fna >> p_aeruginosa_PA14.fna

echo ">e_coli_K12_MG1655" > e_coli_K12_MG1655.fna
grep -v ">" ../2.progenomes/e_coli_K12_MG1655.fna >> e_coli_K12_MG1655.fna

echo ">e_coli_iAi" > e_coli_iAi.fna
grep -v ">" ../2.progenomes/e_coli_iAi.fna >> e_coli_iAi.fna

echo ">b_subtilis_168" > b_subtilis_168.fna
grep -v ">" ../2.progenomes/b_subtilis_168.fna >> b_subtilis_168.fna
