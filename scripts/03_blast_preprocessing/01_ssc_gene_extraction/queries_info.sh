#count of clustering_ids per species
awk '/^>/{sp=$NF; c[sp]++} END{for (s in c) printf "%s\t%d\n", s, c[s]}' \
  /data/pam/team230/sm71/scratch/rp2/run_blast/queries/queries.fa \
  | sort \
  > /data/pam/team230/sm71/scratch/rp2/run_blast/queries/queries_per_species_counts.tsv

#Enterococcus_faecalis   10138636
#Enterococcus_faecium    1696917
#Staphylococcus_argenteus        55017
#Staphylococcus_aureus   1768507
#Staphylococcus_capitis  349773
#Staphylococcus_epidermidis      149985
#Staphylococcus_haemolyticus     795818
#Staphylococcus_pseudintermedius 3166426
#Staphylococcus_sciuri   2334454
#Streptococcus_agalactiae        9945383
#Streptococcus_dysgalactiae      204462
#Streptococcus_equi      276716
#Streptococcus_mitis     105671
#Streptococcus_mutans    1005044
#Streptococcus_pneumoniae        3607001
#Streptococcus_pyogenes  1229927
#Streptococcus_suis      3934755
#Streptococcus_uberis    596099

