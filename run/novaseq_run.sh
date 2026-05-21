#!/usr/bin/env bash

export PATH=/envs/git_env/bin:$PATH

set -e

cd /work

# filter and trim using cutadapt and dada2. ASV inference using dada2.
Rscript novaseq/scripts/project_scripts/loessErrfun_mod4_sol.R -d /envs/git_env/bin/cutadapt
echo "-----------------filter and trim done-----------------"
break

# remove chimeras and singletons
Rscript novaseq/scripts/original_scripts/18S_chimera_singleton_removal_taxonomic_classification.R
echo "-----------------chimera and singleton removal + taxonomic classification done-----------------"

# generate fasta file, count table and ensemble taxonomy for 18S ASVs
Rscript novaseq/scripts/original_scripts/Ensemble_taxonomy_18S_ASVs.R

# negative control correction
Rscript novaseq/scripts/original_scripts/18S_blank_corr.R
echo "-----------------negative control correction done-----------------"

# subset the ASVs remaining after negative control correction
grep -w -A 1 -f novaseq/18S/blank_corr/no_contam_headers_18S.txt novaseq/18S/18S_nochim_nosingle_ASVs.fa  --no-group-separator > novaseq/18S/blank_corr/18S_nochim_nosingle_nocontam.fa 

# generate headers with the abundance of each ASV included
Rscript novaseq/scripts/original_scripts/dereplication_headers.R
echo "-----------------dereplication headers done-----------------"

# replace headers with dereplicated header names
grep -v "^--" novaseq/COI/blank_corr/COI_nochim_nosingle_nopseudo_nocontam.fa | awk 'NR%2==0' | paste -d'\n' novaseq/COI/ASV_dereplicated.txt - > novaseq/COI/COI_dereplicated_ASVs.fa

# change directory for swarm to run 
pushd novaseq/COI/

# cluster ASVs into MOTUs using swarm
/envs/git_env/bin/swarm -d 13 -i swarm/internal.txt -o swarm/output.txt -s swarm/statistics.txt -u swarm/uclust.txt -w swarm/COI_cluster_reps.fa COI_dereplicated_ASVs.fa
echo "-----------------swarm done-----------------"

# move back to root
popd

# generate MOTU tables from swarm output 
Rscript novaseq/scripts/original_scripts/MOTU_tables.R

# remove read abundance line from sequence header
awk -F'_' '{print $1}' novaseq/COI/swarm/COI_cluster_reps.fa > novaseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa

# generate match lists using BLASTn
makeblastdb -in novaseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa -parse_seqids -dbtype nucl
blastn -db novaseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa -outfmt '6 qseqid sseqid pident' -out novaseq/COI/MOTU/match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query novaseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa

# LULU curation
Rscript novaseq/scripts/original_scripts/LULU_curation.R
echo "-----------------LULU curation done-----------------"

# generate fasta files with remaining MOTUs after LULU curation
grep -w -A 1 -f novaseq/COI/MOTU/lulu_curated_headers.txt novaseq/COI/MOTU/COI_cluster_reps_lulu_ready.fa --no-group-separator > novaseq/COI/MOTU/COI_cluster_reps_lulu_curated.fa

# taxonomic assignment using BOLDigger3 on public animal library (--db 1) on exhaustive search mode (--mode 3)
pip install boldigger3==2.1.4
pip install lxml-html-clean==0.4.3
boldigger3 identify novaseq/COI/MOTU/COI_cluster_reps_lulu_curated.fa --db 1 --mode 3

echo "-----------------BOLDigger done-----------------"

# create taxonomic table from boldigger output 
Rscript novaseq/scripts/project_scripts/BOLDigger_tax_table.R

echo "-----------------script finished successfully-----------------"