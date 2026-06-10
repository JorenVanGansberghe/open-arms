## These are the steps of the pipeline:

### Filter and trim using cutadapt 

`Rscript novaseq/scripts/project_scripts/loessErrfun_mod4_sol.R -d CUTADAPT_DIR`

### Remove chimeras and singletons and taxonomic classification of ASVs

`Rscript novaseq/scripts/original_scripts/18S_chimera_singleton_removal_taxonomic_classification.R`

### Generate fasta file, count table and ensemble taxonomy for 18S ASVs

`Rscript novaseq/scripts/original_scripts/Ensemble_taxonomy_18S_ASVs.R`

### Blank correction 

`Rscript novaseq/scripts/original_scripts/18S_blank_corr.R`

`grep -w -A 1 -f novaseq/18S/blank_corr/no_contam_headers_18S.txt novaseq/18S/18S_nochim_nosingle_ASVs.fa  --no-group-separator > novaseq/18S/blank_corr/18S_nochim_nosingle_nocontam.fa`

### Clustering ASVs into MOTUs

`Rscript novaseq/scripts/original_scripts/dereplication_headers.R`

`grep -v "^--" novaseq/18S/blank_corr/18S_nochim_nosingle_nocontam.fa | awk 'NR%2==0' | paste -d'\n' novaseq/18S/ASV_dereplicated.txt - > novaseq/18S/18S_dereplicated_ASVs.fa`

`cd novaseq/18S/`

`PATH/TO/SWARM -d 13 -i swarm/internal.txt -o swarm/output.txt -s swarm/statistics.txt -u swarm/uclust.txt -w swarm/18S_cluster_reps.fa 18S_dereplicated_ASVs.fa`

`popd`

### LULU curation

`Rscript novaseq/scripts/original_scripts/MOTU_tables.R`

`awk -F'_' '{print $1}' novaseq/18S/swarm/18S_cluster_reps.fa > novaseq/18S/MOTU/18S_cluster_reps_lulu_ready.fa`

`makeblastdb -in novaseq/18S/MOTU/18S_cluster_reps_lulu_ready.fa -parse_seqids -dbtype nucl`

`blastn -db novaseq/18S/MOTU/18S_cluster_reps_lulu_ready.fa -outfmt '6 qseqid sseqid pident' -out novaseq/18S/MOTU/match_list.txt -qcov_hsp_perc 80 -perc_identity 90 -query novaseq/18S/MOTU/18S_cluster_reps_lulu_ready.fa`

`Rscript novaseq/scripts/original_scripts/LULU_curation.R`

`grep -w -A 1 -f novaseq/18S/MOTU/lulu_curated_headers.txt novaseq/18S/MOTU/18S_cluster_reps_lulu_ready.fa --no-group-separator > novaseq/18S/MOTU/18S_cluster_reps_lulu_curated.fa`

`grep -w -A 1 -f novaseq/18S/MOTU/lulu_curated_headers.txt novaseq/18S/MOTU/18S_cluster_reps_lulu_ready.fa --no-group-separator > novaseq/18S/MOTU/18S_cluster_reps_lulu_curated.fa`

`Rscript novaseq/scripts/original_scripts/Subset_18S_taxonomy.R`