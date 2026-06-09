#!/usr/bin/env Rscript

library(dada2)

# specify input directory 
input_dir <- file.path("novaseq", "18S")

# specify output directory 
output_dir <- file.path("novaseq", "18S")

# load sequencing batch numbers
batch_dir    <- file.path("novaseq", "18S", "fastq_files")
#batch_list <- list.files(batch_dir, pattern = "Batch")
batch_list <- grep("^(?!.*Batch_5).*Batch", 
                   list.files(batch_dir), 
                   value = TRUE, perl = TRUE)

get_sample_name <- function(fname) strsplit(basename(fname), "_")[[1]][2]
run_18S <- unname(sapply(batch_list, get_sample_name))

# Load the sequence tables of the different sequence runs

rds_list <- list()
for (run in run_18S) {
  run_file <- file.path(input_dir, paste("seqtab_Batch", run, "_mod4.rds", sep = ""))
  rds_run <- readRDS(run_file)
  rds_list[[as.character(run)]] <- rds_run
}

# Merge sequence tables

if (length(rds_list) > 1) {
  # Merge sequence tables
  merged <- do.call(mergeSequenceTables, unname(rds_list))
} else { 
  merged <- as.matrix(rds_list[[1]])
}

saveRDS(merged, file = file.path(output_dir, "merged_seqtab_18S.rds"))


# Remove chimeras #

seqtab.nochim <- removeBimeraDenovo(merged, multithread=T, verbose=TRUE)

# Save sequence table with the non-chimeric sequences as RDS file:

saveRDS(seqtab.nochim, file = file.path(output_dir, "seqtab_nochim_18S.rds"))

# It is possible that a large fraction of the total number of UNIQUE SEQUENCES will be chimeras.
# However, this is usually not the case for the majority of the READS.
# Calculate percentage of the reads that were non-chimeric.

non_chimeric <- (sum(seqtab.nochim)/sum(merged))*100
message(non_chimeric, "% of sequences were non-chimeric")

# Remove singletons from the non-chimeric ASVs

#Transform counts to numeric (as they will most likely be integers)
mode(seqtab.nochim) = "numeric"

# Subset columns with counts of > 1 and save to file
seqtab.nochim.nosingle <- seqtab.nochim[,colSums(seqtab.nochim) > 1]
saveRDS(seqtab.nochim.nosingle, file = file.path(output_dir, "seqtab_nochim_nosingle_18S.rds"))


## Track reads through the entire dada2 pipeline ##

# Track reads through the chimera and singleton removal step.

# Read non-chimeric and non-singleton table again, as it has been modified
seqtab.nochim.nosingle <- readRDS(file = file.path(output_dir, "seqtab_nochim_nosingle_18S.rds"))

track_nochim_nosingle <- cbind(rowSums(seqtab.filtered), rowSums(seqtab.nochim), rowSums(seqtab.nochim.nosingle))
colnames(track_nochim_nosingle) <- c("nonchim", "nosingle")

# Read tracking tables of the single runs and combine

track_list <- list()
for (run in run_18S) {
  track_file <- file.path(input_dir, paste("track_Batch", run, "_mod4.txt", sep = ""))
  txt_track <- read.table(track_file, sep = "\t", header = T, row.names = 1)
  track_list[[as.character(run)]] <- txt_track
}

tracks <- do.call(rbind, unname(track_list))

# Combine all tracking tables

tracks <- tracks[order(match(rownames(tracks), rownames(track_nochim_nosingle))),]
track_18S <- cbind(tracks, track_nochim_nosingle)

# Calculate percentages for each step compared to input

track_18S$cutadapt_perc <- (track_18S$cutadapt / track_18S$input)*100
track_18S$filtered_perc <- (track_18S$filtered / track_18S$input)*100
track_18S$denoisedF_perc <- (track_18S$denoisedF / track_18S$input)*100
track_18S$denoisedR_perc <- (track_18S$denoisedR / track_18S$input)*100
track_18S$merged_perc <- (track_18S$merged / track_18S$input)*100
track_18S$nonchim_perc <- (track_18S$nonchim / track_18S$input)*100
track_18S$nosingle_perc <- (track_18S$nosingle / track_18S$input)*100

# Save final read tracking table to file

write.table(track_18S, file = file.path(output_dir, "track_18S.txt"), sep = "\t", col.names = NA)



## Assign taxonomy ##

# Read non-chimeric, no-singleton sequence table again (in case the processes above have altered it)

seqtab.nochim.nosingle<-readRDS(file = file.path(output_dir, "seqtab_nochim_nosingle_18S.rds"))

## Official Silva v138.2 (prokaryote and eukaryote reference set)

# Download here: https://zenodo.org/records/14169026

silva.ref<-"novaseq/Taxonomy_reference_sets/silva_nr99_v138.2_toSpecies_trainset.fa.gz"

# Set minBoot to 70

taxa_silva <- assignTaxonomy(seqtab.nochim.nosingle, silva.ref, minBoot=70, multithread=T,outputBootstraps = T)

saveRDS(taxa_silva, file = file.path(output_dir, "taxa_18S_silva.rds"))

## Contributed Silva Eukaryote v132 99% clustered reference set

# Download here: https://zenodo.org/record/1447330

silva.euk.ref <- "novaseq/Taxonomy_reference_sets/silva_132.18s.99_rep_set.dada2.fa.gz" 

# Set minBoot to 70.
# Define taxonomic ranks for this specific reference set

taxa_silva_euk <- assignTaxonomy(seqtab.nochim.nosingle, silva.euk.ref, minBoot=70, multithread=T,outputBootstraps = T,taxLevels=c("Domain","Division","Division_X","Subdivision","Class","Order_Family","Species","Strain"))

saveRDS(taxa_silva_euk, file = file.path(output_dir,"taxa_18S_silva_euk.rds"))

## PR2 v5.0.0

# Download here: https://github.com/pr2database/pr2database/releases

pr2.ref <- "novaseq/Taxonomy_reference_sets/pr2_version_5.1.1_SSU_dada2.fasta.gz" 

# Set minBoot to 70.
# Define taxonomic ranks for this specific reference set

taxa_pr2 <- assignTaxonomy(seqtab.nochim.nosingle, pr2.ref, minBoot=70, multithread=T,outputBootstraps = T, taxLevels = c("Domain","Supergroup","Division","Subdivision","Phylum","Class_X","Class_Order_Family","Genus","Species"))

saveRDS(taxa_pr2, file = file.path(output_dir,"taxa_18S_pr2.rds"))

### Continue with next script for 18S ensemble taxonomy ###

# ASV count and taxonomy tables, as well as the fasta file with ASV sequences will be produced there #

