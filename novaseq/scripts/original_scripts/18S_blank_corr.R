#!/usr/bin/env Rscript

### 18S ###

output_dir <- "novaseq/18S/blank_corr"

if(!dir.exists(output_dir)) dir.create(output_dir)

# read 18S sample summary and extract negative control IDs

18S_summary <- read.csv("metadata/generated_meta/18S_batch3.4.5.6-metadata.csv", 
                        sep = ",", header = TRUE)
18S_negative <- c(unique(18S_summary$PCR_negative_control_Code_1), unique(18S_summary$PCR_negative_control_Code_2))

# Read ASV count table (output from dada2)

asv_table <- read.table("novaseq/18S/18S_ASV_counts_nosingle.txt",
                        sep = "\t", header = T, row.names = 1,
                        as.is = T, check.names = F)

# Subset rows where the ASV read count in the blank / negative samples  exceeds 10 % of an ASVs total read count 
# before combining the tables using rbind, add a column with ASV IDs. 
# This is necessary to stop R from introducing new rownames adding zeros to ASV names if duplicate ASVs exist in the newly created blank dataframes

blank_list <- list()

for (id in 18S_negative) {
  subset_id <- asv_table[grep(id, colnames(asv_table))]

  negative_id <- subset(asv_table, subset_id > 0.1 * rowSums(asv_table))
  blank <- cbind(ASV_ID = rownames(negative_id), negative_id)

  blank_list[[id]] <- blank
}

# Combine blank dataframes
blanks <- do.call(rbind, unname(blank_list))

# Remove ASVs if there are duplicates in the "blanks" table (in case an ASV's read count exceeded 10 % of the total read count in more than one blank sample)

blanks <- blanks[!duplicated(blanks$ASV_ID), ]

# Remove the potential contaminant ASVs from the ASV table

asv_table <- cbind(ASV_ID = rownames(asv_table), asv_table)
asv_no_contams <- asv_table[!(asv_table$ASV_ID %in% blanks$ASV_ID),]

# Write table of potential contaminant ASVs and ASV count table devoid of contaminants

write.table(blanks, file = file.path(output_dir, "asv_contaminants_18S.txt"), sep="\t", row.names = F)

write.table(asv_no_contams, file = file.path(output_dir, "asv_no_contaminants_18S.txt"), sep="\t", row.names = F)

# Write headers of non-contaminant asvs to file to subset the corresponding sequences of non-NUMT fasta file later on

no_contam_headers <- paste0(">",asv_no_contams$ASV_ID)
write.table(no_contam_headers, file = file.path(output_dir, "no_contam_headers_18S.txt"), sep="\t", row.names = F,quote=F,col.names = F)
