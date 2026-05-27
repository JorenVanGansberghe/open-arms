library(dplyr)

# specify input directory

input_dir <- file.path("novaseq", "18S")

# specify output directory 

output_dir <- file.path("novaseq", "18S")

# specify 18S/MOTU directory

MOTU_dir <- file.path("novaseq", "18S", "MOTU")

# Read the ensemble taxonomy table generated after the dada2 pipeline

tax_table <- read.table(file = file.path(input_dir, "18S_tax_table.txt"), sep="\t", header=T)

# Read the headers of LULU curated MOTUs

headers<-read.table(file = file.path(MOTU_dir, "lulu_curated_headers.txt"), sep="\t")
headers[,1]<-gsub(">","",headers[,1]) # Remove the ">" from the MOTU names

# Subset the taxonomy assignments of ASVs which are the representative ASV sequences of LULU curated MOTUs

taxa_lulu<-tax_table %>% filter(ASV %in% headers[,1])
colnames(taxa_lulu)[1]<-"MOTU"

# Write to file

write.table(taxa_lulu, file = file.path(output_dir, "18S_final_tax_table.txt"), sep="\t", row.names = F)
