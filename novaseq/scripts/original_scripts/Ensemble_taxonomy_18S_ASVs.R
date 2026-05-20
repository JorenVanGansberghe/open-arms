#!/usr/bin/env Rscript


library(devtools)
devtools::install_github("dcat4/ensembleTax", build_manual = FALSE, build_vignettes = TRUE) # Installation issue when build_manual was set to TRUE (problem with pdflatex) 
library(ensembleTax)
library(tidyr)
library(dplyr)
library(Biostrings)
library(dada2)
library(stringr)



# specify input directory

input_dir <- file.path("novaseq", "18S")

# specify output directory 

output_dir <- file.path("novaseq", "18S")

# Read non-chimeric, non-singleton sequence table

seqtab.nochim.nosingle <- readRDS(file = file.path(output_dir, "seqtab_nochim_nosingle_18S.rds"))

# Read taxonomy objects (dada2 output)

silva <- readRDS(file = file.path(output_dir, "taxa_18S_silva.rds"))
silva.euk <- readRDS(file = file.path(output_dir,"taxa_18S_silva_euk.rds"))
pr2 <- readRDS(file = file.path(output_dir,"taxa_18S_pr2.rds"))

## Combine the two Silva taxonomies ##

# remove the Strain column of silva.euk

silva.euk$tax <- silva.euk$tax[,-8]
silva.euk$boot <- silva.euk$boot[,-8]

# Make dataframes of the tax objects of the two Silva taxonomies 

silva_tax <- as.data.frame(silva$tax)
silva_euk_tax <- as.data.frame(silva.euk$tax)

# Create two vectors for Class and Order_Family level
# For ASVs where Silva classification was not NA, this information was kept. 
# If Silva classification was NA, the classification of the Silva eukaryote reference set was kept.
# If both were NA, NA was set. 

class <- ifelse(!is.na(silva_tax$Class),silva_tax$Class,silva_euk_tax$Class)

order_family <- ifelse(!is.na(silva_tax$Order),silva_tax$Order,silva_euk_tax$Order_Family)

# Combine columns of Silva and Silva eukaryote taxonomies

silva_taxonomy <- cbind(silva_euk_tax[,1:4],silva_tax[,1:2],class,order_family,silva_euk_tax[,7])
colnames(silva_taxonomy)[7:9] <- c("Class","Order_Family","Species") # Set some column names

# One ASV was classified as Bacteria in the Silva classification. Taxonomy will be set to NA for all ranks

silva_taxonomy[which(silva_taxonomy$Kingdom=="Bacteria"),] <- NA

# remove Kingdom column

silva_taxonomy <- silva_taxonomy[,-5]

# Split Species column strings into separate columns
# First. make rownames as a column, as separate_wider_functions remove rownames (probably bug)

silva_taxonomy$seqs <- rownames(silva_taxonomy)
silva_taxonomy <- silva_taxonomy[,c(9,1:8)] # quick re-arranging of columns
silva_taxonomy <- silva_taxonomy %>% separate_wider_delim(Species, delim = "_", names_sep="",too_few = "align_start")

# Keep only first two columns of Species strings

silva_taxonomy <- silva_taxonomy[,-(11:ncol(silva_taxonomy))]

# Set second column of species strings to NA if it says sp. or cf. or environmental

silva_taxonomy[which(silva_taxonomy[,10]=="sp." | silva_taxonomy[,10]=="cf." | silva_taxonomy[,10]=="environmental"),10] <- NA

# Set species strings to NA if the first species string column contains "U/uncultured" or "unidentified"

silva_taxonomy[which(silva_taxonomy[,9]=="uncultured" | silva_taxonomy[,9]=="Uncultured" | silva_taxonomy[,9]=="unidentified"),9:10] <- NA

# Add a column with genus and species as one string

silva_taxonomy <- as.data.frame(silva_taxonomy)
silva_taxonomy$Species <- paste(silva_taxonomy[,9],silva_taxonomy[,10],sep="_")

# Set Species column to NA if the 10th column is NA

silva_taxonomy$Species <- ifelse(is.na(silva_taxonomy$Species2),NA,silva_taxonomy$Species)

# Make a genus level column

colnames(silva_taxonomy)[9] <- "Genus"

# Remove the remaining column of the second species string

silva_taxonomy <- silva_taxonomy[,-10]

# Save file

saveRDS(silva_taxonomy, file = file.path(output_dir, "silva_taxonomy.rds"))

# Generate a new Silva taxonomy object in the form of dada2's assignTaxonomy output

silva_taxonomy <- as.matrix(silva_taxonomy) # create character matrix
rownames(silva_taxonomy) <- silva_taxonomy[,1] # make sequences rownames
silva_taxonomy <- silva_taxonomy[,-1] # remove sequences column
# Make dummy bootstrap table
boots <- matrix(rep(100,nrow(silva_taxonomy)*9), nrow = nrow(silva_taxonomy), ncol = 9, byrow = TRUE,dimnames = list(rownames(silva_taxonomy),colnames(silva_taxonomy)))
silva_new <- silva
silva_new$tax <- silva_taxonomy
silva_new$boot <- boots

## ensembleTax workflow ##

# Make rubric with sequences mapped to short ASV IDs
rubric <- DNAStringSet(getSequences(seqtab.nochim.nosingle))
# this creates names (ASV1, ASV2, ..., ASVX) for each ASV
snam <- vector(mode = "character", length = length(rubric))
for (i in 1:length(rubric)) {
  snam[i] <- paste0("ASV", as.character(i))
}
names(rubric) <- snam

# Pre-process with the bayes2taxdf function
# ensembleTax supports Silva v138 and PR2 v4.14.0
# because we used PR2 v5.0.0 and customized a Silva taxonomy here, we will set db = NULL and provide the ranks present in our taxonomy assignments

silva.pretty <- bayestax2df(silva_new, 
                            db = NULL, 
                            ranks = colnames(silva_new$tax),
                            boot = 70,
                            rubric = rubric,
                            return.conf = FALSE)

pr2.pretty <- bayestax2df(pr2, 
                          db = NULL, 
                          ranks = colnames(pr2$tax),
                          boot = 70,
                          rubric = rubric,
                          return.conf = FALSE)

# Translate Silva taxonomic assignments onto the taxonomic nomenclature pf the PR2 assignments

# Because ensembleTax does not yet support PR2 v5.0.0., we create a custom tax2map2 object for use with taxmapper function
ff <- "novaseq/Taxonomy_reference_sets/pr2_version_5.1.1_SSU_dada2.fasta.gz"  # read PR2 reference fasta
fastaFile <- readDNAStringSet(ff)
seq_name = names(fastaFile)
taxmap <- str_split(seq_name, pattern = ";", simplify = TRUE)
taxmap <- as.data.frame(taxmap[, -ncol(taxmap)], stringsAsFactors = FALSE)
colnames(taxmap) <- colnames(pr2$tax)
taxmap <- unique(taxmap, MARGIN = 1)
any(is.na(taxmap)) # output should be FALSE

silva.mapped2pr2 <- taxmapper(silva.pretty,
                              tt.ranks = colnames(silva.pretty)[3:ncol(silva.pretty)],
                              tax2map2 = taxmap,
                              exceptions = c("Archaea", "Bacteria"),
                              ignore.format = TRUE,
                              synonym.file = "default",
                              streamline = TRUE,
                              outfilez = NULL)

## End ensembleTax workflow here. assign.ensembleTax function will not be used.
# assign.ensembleTax sets lower ranks to NA when a rank above is NA, even though two taxonomies may agree at this lower rank
# We therefore used alternative commands instead. See below.

# Continue...

# The info in the Order_Family column of the original Silva dada2 assignments are usually not present in the PR2 assignments
# Add this info as new column to the PR2 and mapped Silva assignments if the string in silva.pretty$Class matches the string in either pr2.pretty$Class_X or pr2.pretty$Class_Order_Family / silva.mapped2pr2$Class_X or silva.mapped2pr2$Class_Order_Family

order_family_x_1 <- ifelse(silva.pretty$Class==pr2.pretty$Class_X | silva.pretty$Class==pr2.pretty$Class_Order_Family,silva.pretty$Order_Family,NA)
pr2.pretty <- cbind(pr2.pretty[,1:9],order_family_x_1,pr2.pretty[,10:11])
colnames(pr2.pretty)[10] <- "Order_Family_X"
order_family_x_2 <- ifelse(silva.pretty$Class==silva.mapped2pr2$Class_X | silva.pretty$Class==silva.mapped2pr2$Class_Order_Family,silva.pretty$Order_Family,NA)
silva.mapped2pr2 <- cbind(silva.mapped2pr2[,1:9],order_family_x_2,silva.mapped2pr2[,10:11])
colnames(silva.mapped2pr2)[10] <- "Order_Family_X"

# Perform some formatting of both the Silva and PR2 assignments

is.na(silva.mapped2pr2) <- array(grepl('uncultured', as.matrix(silva.mapped2pr2)), dim(silva.mapped2pr2)) # Set as NA when "uncultured_" appears in Silva assignments
is.na(pr2.pretty) <- array(grepl('_X', as.matrix(pr2.pretty)), dim(pr2.pretty)) # Set as NA when "_X" appears in PR2 assignments
is.na(silva.mapped2pr2) <- array(grepl('_X', as.matrix(silva.mapped2pr2)), dim(silva.mapped2pr2)) # Set as NA when "_X" appears in Silva assignments
is.na(silva.mapped2pr2) <- array(grepl('_sp.', as.matrix(silva.mapped2pr2)), dim(silva.mapped2pr2)) # Set as NA when "_sp." appears in Silva species level assignments
is.na(pr2.pretty) <- array(grepl('_sp.', as.matrix(pr2.pretty)), dim(pr2.pretty))# Set as NA when "_sp." appears in PR2 species level assignments

# Set NAs to character string "NA" to enable the following procedures

silva.mapped2pr2[is.na(silva.mapped2pr2)] <- "NA"
pr2.pretty[is.na(pr2.pretty)] <- "NA"

## Generate final taxonomy
# For ASVs where PR22 and Silva agreed at the respective level, this information was kept. 
# For levels of an ASV were the two disagreed, the assignments were set to NA. 
# Where one of the two taxonomies was NA, but the other one was not, the latter's assignment was kept.

Domain <- ifelse(pr2.pretty$Domain==silva.mapped2pr2$Domain,pr2.pretty$Domain,ifelse(pr2.pretty$Domain=="NA",silva.mapped2pr2$Domain,ifelse(silva.mapped2pr2$Domain=="NA",pr2.pretty$Domain,"NA")))
Supergroup <- ifelse(pr2.pretty$Supergroup==silva.mapped2pr2$Supergroup,pr2.pretty$Supergroup,ifelse(pr2.pretty$Supergroup=="NA",silva.mapped2pr2$Supergroup,ifelse(silva.mapped2pr2$Supergroup=="NA",pr2.pretty$Supergroup,"NA")))
Division <- ifelse(pr2.pretty$Division==silva.mapped2pr2$Division,pr2.pretty$Division,ifelse(pr2.pretty$Division=="NA",silva.mapped2pr2$Division,ifelse(silva.mapped2pr2$Division=="NA",pr2.pretty$Division,"NA")))
Subdivision <- ifelse(pr2.pretty$Subdivision==silva.mapped2pr2$Subdivision,pr2.pretty$Subdivision,ifelse(pr2.pretty$Subdivision=="NA",silva.mapped2pr2$Subdivision,ifelse(silva.mapped2pr2$Subdivision=="NA",pr2.pretty$Subdivision,"NA")))
Phylum <- ifelse(pr2.pretty$Phylum==silva.mapped2pr2$Phylum,pr2.pretty$Phylum,ifelse(pr2.pretty$Phylum=="NA",silva.mapped2pr2$Phylum,ifelse(silva.mapped2pr2$Phylum=="NA",pr2.pretty$Phylum,"NA")))
Class_X <- ifelse(pr2.pretty$Class_X==silva.mapped2pr2$Class_X,pr2.pretty$Class_X,ifelse(pr2.pretty$Class_X=="NA",silva.mapped2pr2$Class_X,ifelse(silva.mapped2pr2$Class_X=="NA",pr2.pretty$Class_X,"NA")))
Class_Order_Family <- ifelse(pr2.pretty$Class_Order_Family==silva.mapped2pr2$Class_Order_Family,pr2.pretty$Class_Order_Family,ifelse(pr2.pretty$Class_Order_Family=="NA",silva.mapped2pr2$Class_Order_Family,ifelse(silva.mapped2pr2$Class_Order_Family=="NA",pr2.pretty$Class_Order_Family,"NA")))
Order_Family_X <- ifelse(pr2.pretty$Order_Family_X==silva.mapped2pr2$Order_Family_X,pr2.pretty$Order_Family_X,ifelse(pr2.pretty$Order_Family_X=="NA",silva.mapped2pr2$Order_Family_X,ifelse(silva.mapped2pr2$Order_Family_X=="NA",pr2.pretty$Order_Family_X,"NA")))
Genus <- ifelse(pr2.pretty$Genus==silva.mapped2pr2$Genus,pr2.pretty$Genus,ifelse(pr2.pretty$Genus=="NA",silva.mapped2pr2$Genus,ifelse(silva.mapped2pr2$Genus=="NA",pr2.pretty$Genus,"NA")))
Species <- ifelse(pr2.pretty$Species==silva.mapped2pr2$Species,pr2.pretty$Species,ifelse(pr2.pretty$Species=="NA",silva.mapped2pr2$Species,ifelse(silva.mapped2pr2$Species=="NA",pr2.pretty$Species,"NA")))

tax_table <- cbind(pr2.pretty[,1:2],Domain,Supergroup,Division,Subdivision,Phylum,Class_X,Class_Order_Family,Order_Family_X,Genus,Species)


## Write a fasta file of the non-chimeric and non-singleton sequences with >ASV... type headers

asv_seqs <- tax_table$ASV
asv_headers <- paste0(">",tax_table$svN)
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file = file.path(output_dir, "18S_nochim_nosingle_ASVs.fa"))

# Write an ASV count table of the non-chimeric and non-singleton sequences with short >ASV... type names

ASV_counts <- t(seqtab.nochim.nosingle) # transposing table

# Sort count table based on sequence order in tax_table and write file with ASV names

ASV_counts <- ASV_counts[order(match(rownames(ASV_counts),tax_table[,2])),]

rownames(ASV_counts) <- tax_table$svN

write.table(ASV_counts,file = file.path(output_dir, "18S_ASV_counts_nosingle.txt"), sep="\t", quote=F, col.names=NA)

# write final taxonomy table to file

tax_table <- tax_table[,-2]
colnames(tax_table)[1] <- "ASV"
write.table(tax_table, file = file.path(output_dir, "18S_tax_table.txt"), sep="\t",row.names = F)