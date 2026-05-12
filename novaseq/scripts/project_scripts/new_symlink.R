#!/usr/bin/env Rscript

# =============================================================================
# SYMLINK FASTQ FILES FOR COI/18S METABARCODING PIPELINE
# =============================================================================
# This script creates symbolic links from the Genoscope delivery directory
# to a structured batch directory for downstream processing.
#
# USAGE:
#   - First run:  Run the full script as-is
#   - Adding new batches: Update the CSV path below (e.g., COI_batch3.4.5.6.csv)
#                         and re-run the full script. It will skip existing 
#                         symlinks and preserve batch numbering automatically.
# =============================================================================

#library(dplyr)

# ---- CONFIGURATION (edit these as needed) ------------------------------------

# Source directory (Genoscope delivery)
source_dir <- "/cfs/klemming/projects/supr/naiss2025-23-46/ARMS_Genoscope_data_full/projet_DBB/ARMS-Macro/" # directory where eDNA fastq data is located
pcr1_dir   <- file.path(source_dir, "NEGATIVE_CONTROLS", "pcr1")
pcr2_dir   <- file.path(source_dir, "NEGATIVE_CONTROLS", "pcr2")
code_dir   <- file.path(source_dir, "18S_V1V2", "18S_V1V2_Metazoaire_SSUF04_SSURmod") # directory with the primers

# Target directory (your pipeline input)
fastq_dir  <- "/cfs/klemming/projects/supr/naiss2025-23-46/Joren/Project_18S/novaseq/18S/fastq_files"

# Metadata CSV — UPDATE THIS WHEN ADDING NEW BATCHES
csv_path   <- "metadata/generated_meta/18S_batch3.4.5.6-metadata.csv"

# Batch translations file (written to fastq_dir)
batch_file <- file.path(fastq_dir, "Batch_name_translations_18S.csv")

# ---- HELPER FUNCTIONS --------------------------------------------------------

# Safe symlink: skips if target already exists
safe_symlink <- function(source, target_path){
    if (!file.exists(target_path)){
        file.symlink(source, target_path)
    }
}

# Find the correct batch directory for a file based on flowcell ID
symlink_target_function <- function(file){
    target_dir <- NULL
    for (i in 1:nrow(sample_batch)){
        if (grepl(sample_batch[i, 2], file) == TRUE){
            target_dir <- file.path(fastq_dir, paste0("Batch_", sample_batch[i, 1]))
            if(!dir.exists(target_dir)) dir.create(target_dir)
        }
    }
    if (is.null(target_dir)) {
        message("WARNING: No batch match for file: ", file)
    }
    return(target_dir)
}

# Ensure target directory exists
dir.create(fastq_dir, recursive = TRUE, showWarnings = FALSE)

# ---- READ AND FILTER METADATA ------------------------------------------------

samples <- read.csv(csv_path, sep = ",", header = TRUE, fileEncoding = "latin1")

# Keep only rows with valid codes
samples <- samples[nchar(trimws(samples$Code)) > 0 & grepl("^DBQ_", samples$Code), ]

# List source files (filter to fastq only)
pcr1_files <- list.files(pcr1_dir)
pcr1_files <- pcr1_files[grepl("\\.fastq", pcr1_files)]
pcr2_files <- list.files(pcr2_dir)
pcr2_files <- pcr2_files[grepl("\\.fastq", pcr2_files)]
code_files <- list.files(code_dir)

# ---- PARSE SAMPLE CODES -----------------------------------------------------

# Extract PCR negative control codes
pcr_control1 <- unique(samples$PCR_negative_control_Code_1)
pcr_control2 <- unique(samples$PCR_negative_control_Code_2)

# Split Code column: DBQ_AAABOSTA_1_HKYCCDRX3.UDI366-BID16
#   -> col 1: DBQ, col 2: AAABOSTA, col 3: 1, col 4: HKYCCDRX3.UDI366-BID16
#   -> col 5: HKYCCDRX3 (flowcell ID), col 6: UDI366-BID16
code_split <- do.call(rbind, strsplit(samples$Code, "_"))
code_split <- cbind(code_split, do.call(rbind, strsplit(code_split[, 4], "\\.")))

# Sample names (first two parts of Code)
sample_names <- unique(paste(code_split[,1], code_split[,2], sep = "_"))

# ---- ASSIGN BATCH NUMBERS (preserves existing numbering) --------------------

all_flowcells <- unique(code_split[, 5])

if (file.exists(batch_file)){
    # Incremental run: keep existing batch numbers, add new flowcells
    existing_batches <- read.csv(batch_file)
    existing_batches$Batch_X <- as.integer(trimws(existing_batches$Batch_X))
existing_batches$seq_run_code <- trimws(existing_batches$seq_run_code)
    new_flowcells <- all_flowcells[!all_flowcells %in% existing_batches$seq_run_code]
    if (length(new_flowcells) > 0){
        start <- max(existing_batches$Batch_X) + 1
        new_rows <- data.frame(
            Batch_X = start:(start + length(new_flowcells) - 1),
            seq_run_code = new_flowcells
        )
        sample_batch <- rbind(existing_batches, new_rows)
        message("Added ", length(new_flowcells), " new flowcell(s): ", 
                paste(new_flowcells, collapse = ", "))
    } else {
        sample_batch <- existing_batches
        message("No new flowcells found")
    }
} else {
    # First run: assign batch numbers from scratch
    sample_batch <- data.frame(
        Batch_X      = seq_len(length(all_flowcells)),
        seq_run_code = all_flowcells
    )
    colnames(sample_batch) <- c("Batch_X", "seq_run_code")
    message("First run: assigned ", nrow(sample_batch), " batch(es)")
}

message("Batch assignments:")
print(sample_batch)

# ---- SYMLINK PCR1 NEGATIVE CONTROLS -----------------------------------------

for (ctrl in pcr_control1){
    pcr1_match <- grep(ctrl, pcr1_files, value = TRUE)
    if (length(pcr1_match) >= 1){
        for (f in pcr1_match){
            target <- symlink_target_function(f)
            if (!is.null(target)){
                safe_symlink(file.path(pcr1_dir, f), file.path(target, f))
            }
        }
    }
}
message("PCR1 controls symlinked")

# ---- SYMLINK PCR2 NEGATIVE CONTROLS -----------------------------------------

for (ctrl in pcr_control2){
    pcr2_match <- grep(ctrl, pcr2_files, value = TRUE)
    if (length(pcr2_match) >= 1){
        for (f in pcr2_match){
            target <- symlink_target_function(f)
            if (!is.null(target)){
                safe_symlink(file.path(pcr2_dir, f), file.path(target, f))
            }
        }
    }
}
message("PCR2 controls symlinked")

# ---- SYMLINK SAMPLE FILES ---------------------------------------------------

# Find subdirectories matching our flowcell IDs
batch_list <- list()
for (i in 1:nrow(sample_batch)){
    match <- grep(sample_batch[i, 2], code_files, value = TRUE)
    if (length(match) > 0) {
        batch_list <- c(batch_list, match)
    }
}

# Symlink sample fastq files
for (seq_run in batch_list){
    subdir_path <- file.path(code_dir, seq_run)
    subdir_files <- list.files(subdir_path)
    for (name in sample_names){
        code_match <- grep(name, subdir_files, value = TRUE)
        if (length(code_match) >= 1){
            for (f in code_match){
                target <- symlink_target_function(f)
                if (!is.null(target)){
                    safe_symlink(file.path(subdir_path, f), file.path(target, f))
                }
            }
        }
    }
}
message("Sample files symlinked")

# ---- EXPORT BATCH TRANSLATIONS -----------------------------------------------

write.csv(sample_batch, file = batch_file, row.names = FALSE, quote = FALSE)

message("Script finished successfully")