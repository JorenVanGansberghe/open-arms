------------------------------------------------------------------------

editor_options: markdown: wrap: 72 ---

# ARMS-MBON 18S Processing Pipeline

This repository contains the scripts and workflow used to process ARMS-MBON 18S metabarcoding data from 2022-2025. The processed results are used to generate summary tables, visualizations, taxonomic assignments and supporting material for a data paper.

The workflow is designed to run on the Dardel HPC system using SLURM and an Apptainer sandbox. The main pipeline combines R scripts and command-line bioinformatics tools, and uses an ensemble approach to taxonomic assignment combining three reference databases.

## 1. Project overview

The aim of this project is to process a large ARMS-MBON dataset using an existing R-script pipeline. The pipeline starts from paired-end 18S V1V2 FASTQ files and produces curated MOTU tables and taxonomy tables that can be used for summaries, graphs and a final data paper.

The workflow does the following:

Raw paired-end 18S FASTQ files

-\> Primer removal with cutadapt

-\> Quality filtering and ASV inference with DADA2

-\> Chimera and singleton removal

-\> Taxonomic classification using three reference databases (Silva v138.2, Silva Euk v132, PR2 v5.1.1)

-\> Ensemble taxonomy combining the three classifier outputs

-\> Negative control / blank correction

-\> ASV dereplication and abundance-labelled FASTA creation

-\> Clustering ASVs into MOTUs using swarm

-\> MOTU table creation

-\> LULU curation

-\> Subsetting final taxonomy to curated MOTUs

-\> Final count and taxonomy tables

## 2. Software and tools

The pipeline uses the following tools:

HPC / container system:

- Dardel HPC
- SLURM
- Apptainer
- PDC module environment

R packages:

- dada2
- ShortRead
- Biostrings
- ggplot2
- ensembleTax
- lulu
- tidyr
- dplyr
- stringr
- devtools
- argparse

Command-line tools: - cutadapt - swarm - BLAST+

## 3. Making SSH keys to log into Dardel server

#### 3.1 Make an account on NAISS SUPR:

- <https://supr.naiss.se/login/>
- Register New Person
- Register without Federal Identity
- Fill in form:
  - Personal email is okay, GU email also possible
  - Check “Swedish Academic Employee?”
  - Postal Code and City should be the one from your residence in Sweden
- Need to accept NAISS User Agreement
  - Insert GU mail (only necessary if personal email was used to make account)
- Ask Matthias to include you in the projects:
  - NAISS 2025/22-761 (PDC)
  - NAISS 2025/23-703 (PDC)

#### 3.2 Request your PDC account:

- Go to “accounts” tab in sidebar on SUPR
- Request your PDC account
- Wait for Centre to create your account
  - If it takes longer than 3 days for them to activate your account, you can ask them where in the activation process they are via following link: <https://supr.naiss.se/support/?centre_resource=c7>
- A mail containing your username will be send when your account is activated

#### 3.3 Generate SSH keys:

It is possible to work in linux, windows or both (macOS is also possible). I worked mostly in Visual Studio Code on Windows

Linux:

- Open terminal

```{bash}
cd
ls -a 
# If .ssh is already present, skip next line
mkdir .ssh
ssh-keygen -t ed25519 -f ~/.ssh/id-ed25519-pdc
```

- push enter twice to skip setting a passphrase

Windows:

- Open "Command prompt"

```{bash}
ssh-keygen -t ed25519 -f %USERPROFILE%.ssh\id-ed25519-pdc
```

- push enter twice to skip setting a passphrase

#### 3.4 Add new keys to your PDC account

- Go to PDC login portal: <https://loginportal.pdc.kth.se/>
- Log into SUPR and click "Prove My Identity to PDC"
- You are now redirected to login portal
- Click on "Add new key"
- Upload public key file or paste public key in "SSH public key"-field
  - Public key file is called "id-ed25519-pdc.pub"
    - In windows, the file name extensions have to be shown to see which file is the public key:
    - "view" -\> "show" -\> 'File name extension"
- Choose key name
- Save

Add a public key for each operating system you will work in.

You also have to "Add address" each time you work from a different WiFi or if the WiFi cycles through IP-addresses. The WiFi at Gothenburg University can cycle through multiple IP-addresses per day.

## 4. Log into dardel

Linux:

```{bash}
ssh -i ~/.ssh/id-ed25519-pdc your_pdc_username@dardel.pdc.kth.se
```

Windows:

In windows, it is easiest to use Visual Studio Code to work on the dardel cluster.

Visual Studio Code:

Setting up connection:

- Install the extension "Remote - SSH" by Windows

- Click "\>\<" at left bottom of VSC window

- Select "Connect to Host..."

- Select "Configure SSH Host..."

- Select "C:\\Users\\username.ssh\\config"

- Insert following text in config file:

  ::: {}
  | Host dardel.pdc.kth.se
  |   HostName dardel.pdc.kth.se
  |   IdentityFile \~/.ssh/id-ed25519-pdc
  |   User dardel_username
  :::

Logging in:

- Click "\>\<" at left bottom of VSC window
- Select "Connect to Host..."
- Click "dardel.pdc.kth.se"

## 5. Clone Github directory on Dardel

To clone the forked repository containing the scripts, first a connection with GitHub has to be set up.

```{bash}
cd
ssh-keygen -t ed25519 -C "githubaccount@student.howest.be"
eval $(ssh-agent -s)
ssh-add ~/.ssh/id_github_ed25519
ssh -T git@github.com

```

After setting up the connection, your public key has to be added to your GitHub.

```{bash}
cat ~/.ssh/id_github_ed25519.pub
```

Copy the output and add it in GitHub -\> Settings -\> SSH and GPG keys -\> New SSH key

It is best to make a new directory in the project directory, because your working directory does not have sufficient storage

```{bash}
cd /cfs/klemming/projects/supr/naiss2025-23-46
mkdir new_directory_name
```

Now clone the repository in the new directory

```{bash}
cd new_directory_name
git clone git@github.com:GitHub_account_name/repository_name.git
```

## 6. Setting up the sandbox on dardel

To be able to run your pipeline, a sandbox has to be set up. Apptainer is used, which is similar to a docker container but is optimized for high-performance computing (HPC). The sandbox only has to be set up once.

```{bash}
# makes it possible to use PDC commands, needed to load apptainer
module load PDC/24.11

# load apptainer
module load apptainer/1.4.0-cpeGNU-24.11

# build container in temp directory
  # $PDC_TMP is your home directory
apptainer build --sandbox $PDC_TMP/temparms-sandbox docker://continuumio/miniconda3:25.3.1-1

# make directory for the conda environment
mkdir $PDC_TMP/conda_envs

# creates conda environment from the repository environment
apptainer exec --bind /cfs/klemming/projects/supr/naiss2025-23-46/username/Project_18S:/work --bind $PDC_TMP/conda_envs:/envs  $PDC_TMP/temparms-sandbox conda env create  -f /work/environment/git_env.yml -p /envs/git_env
```

## 7. Other important bash commands on dardel

```{bash}
# information about the projects on the server
projinfo

# information about the different partitions on the server
sinfo

# run script on dardel
sbatch path/to/novaseq_slurm.sh

# show jobs from user that are queueing to be run on Dardel and how long they have been running
squeue -u pdcUsername

# show start time of jobs of user
squeue -u pdcUsername --start

# show start time specific job queueing
  # jobid is shown in first column of previous 2 commands
squeue --start -j jobid

# cancel a job
scancel jobid
```

## 8. Unable to log into dardel with VSCode after a while?

When logging into dardel using visual studio code, it installs a server component in the .vscode-server directory. Every time visual studio code is updated, a new server component is installed on dardel without deleting the old one. This results in .vscode-server growing in size over time until your homefolder exceeds its allowed size, resulting in the inability to log into dardel using vscode.

If this occurs, the best way to solve this is to log into dardel on your virtual machine and manually remove the .vscode-server directory, as it will be reinstalled upon login into dardel with vscode.

```{bash}
# remove .vscode-server directory on dardel
rm -rf ~/.vscode-server
```

The easiest way to prevent this problem is to insert the following bash script in a .bashrc file in your home directory. It deletes the older server components, preventing the .vscode-server to balloon in size.

```{bash}
# Auto-cleanup old VS Code server versions (keeps most recent)
cleanup_vscode_server() {
    local servers_dir="$HOME/.vscode-server/cli/servers"
    # Executes code if the directory at server_dir exists
    if [ -d "$servers_dir" ]; then
        # keep=2 will keep lru.json and the newest server installation 
        local keep=2
        # Earlier server installations are deleted
        ls -t "$servers_dir" | tail -n +$((keep + 1)) | while read -r old; do
            echo "Removing old VS Code server: $old"
            rm -rf "${servers_dir:?}/$old"
        done
    fi
}
cleanup_vscode_server
```

## 9. Input data organisation

The first R script expects FASTQ files to be organised by batch under:

`novaseq/18S/fastq_files/`

The data files can be ordered by different batches inside the /fastq_files directory

For example: `novaseq/18S/fastq_files/Batch_1`

To get the FASTQ files in the correct place, use the Rscript `new_symlink.R`

### 9.2 Reference databases

Download the following reference databases and place them in `novaseq/Taxonomy_reference_sets/`:

| Database | Version | Download |
|----|----|----|
| Silva (prokaryote + eukaryote) | v138.2 | [Zenodo](https://zenodo.org/records/14169026) |
| Silva Eukaryote | v132, 99% clustered | [Zenodo](https://zenodo.org/record/1447330) |
| PR2 | v5.1.1 | [GitHub releases](https://github.com/pr2database/pr2database/releases) |

Expected filenames:

```         
metadata/Taxonomy_reference_sets/silva_nr99_v138.2_toSpecies_trainset.fa.gz
metadata/Taxonomy_reference_sets/silva_132.18s.99_rep_set.dada2.fa.gz
metadata/Taxonomy_reference_sets/pr2_version_5.1.1_SSU_dada2.fasta.gz
```

## 10. Metadata

The metadata file used for the current 2022-2025 ARMS-MBON processing is:

`metadata/generated_meta/18S_batch3.4.5.6-metadata.csv`

This metadata is used later in the pipeline, especially for blank correction.

The metadata also contains sample identifiers that can be used for downstream plotting, for example to extract deployment and retrieval dates from MaterialSampleID.

## 11. Pipeline steps

### Step 1: Symlink FASTQ files

Script: `new_symlink.R`

Run command: `Rscript new_symlink.R`

This script organises raw FASTQ files from the Genoscope delivery directory into a structured batch directory using symbolic links. Batch numbers are assigned based on flowcell IDs and preserved across incremental runs. Edit `source_dir`, `fastq_dir` and `csv_path` at the top of the script before running.

Main outputs include:

`novaseq/18S/fastq_files/Batch_X/` (symlinked FASTQ files per batch)

`novaseq/18S/fastq_files/Batch_name_translations_18S.csv`

### Step 2: Primer removal, filtering, trimming and ASV inference

Script: `loessErrfun_mod4_sol.R`

Run command: `Rscript loessErrfun_mod4_sol.R -d /envs/git_env/bin/cutadapt`

This script performs:

- batch detection
- paired FASTQ detection - primer checking
- primer removal with cutadapt
  - (FWD: `GCTTGTCTCAAAGATTAAGCC`, REV: `CCTGCTGCCTTCCTTRGA`)
- quality filtering and trimming with DADA2
  - (maxEE = 2/4, minLen = 50)
- error modelling using a modified loess error function suited for NovaSeq binned quality scores
- ASV inference with pseudo-pooling
- read merging
  - (minOverlap = 10, maxMismatch = 1)
- ASV sequence table creation

Specify which batches to process by editing the `batch_list` variable at the bottom of the script.

Main outputs include:

`novaseq/18S/seqtab_Batch..._mod4.rds`

`novaseq/18S/track_Batch..._mod4.txt`

`novaseq/18S/18S_..._quality_forward.jpg`

`novaseq/18S/18S_..._quality_reverse.jpg`

### Step 3: Chimera and singleton removal and taxonomic classification

Script: `18S_chimera_singleton_removal_taxonomic_classification.R`

Run command: `Rscript 18S_chimera_singleton_removal_taxonomic_classification.R`

This script merges batch sequence tables, removes chimeras using `removeBimeraDenovo`, removes singleton ASVs (total count ≤ 1), and assigns taxonomy using three reference databases at minBoot = 70.

Main outputs include:

`novaseq/18S/seqtab_nochim_nosingle_18S.rds`

`novaseq/18S/taxa_18S_silva.rds`

`novaseq/18S/taxa_18S_silva_euk.rds`

`novaseq/18S/taxa_18S_pr2.rds`

`novaseq/18S/track_18S.txt`

### Step 4: Ensemble taxonomy

Script: `Ensemble_taxonomy_18S_ASVs.R`

Run command: `Rscript Ensemble_taxonomy_18S_ASVs.R`

This step builds a consensus taxonomy by merging the two Silva reference databases into a single combined Silva assignment, translating this onto the PR2 taxonomic framework using `ensembleTax::taxmapper`, and generating final rank-by-rank assignments. Where PR2 and Silva agree, that assignment is kept. Where one is NA the other is used. Where they conflict, NA is set.

Main outputs include:

`novaseq/18S/18S_nochim_nosingle_ASVs.fa`

`novaseq/18S/18S_ASV_counts_nosingle.txt`

`novaseq/18S/18S_tax_table.txt`

### Step 5: Negative control / blank correction

Script: `18S_blank_corr.R`

Run command: `Rscript 18S_blank_corr.R`

This step uses metadata to identify PCR negative controls and removes ASVs that are likely contaminants. An ASV is flagged if its read count in any negative control exceeds 10% of its total read count across all samples.

Main outputs include:

`novaseq/18S/blank_corr/asv_contaminants_18S.txt`

`novaseq/18S/blank_corr/asv_no_contaminants_18S.txt`

`novaseq/18S/blank_corr/no_contam_headers_18S.txt`

### Step 6: Dereplication headers

Script: `dereplication_headers.R`

Run command: `Rscript dereplication_headers.R`

This step creates FASTA headers containing ASV abundance information, which is needed for swarm clustering. The shell script then combines these headers with the blank-corrected FASTA sequences:

`grep -v "^--" 18S_nochim_nosingle_nocontam.fa | awk 'NR%2==0' | paste -d'\n' ASV_dereplicated.txt - > 18S_dereplicated_ASVs.fa`

Main outputs include:

`novaseq/18S/ASV_dereplicated.txt`

`novaseq/18S/18S_dereplicated_ASVs.fa`

### Step 7: MOTU clustering with swarm

Command: `/envs/git_env/bin/swarm \   -d 13 \   -i swarm/internal.txt \   -o swarm/output.txt \   -s swarm/statistics.txt \   -u swarm/uclust.txt \   -w swarm/18S_cluster_reps.fa \   18S_dereplicated_ASVs.fa`

This clusters ASVs into MOTUs using swarm with distance parameter -d 13.

Main outputs include:

`novaseq/18S/swarm/output.txt`

`novaseq/18S/swarm/uclust.txt`

`novaseq/18S/swarm/18S_cluster_reps.fa`

### Step 8: MOTU table creation

Script: `MOTU_tables.R`

Run command: `Rscript MOTU_tables.R`

This script converts ASV-level count information and swarm clustering output into a MOTU table by mapping ASVs to their representative MOTU and aggregating read counts.

Main outputs include:

`novaseq/18S/MOTU/motu_table_18S.txt`

### Step 9: BLAST match list for LULU

Command: `makeblastdb -in novaseq/18S/MOTU/18S_cluster_reps_lulu_ready.fa -parse_seqids -dbtype nucl`

`blastn -db novaseq/18S/MOTU/18S_cluster_reps_lulu_ready.fa -outfmt '6 qseqid sseqid pident' -out novaseq/18S/MOTU/match_list.txt -qcov_hsp_perc 80 -perc_identity 90 -query novaseq/18S/MOTU/18S_cluster_reps_lulu_ready.fa`

This creates a sequence similarity match list required by LULU.

Main outputs include:

`novaseq/18S/MOTU/match_list.txt`

### Step 10: LULU curation

Script: `LULU_curation.R`

Run command: `Rscript LULU_curation.R`

LULU reduces likely erroneous MOTUs using sequence similarity and co-occurrence information. Parameters used: `minimum_match = 0.90`, `minimum_ratio = 100`, `minimum_relative_cooccurence = 0.95`.

Main outputs include:

`novaseq/18S/MOTU/lulu_motu_table_18S.txt`

`novaseq/18S/MOTU/motu_map_lulu_18S.txt`

`novaseq/18S/MOTU/lulu_curated_headers.txt`

`novaseq/18S/MOTU/18S_cluster_reps_lulu_curated.fa`

### Step 11: Taxonomy subsetting

Script: `Subset_18S_taxonomy.R`

Run command: `Rscript Subset_18S_taxonomy.R`

This script subsets the ensemble taxonomy table to retain only taxonomic assignments for MOTUs that survived LULU curation.

Main final outputs include:

`novaseq/18S/18S_final_tax_table.txt`

## 7. Running the full pipeline

Submit the job: `sbatch run/novaseq_slurm.sh`

Check job status: `squeue -u <PDC-username>`
