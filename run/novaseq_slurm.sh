#!/bin/bash -l

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A naiss2026-4-900

# The name of the script is myjob
#SBATCH -J openarms-joren-18S-2022-2025

# The partition
#SBATCH -p main ##main is the part of the cluster where you run it, main is 24h max

# 24 hours wall clock time will be given to this job
#SBATCH -t 24:00:00

# Number of nodes
#SBATCH --nodes=1

#SBATCH --mail-user=Joren.VAN.GANSBERGHE@student.howest.be
#SBATCH --mail-type=ALL

## Set the names for the error and output files. 
## It can be smart to set a path to these to your project directory, which you can do by adding that path right after the '=' sign
#SBATCH --error=/cfs/klemming/projects/supr/naiss2025-23-46/Joren/Project_18S/sbatch_jobs/job.%J.err
#SBATCH --output=/cfs/klemming/projects/supr/naiss2025-23-46/Joren/Project_18S/sbatch_jobs/job.%J.out

## module purge, then load necessary modules

module purge 
module load PDC/24.11
module load apptainer/1.4.0-cpeGNU-24.11

WORKDIR=/cfs/klemming/projects/supr/naiss2025-23-46/Joren/Project_18S;
TMPDIR=/cfs/klemming/scratch/j/jorenvgb;
SANDBOX="$TMPDIR/temparms-sandbox"

apptainer exec --bind "$WORKDIR":/work --bind "$TMPDIR/conda_envs":/envs --bind /cfs/klemming/projects/supr/naiss2025-23-46/ARMS_Genoscope_data_full/projet_DBB/ARMS-Macro:/cfs/klemming/projects/supr/naiss2025-23-46/ARMS_Genoscope_data_full/projet_DBB/ARMS-Macro "$SANDBOX" bash /work/run/novaseq_run.sh



