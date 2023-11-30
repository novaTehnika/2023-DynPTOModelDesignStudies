#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=12500M
#SBATCH -t 16:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simmo536@umn.edu
#SBATCH -p small
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

cd ~/2023-DynPTOModelDesignStudies
module load matlab
matlab -nodisplay -r \
"iVar = ${SLURM_ARRAY_TASK_ID}; \
display(['iVar = ',num2str(iVar)]); \
SS = ${SS}; \
display(['SS = ',num2str(SS)]); \
addpath('Utilities'); \
parSafeStartSlurm; \
study_parPTO_LPaccum; \
rmdir(storage_folder)"

# Commands to use
# sbatch --export=SS=1 --array=1-1800 ~/2023-DynPTOModelDesignStudies/study_parPTO_LPaccum.sh
# dos2unix  study_parPTO_LPaccum.sh
