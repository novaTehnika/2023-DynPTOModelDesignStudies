#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=12500M
#SBATCH -t 24:00:00 #48:00:00 # 4 workers, 8hs per sim, 14 sims per job, 1.5 S.F.
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
nWorkers = ${SLURM_NTASKS}-1; \
parSafeStartSlurm; \
study_parPTO_accum_woRV; \
rmdir(storage_folder)"

# Commands to use
# sbatch --export=SS=2 --array=1-500 ~/2023-DynPTOModelDesignStudies/study_parPTO_accum_woRV.sh
# dos2unix  study_parPTO_accum_woRV.sh

