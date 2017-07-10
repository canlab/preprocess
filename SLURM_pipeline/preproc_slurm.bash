#!/bin/bash
#SBATCH -p blanca-ics
#SBATCH -n 1
#SBATCH -t 360
#SBATCH --mem=16G

cd /work/ics/data/projects/wagerlab/labdata/collab/Mackey_Duloxetine

matlab -nodisplay -nosplash -nodesktop -r "addpath('/work/ics/data/projects/wagerlab/Resources/spm12'); spm('fmri'); load('preproc_batches/duloxetine_preproc_batch_${SLURM_ARRAY_TASK_ID}.mat'); spm_jobman('initcfg'); spm_jobman('run',{matlabbatch'}); exit;"


