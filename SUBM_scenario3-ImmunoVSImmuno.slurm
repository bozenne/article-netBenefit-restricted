#!/bin/bash

#SBATCH --job-name=scenario3-ImmunoVSImmuno  # Job name
#SBATCH --array=1-20     # Number of jobs
#SBATCH --partition=long    # Queue. Use long for run time >8h and standard otherwise
#SBATCH --mem=1G # Memory limit, e.g. reserve 1 GB memory 
#SBATCH --output=output/scenario3-ImmunoVSImmuno/slurm-%x-%a-%J.out # Direct output to subdirectory
#SBATCH --error=output/scenario3-ImmunoVSImmuno/slurm-%x-%a-%J.out # Same file for output and error

R CMD BATCH --vanilla BATCH_scenario3-ImmunoVSImmuno.R output/scenario3-ImmunoVSImmuno/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout

## go to directory    ## cd /projects/biostat01/people/hpl802/GPC/article-restricted/
## clean outputs      ## rm -r ./output/scenario3-ImmunoVSImmuno/*
## clean results      ## rm -r ./Results/scenario3-ImmunoVSImmuno/*
## submission command ## sbatch SUBM_scenario3-ImmunoVSImmuno.slurm

## submission output  ## Submitted batch job 1512432 (time using C-c t: fredag 22-06-17 at 18:20) 

## init               ## module module load gcc/11.2.0 R/4.1.2
                      ## cat .bash_profile (in the root directory on the server to visualize the current configuration)
                      ## vi .bash_profile (in the root directory on the server to modify the current configuration)

## commands           ## squeue           : view current jobs 
                      ## squeue -u id     : view current jobs for user id
                      ## scancel          : delete the job with job id 1034 from the queue type
                      ## sinfo            : view state of the queues
                      ## sshare           : check own usage of the server
                        
## documentation      ## from SLURM: https://slurm.schedmd.com/pdfs/summary.pdf: key functions
                      ##           : https://slurm.schedmd.com/sbatch.html (filename pattern)
                      ## from KU   : https://hpc.ku.dk/documentation/slurm.html
                      ##           : https://github.com/aejensen/computing-biostat-ucph/blob/main/README.md
