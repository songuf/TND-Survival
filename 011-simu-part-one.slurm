#!/bin/sh
#SBATCH --job-name=simu               # Job name
#SBATCH --array=1-100
#SBATCH --account=xxx       
#SBATCH --qos=xxx-b          
#SBATCH --mail-type=ALL               # Mail events
#SBATCH --mail-user=xxx@xxx.edu       # Where to send email 
#SBATCH --ntasks=1                    # Run on a single machine (node)
#SBATCH --cpus-per-task 1             # Run on a single CPU
#SBATCH --mem=2gb                     # Memory limit
#SBATCH --time=48:00:00               # Time: hrs:min:sec
#SBATCH --output=xxx.Rout             # Output and error log 

pwd; hostname; date 

module load R/4.1 

echo "Running" 

R CMD BATCH --vanilla 031-dissert-simu.R ../Rout/031-dissert-simu.$SLURM_ARRAY_TASK_ID.Rout

date