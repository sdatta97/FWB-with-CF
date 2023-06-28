#!/bin/bash

#SBATCH --mem=4GB
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=./logs/%a_%A.out
#SBATCH --error=./logs/%a_%A.err


#SBATCH --job-name=FWB
#SBATCH --mail-type=END
#SBATCH --mail-user=sd3927@nyu.edu


module purge
module load matlab/2022a

cd /scratch/sd3927/FWB-CF-v2/Simulation

echo "$SLURM_ARRAY_TASK_ID"

matlab -batch "SimulationMain, exit"

