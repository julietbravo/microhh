#!/bin/bash
#SBATCH --account=MAQ440
#SBATCH --time=01:00:00
#SBATCH --mem=4096
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --output=out_%j.txt
#SBATCH --error=err_%j.txt
#SBATCH --job-name=case_name
#SBATCH --partition=ESG_Std
#SBATCH --mail-type=END
#SBATCH --mail-user=dummy@dummy.org	

mpirun ./microhh init patch
mpirun ./microhh run patch
