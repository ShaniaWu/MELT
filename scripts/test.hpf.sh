#!/bin/bash
#SBATCH --job-name=report_test
#SBATCH --time=100:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out
#SBATCH --constraint=CentOS7
#SBATCH --mail-type=ALL

module load anaconda3

source activate /hpf/largeprojects/ccmbio/swu/MELT/scripts/env


python MELT_report_mcouse.py
