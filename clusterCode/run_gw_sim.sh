#!/bin/bash
#$ -cwd -V -b y
#$ -N simulation
#$ -o out.txt
#$ -e err.txt
#$ -pe threads 1-10
#$ -A mytestParallelR
#$ -l mem_reserve=3G -l h_vmem=3G -l fourperhost=4
#$ -P long

# source /SFS/product/Modules/default/init/bash
module purge
module load R/3.4.3
module load openmpi/1.10.7_sge
module list
which orted

Rscript gw_run_sim.R $1 $2