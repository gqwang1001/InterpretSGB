#!/bin/bash
#$ -cwd -V -b y
#$ -N simulation
#$ -o out.txt
#$ -e err.txt
#$ -pe threads 1
#$ -A mytestParallelR
#$ -l mem_reserve=3G -l h_vmem=3G
#$ -P long

source /SFS/product/Modules/default/init/bash
module purge
module load R/3.4.3
module load openmpi/1.10.7_sge
module list
which orted
cd /SFS/scratch/zhapingy/sim32/code
pwd
R CMD BATCH --no-save --no-restore $1 ${SGE_TASK_ID}.Rout