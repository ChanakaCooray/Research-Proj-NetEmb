#!/bin/bash
#PBS -q default
#PBS -N Job_GenerateSumMatrix_500k
##serial jobs: only 1 processor core is requested
#PBS -l select=1:mem=8gb:ncpus=4
#PBS -l walltime=08:00:00
##replace "x-ccast-prj" below with "x-ccast-prj-[your sponsor's project group]"
#PBS -W group_list=x-ccast-prj-luliu
#
#
cd $PBS_O_WORKDIR
#
python ../src/preprocessing/GenerateSumMatrix-new.py --data /gpfs1/scratch/chanaka.cooray/jobs/work/output/1CDX1 --output ../test/output/1CDX1 --bin-size "500k" --config-file ../metadata/mm9_chrom_sizes.txt
python ../src/preprocessing/GenerateSumMatrix-new.py --data /gpfs1/scratch/chanaka.cooray/jobs/work/output/1CDX2 --output ../test/output/1CDX2 --bin-size "500k" --config-file ../metadata/mm9_chrom_sizes.txt
python ../src/preprocessing/GenerateSumMatrix-new.py --data /gpfs1/scratch/chanaka.cooray/jobs/work/output/1CDX3 --output ../test/output/1CDX3 --bin-size "500k" --config-file ../metadata/mm9_chrom_sizes.txt
python ../src/preprocessing/GenerateSumMatrix-new.py --data /gpfs1/scratch/chanaka.cooray/jobs/work/output/1CDX4 --output ../test/output/1CDX4 --bin-size "500k" --config-file ../metadata/mm9_chrom_sizes.txt
#

