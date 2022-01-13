#!/bin/bash

## select the pool
#SBATCH -p pool1

## select the walltime
#SBATCH -t 71:59:59

## select the output dir
#SBATCH -o /home/haoyunl/projects/SeqSimulator/code/kjwcrae/%j.out

## exclude those nodes
#SBATCH --exclude=compute-0-[22-28,30]
#SBATCH --exclude=compute-1-[7-8,20-22,31-34,36]

#SBATCH -c 10
#SBATCH --mem-per-cpu 4G


python python_pipeline_v1.py 0 0 0 10 freebayes