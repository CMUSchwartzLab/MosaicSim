#!/bin/bash

## select the pool
#SBATCH -p pool1

## select the walltime
#SBATCH -t 71:59:59

## select the output dir
#SBATCH -o /home/haoyunl/projects/SeqSimulator/code/bgkbpsu/%j.out

## exclude those nodes
#SBATCH --exclude=compute-0-[22-28,30]
#SBATCH --exclude=compute-1-[7-8,20-22,31-34,36]

#SBATCH -c 8
#SBATCH --mem-per-cpu 5G


python python_pipeline_v1.py bgkbpsu 0 0 0 8 strelka
