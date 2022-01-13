srun -p rs2 -c 40 --mem-per-cpu 12G python3 python_pipeline_v1.py lowcovnosnv 0 0 0 40 strelka cnvkit delly
srun -p pool3-bigmem -c 8 --mem-per-cpu 14G python3 python_pipeline_v1.py lowcovnosnv 1 0 0 8 strelka cnvkit delly
