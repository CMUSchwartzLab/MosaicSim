import os
from parameters import *

os.system("srun -p pool1 --mem-per-cpu 42G -t 2-23:00:00 python3 /home/assrivat/MosaicSim/sim.py; srun -p pool1 -c 16 --mem-per-cpu 3000 -t 2-23:00:00 python3 /home/assrivat/MosaicSim/pipeline/python_pipeline_v1.py {} 1 0 0 16 strelka None delly hg38_no_alt; srun -p pool1 --mem-per-cpu 45G python3 /home/assrivat/MosaicSim/fullAnalysis.py {}".format(randomid, randomid))
