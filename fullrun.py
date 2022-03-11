import os
from parameters import *


os.system("srun -p pool1 --mem-per-cpu 45G -t 3-00:00:00 python3 /home/assrivat/MosaicSim/sim.py; srun -p pool1 -c 7 --mem-per-cpu 2200 -t 3-00:00:00 python3 /home/assrivat/MosaicSim/pipeline/python_pipeline_v1.py {} 1 0 0 7 strelka None delly hg38_no_alt".format(randomid))
#os.system("srun -p interactive --mem-per-cpu 10G -t 5:00:00 python3 /home/assrivat/MosaicSim/sim.py; srun -p interactive -c 4 --mem-per-cpu 2200 -t 5:00:00 python3 /home/assrivat/MosaicSim/pipeline/python_pipeline_v1.py {} 1 0 0 4 None None delly fakegenome; srun -p interactive --mem-per-cpu 10G python3 /home/assrivat/MosaicSim/fullAnalysis.py {}".format(randomid, randomid))

#os.system("srun -p rs2 --mem-per-cpu 150G -t 5-22:00:00 python3 /home/assrivat/MosaicSim/sim.py; srun -p rs2 -c 128 --mem-per-cpu 2200 -t 7-22:00:00 python3 /home/assrivat/MosaicSim/pipeline/python_pipeline_v1.py {} 1 0 0 128 strelka None delly hg38_no_alt; srun -p rs2 --mem-per-cpu 100G python3 /home/assrivat/MosaicSim/fullAnalysis.py {}".format(randomid, randomid))
