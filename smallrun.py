import os
import sys

useid = sys.argv[1]
os.system("srun -p rs2 --mem-per-cpu 55G -t 7-23:00:00 python3 /home/assrivat/MosaicSim/sim.py {}; srun -p rs2 -c 20 --mem-per-cpu 3000 -t 7-23:00:00 python3 /home/assrivat/MosaicSim/pipeline/python_pipeline_v1.py {} 1 0 0 20 strelka None delly hg38_no_alt; srun -p rs2 --mem-per-cpu 32G python3 /home/assrivat/MosaicSim/fullAnalysis.py {}".format(useid, useid,useid))
