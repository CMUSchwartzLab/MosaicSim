srun -p pool1 --mem-per-cpu 46G -t 2-23:00:00 python3 /home/assrivat/MosaicSim/sim.py
srun -p pool3-bigmem --mem-per-cpu 125G -t 3-23:00:00 python3 /home/assrivat/MosaicSim/sim.py
srun -p rs2 --mem-per-cpu 100G -t 3-22:00:00 python3 /home/assrivat/MosaicSim/sim.py
