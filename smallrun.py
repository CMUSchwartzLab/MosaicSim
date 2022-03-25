import os
import subprocess
import sys
from multiprocessing import Pool
useid = sys.argv[2]
partition = sys.argv[1]
nodes = int(sys.argv[3])
if(partition == 'rs2'):
	mem = 100
	threads = 128
	mempercpu = 2300
	timestring = "7-23:00:00"
if(partition == 'bridges'): 
	mem = 100
	threads = 128
	mempercpu = 2500
	timestring = "2-23:00:00"
if(partition == 'pool1'): 
	mem = 45
	threads = 16
	mempercpu = 2500
	timestring = "2-23:00:00"


def oneRun(runnum): 
	step1= subprocess.run("srun -p {} --mem-per-cpu {}G -t {} python3 /home/assrivat/MosaicSim/sim.py {}{}".format(partition,mem, timestring, useid,runnum), shell=True, check =True)
	step2 = subprocess.run(" srun -p {} -c {} --mem-per-cpu {} -t {} python3 /home/assrivat/MosaicSim/pipeline/python_pipeline_v1.py {}{} 1 0 0 {} strelka None delly hg38_no_alt".format(partition, threads, mempercpu,timestring, useid, runnum,threads), shell = True, check = True)
	step3 = subprocess.run("srun -p {} --mem-per-cpu 32G -t {} python3 /home/assrivat/MosaicSim/fullAnalysis.py {}{}".format(partition,timestring, useid,runnum), shell = True, check = True)

if __name__ == '__main__':
	total_list = list(range(nodes))
	with Pool(processes = nodes) as p:
		p.map(oneRun, total_list)
	sys.exit()
