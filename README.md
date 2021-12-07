# MosaicSim
# Dependencies
The simulator requires the following python packages to run: pandas, numpy, biopython, tskit, msprime, glob

Each of the packages can be easily installed using pip, and running "pip install PACKAGE" for each of the packages above.

Alternatively, one could install these packages in a conda environment
# Execution
The user can adjust the parameters of the simulation in the "parameters.py" file. Be sure to change the base_working_dir to the directory you want to store your data, and change the other parameters as desired. Descriptions of each of the parameters are included as comments in the parameters.py script and in the paper.

To run the simulation, simply execute "python3 sim.py" in the repository.
# System Requirements
The program should take less than 30G of memory and only one core per execution. The storage requirements depend on the parameters of the simulation. Run time was approximately 4 and a half hours for 3, 30x WGS samples.
