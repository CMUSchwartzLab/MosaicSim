# MosaicSim
# Dependencies
The user should first download python3 if they do not already have it.
The simulator requires the following python packages to run: pandas, numpy, biopython, tskit, msprime, glob

Each of the packages can be easily installed using pip, and running "pip3 install PACKAGE" for each of the packages above.

Alternatively, one could install these packages in a conda environment
# Data Dependency
The simulator requires a genome from which to simulate from. We use the hg38 reference genome, which can be downloaded by first running "wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" and then "gunzip hg38.fa.gz". Be sure the unzipped hg38.fa file is in the MosaicSim/data directory to avoid a runtime error. The user could also swap this genome for another reference genome by simply updating the "full_genome" parameter in the parameters.py file, as long as it is in standard FASTA format. If the alternative FASTA file does not follow the standard labeling conventions of hg38, the user may have to edit the chromosome preprocessing code in sim.py as desired.

# Execution
Each run of the simulator is defined by first adjusting the study design and simulation parameters in "parameters.py" and then running the "sim.py" script. In particular, be sure to change the storage_dir variable to the directory you want to store your data, and change the other parameters as desired to your biological and study design choices. Descriptions of each of the parameters are included as comments in the parameters.py script and in the paper.

To run the simulation, simply execute "python3 sim.py [DIR_NAME]" in the repository. DIR_NAME defines the name of the directory under storage_dir where all your data will be stored for the simulation run.

# System Requirements
The program should, for standard sequencing settings, take less than 40G of memory and only one core per execution. The storage requirements depend on the parameters of the simulation. Run time was approximately 3 and a half hours for 3, 30x WGS samples.
