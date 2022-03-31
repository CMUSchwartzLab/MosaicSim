#parameters for simulation
import random
import pandas as pd
import string
tab = str.maketrans("ACTG", "TGAC")

#signature parameters
list_of_bases = ['A', 'C', 'T', 'G']
list_of_pairs = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
sig_df = pd.read_csv('./data/signatures.txt', sep = '\t', header = None)
sig_df = sig_df.iloc[1:]
sig_df = sig_df.iloc[:,1:]
signatures_matrix = sig_df.to_numpy()
num_signatures = 78
use_signatures = False
signature_alpha = 10
signature_distributions = [float(1/num_signatures)]*num_signatures

#genome and exome delimiters
full_genome = './data/hg38.fa'
EXON_FILE = './data/exonsegments.txt'

#parameters for the full simulation. The program will sample randomly from each list
num_samples_list = [1] 
num_tumors_list = [1]
# frag length should be at least 100, and the read length should be less than the frag length
read_len_list = [75, 150, 500, 2000]
frag_len_list = [150,200, 1000, 3000]
#concentration of the clones
alpha_list = [10]
#can be True or False for these
paired_list = [True]
WES_list = [False, True]

use_leaf_only = False
error_rate_list = [0.0, 0.0, 0.0, 0.00001,0.00001,0.0001, 0.001]
clone_list = [5]
#rate list define the tumor mutations, can define or change however you want. 
ultralow_rates_list = [1e-15]
low_rates_list = [1e-15, 1e-13, 1e-10]
medium_rates_list = [5e-9, 4e-10, 9e-10, 3e-9]
high_rates_list = [1e-8, 9e-9, 5e-9]
ultrahigh_rates_list = [5e-8,5e-7,6e-7,5e-8,1e-7]
coverage_list = [2, 5, 10, 15,15, 25,25, 30,30]
pop_list = [8e8] #recommend to keep close to this value
num_single_cell_list = [0]
liquid_biopsy_list = [False]
batch_size = 1
subblock_size = 1
LSH_hash = False #If false, doesn't explore alternate reads for WES and uses breakpoints in EXON_FILE
LSHStringNum = 10000000 #G_N in the paper -- the sampling number
kmer_len=50
num_perm=64
thresh=0.33
ctdna_frac_list = [0.96] #fraction of real tumor reads in liquid biopsy

#CHANGE THIS TO  WHERE YOU WANT TO STORE THE DATA!
storage_dir = '/home/assrivat/simulation_results/'

#reference/normal sequencing parameters
ref_coverage = 30 #random.choice(coverage_list)
ref_clones = 5  # this doesnt matter
read_length_index = random.randint(0, len(read_len_list)-1)
ref_read_len = read_len_list[read_length_index]
ref_frag_len = frag_len_list[read_length_index]
ref_tot_nodes = 2*ref_clones-1
ref_root_node = ref_tot_nodes-1
ref_int_nodes = ref_root_node-1
ref_alpha = random.choice(alpha_list)
ref_paired = True
ref_WES = False
ref_erate = 0.0 #random.choice(error_rate_list)

if(random.random() < 0.5): 
	random_list = ultralow_rates_list
else: 
	random_list = ultrahigh_rates_list
#SNV, CNV, DEL, DELSMALL, INVERSION, TRANSLOCATION, BFB, CHROMOTHRIP, CHROMOPLEX, INSERTIONSMALL, KATAEGIS, ANEUPLOIDY
list_of_rates = [high_rates_list, ultralow_rates_list, medium_rates_list, ultralow_rates_list, medium_rates_list, ultralow_rates_list,ultralow_rates_list, ultralow_rates_list, ultralow_rates_list, ultralow_rates_list, ultralow_rates_list, ultralow_rates_list]
 

