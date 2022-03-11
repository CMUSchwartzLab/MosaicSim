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
num_samples_list = [8] 
num_tumors_list = [2]
# frag length should be at least 100, and the read length should be less than the frag length
read_len_list = [75,125,250,1000,10000]
frag_len_list = [150,200,400,2000,20000]
#concentration of the clones
alpha_list = [10]
#can be True or False for these
paired_list = [True,False]
WES_list = [True,False]

use_leaf_only = False
error_rate_list = [0.0, 0.0, 0.000001, 0.00001, 0.0001, 0.001]
clone_list = [5]
#rate list define the tumor mutations, can define or change however you want. 
ultralow_rates_list = [1e-15]
low_rates_list = [1e-15, 1e-13, 1e-10]
medium_rates_list = [5e-9, 4e-10, 9e-10, 3e-9]
high_rates_list = [1e-8, 9e-9, 5e-9]
ultrahigh_rates_list = [1e-8,5e-7,9e-7,5e-8,2e-7]
coverage_list = [2,5,10, 15, 30, 50]
pop_list = [8e8] #recommend to keep close to this value
num_single_cell_list = [0,1]
liquid_biopsy_list = [False]
batch_size = 1
LSH_hash = False
kmer_len=50
num_perm=64
thresh=0.33
ctdna_frac_list = [0.96]

#CHANGE THIS!
randomid = 'simulation_grid_search'
#CHANGE THIS TO  WHERE YOU WANT TO STORE THE DATA!
base_working_dir = f'/home/assrivat/simulation_results/{randomid}/'

#reference/normal sequencing parameters
reference_working_dir = base_working_dir + 'reference/'
ref_coverage = random.choice(coverage_list)
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
ref_erate = random.choice(error_rate_list)

#SNV, CNV, DEL, DELSMALL, INVERSION, TRANSLOCATION, BFB, CHROMOTHRIP, CHROMOPLEX, INSERTIONSMALL, KATAEGIS, ANEUPLOIDY
list_of_rates = [ultrahigh_rates_list, ultralow_rates_list, high_rates_list, ultralow_rates_list, high_rates_list, ultralow_rates_list,ultralow_rates_list, ultralow_rates_list, ultralow_rates_list, ultralow_rates_list, ultralow_rates_list, ultralow_rates_list]
 

