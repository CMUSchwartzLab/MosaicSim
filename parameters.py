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
full_genome = './data/fakegenome.fa'
EXON_FILE = './data/exonsegments.txt'

#parameters for the full simulation. The program will sample randomly from each list 
num_samples_list = [5] 
num_tumors_list = [1]
# frag length should be at least 100, and the read length should be less than the frag length
read_len_list = [125]
frag_len_list = [200]
#concentration of the clones
alpha_list = [10]
#can be True or False for these
paired_list = [True]
WES_list = [False]

use_leaf_only = True
error_rate_list = [0.0]
clone_list = [5]
#rate list define the tumor mutations, can define or change however you want. 
ultralow_rates_list = [1e-15]
low_rates_list = [1e-15, 1e-13, 1e-10]
medium_rates_list = [5e-9, 4e-10, 9e-10, 3e-9]
high_rates_list = [1e-8, 9e-9, 5e-9]
coverage_list = [100]
pop_list = [8e8] #recommend to keep close to this value
num_single_cell_list = [0]
liquid_biopsy = False
batch_size = 1
ctdna_frac_list = [0.96]
randomid = ''.join(random.choices(string.ascii_lowercase, k=7))

#CHANGE THIS!
randomid = 'faketestdata'
#CHANGE THIS TO  WHERE YOU WANT TO STORE THE DATA!
base_working_dir = f'/home/assrivat/simulation_results/{randomid}/'

#reference sequencing parameters
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
list_of_rates = [medium_rates_list, ultralow_rates_list, ultralow_rates_list, ultralow_rates_list, ultralow_rates_list, ultralow_rates_list,ultralow_rates_list, ultralow_rates_list, ultralow_rates_list, ultralow_rates_list, ultralow_rates_list, ultralow_rates_list]
 

