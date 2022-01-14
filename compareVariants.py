import os
import sys
import ast
import math
import gzip
clist = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13','chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19','chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8','chr9', 'chrX', 'chrY']
def compareVcftoInfo(vcf_file, info_file):
	f = open(info_file, 'r')
	d = f.readlines()[0]
	actual_d = ast.literal_eval(d)
	all_actual_muts = set()
	counter = 0
	for i in actual_d:
		for k in actual_d[i]:
			tup = (k[-1], k[-2])
			all_actual_muts.add(tup)
	print(all_actual_muts)
	f2 = gzip.open(vcf_file, 'rt')
	d2 = f2.readlines()
	called_muts = set()
	for i in d2: 
		z = i.split('\t')
		if(len(z) == 11 and z[0] in clist):
			tup = (z[0], int(z[1]))
			called_muts.add(tup)
	print(called_muts)
	intersection_set = called_muts.intersection(all_actual_muts)
	uncalled_muts = all_actual_muts.difference(called_muts)
	false_called_muts = called_muts.difference(all_actual_muts)
	true_positives = len(intersection_set)
	missed_muts = len(uncalled_muts)
	false_positives = len(false_called_muts)
	return true_positives, false_positives, missed_muts
def getTumorDirectories(data_name): 
	data_path = '/home/assrivat/simulation_results/{}/'.format(data_name)
	total_num_tumors = sum([os.path.isdir(data_path+i) for i in os.listdir(data_path)])-1
	list_of_samples = []
	for i in range(total_num_tumors): 
		current_tumor_path = data_path+'tumor_{}/'.format(i)
		num_samples_i = sum(os.path.isdir(current_tumor_path+j) for j in os.listdir(current_tumor_path))
		list_of_samples.append(num_samples_i)
	return list_of_samples

def main():
	data_name = sys.argv[1]
	result_file = open('/home/assrivat/simulation_results/results/{}/SNPcallerstats.txt'.format(data_name), 'w')
	samples = getTumorDirectories(data_name)
	for i in range(len(samples)): 
		for j in range(samples[i]):
			results_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}/snv/strelka/results/variants/somatic.snvs.vcf.gz'.format(data_name, i,j)
			f = '/home/assrivat/simulation_results/{}/tumor_{}/information_list.txt'.format(data_name, i)   
			tp,fp,mm =compareVcftoInfo(results_dir, f)
			result_file.write('tumor {}, sample {}, values: {}, {}, {}\n'.format(i,j,tp,fp,mm))
if __name__ == '__main__':
	main()
 
