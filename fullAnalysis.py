import os
import re
import sys
import glob
import ast
import math
import pandas as pd
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
	#print(all_actual_muts)
	if(isinstance(vcf_file, str)):
		f2 = gzip.open(vcf_file, 'rt')
		d2 = f2.readlines()
		called_muts = set()
		for i in d2: 
			z = i.split('\t')
			if(len(z) == 11 and z[0] in clist):
				tup = (z[0], int(z[1]))
				called_muts.add(tup)
		#print(called_muts)
	else: 
		called_muts = set()
		for ind_file in vcf_file: 
			f2 = gzip.open(ind_file, 'rt')
			d2 = f2.readlines()
			for i in d2: 
				z = i.split('\t')
				if(len(z) == 11 and z[0] in clist):
					tup = (z[0], int(z[1]))
					called_muts.add(tup)
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
		current_tumor_path = data_path+'tumor_{}/*/'.format(i)
		subsamples = glob.glob(current_tumor_path)
		list_of_samples.append(subsamples)
	return list_of_samples
def generateParameters(search_dir):
	info_file = search_dir+'/parameter_list.txt'
	#TODO!
	list_of_parameters = []
	tracker = 0
	with open(info_file, 'rb') as txtfile: 
		for line in txtfile: 
			key, val = line.decode("utf-8").split("\n")[0].split(":")
			if(tracker != 9):
				list_of_parameters.append(val)
			tracker += 1
	return list_of_parameters
def main():
	data_name = sys.argv[1]
	result_file = open('/home/assrivat/simulation_results/results/{}/SNPcallerstats.txt'.format(data_name), 'w')
	samples = getTumorDirectories(data_name)
	list_of_list_vals = []
	for i in range(len(samples)): 
		f = '/home/assrivat/simulation_results/{}/tumor_{}/information_list.txt'.format(data_name, i)  
		aggregate_vcfs = []
		subtuples_list = []
		max_sample_num = 0
		for j in samples[i]:
			r1 = re.compile('samplenum_([0-9]*)')
			sample_num = r1.findall(j)[0]
			the_number = int(sample_num)
			search_dir = '/home/assrivat/simulation_results/{}/tumor_{}/samplenum_{}'.format(data_name, i,j)
			if(int(sample_num) > max_sample_num): 
				max_sample_num = int(sample_num)
			if 'singlecell' in j:
				single_cell_flag = True
				regex = re.compile('singlecell_([0-9]*)')
				single_cell_num = regex.findall(j)[0]
				subtuples_list.append((sample_num,single_cell_num))	
				results_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}_singlecell_{}/snv/strelka/results/variants/somatic.snvs.vcf.gz'.format(data_name, i,sample_num, single_cell_num)
			else: 
				single_cell_flag = False
				single_cell_num = 0
				subtuples_list.append((sample_num,-1))
				results_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}/snv/strelka/results/variants/somatic.snvs.vcf.gz'.format(data_name, i,sample_num)
			aggregate_vcfs.append(results_dir)
			tp,fp,mm =compareVcftoInfo(results_dir, f)
			result_file.write('tumor {}, sample {}, values: {}, {}, {}\n'.format(i,j,tp,fp,mm))
		print(max_sample_num)
		print(subtuples_list)
		for k in range(max_sample_num+1): 
			total_run = [item for item in subtuples_list if item[0] == str(k)]
			subaggregate_list = []
			print(total_run)
			search_dir = '/home/assrivat/simulation_results/{}/tumor_{}/samplenum_{}/'.format(data_name,i,k)
			#TOFIX!
			list_of_parameters = generateParameters(search_dir)
			for item in total_run: 
				if(item[1] == -1): 	
					the_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}/snv/strelka/results/variants/somatic.snvs.vcf.gz'.format(data_name, i,item[0])
				else: 
					the_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}_singlecell_{}/snv/strelka/results/variants/somatic.snvs.vcf.gz'.format(data_name, i,item[0], item[1])
				subaggregate_list.append(the_dir)
			tp,fp,mm= compareVcftoInfo(subaggregate_list,f)
			list_of_parameters.extend([i,k,str(-1),tp,fp,mm])
			list_of_list_vals.append(list_of_parameters)
			result_file.write('tumor {}, sample {}, values: {}, {}, {}\n'.format(i,k,tp,fp,mm))
		tp,fp,mm = compareVcftoInfo(aggregate_vcfs,f)
		new_write_list = ['-1']*13
		new_write_list.extend([i,tp,fp,mm])
		list_of_list_vals.append(new_write_list)
		result_file.write('AGGREGATE TUMOR {}, values: {}, {}, {}\n'.format(i,tp,fp,mm))
		aggregate_vcfs.clear()
		subtuples_list.clear()
	print(list_of_list_vals)
	final_results = pd.DataFrame(list_of_list_vals, columns=['num_leaves', 'dir_conc','cell_pop', 'coverage', 'num_single_cells', 'read_len', 'frag_len', 'paired', 'exome','poisson_time', 'error_rate','tumor_num', 'sample_num', 'aggregate_tumor','true_positive','false_positive','missed_mutation'])
	final_results.to_csv('/home/assrivat/simulation_results/results/{}/resultsSNP.csv'.format(data_name), sep='\t')
if __name__ == '__main__':
	main()
 
