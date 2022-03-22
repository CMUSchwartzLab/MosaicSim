import os
import re
import sys
import glob
import ast
import math
import pandas as pd
import gzip
from collections import defaultdict
from operator import itemgetter
from intervaltree import Interval, IntervalTree
clist = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13','chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19','chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8','chr9', 'chrX', 'chrY']
def intervalSetMetric(interval_a, interval_b): 
	#takes two interval sets and computes a difference between then
	groups = defaultdict(list)
	for i in clist: 
		interval_a.append((-2,-1,i))
		interval_b.append((-2,-1,i))
	for t in interval_a:
		groups[t[-1]].append(t)
	lol1 = list(groups.values())
	groups2 = defaultdict(list)
	for z in interval_b: 
		groups2[z[-1]].append(z)
	lol2 = list(groups2.values())
	lol2.sort(key= lambda x: x[0][2])
	lol1.sort(key=lambda x: x[0][2])
	metric_list = []
	event_num = 0
	for i in range(len(lol1)):
		comp1 = lol1[i]
		comp2 = lol2[i]
		comp1.sort(key = lambda x: x[0])
		comp2.sort(key = lambda x: x[0])
		intervals1 = []
		intervals2 = []
		for j in comp1[1:]: 
			intervals1.append((j[0], j[1]))
		for k in comp2[1:]: 
			intervals2.append((k[0], k[1]))
		intervals1.sort(key = lambda x: x[0])
		intervals2.sort(key = lambda x: x[0])
		t = IntervalTree.from_tuples(intervals1)
		total_intervals = len(intervals2)
		current_count = 0
		for k in intervals2:
			if(len(t[k[0]:k[1]]) > 0):
				current_count += 1
		num1 = len(intervals1)
		num2 = len(intervals2)
		total_comparate = num1+num2
		event_num += total_comparate
		if(total_intervals == 0): 
			metric_list.append([0,0])
		else:
			metric_list.append([current_count/total_intervals, total_comparate])
	run_metric = 0 
	for i in metric_list: 
		new_value = i[0]*(i[1]/event_num)
		run_metric += new_value
	return run_metric
		
def compareVcftoInfoSV(vcf_file, info_file, mut_file, sv_type): 
	f = open(info_file, 'r')
	m = open(mut_file, 'r')
	m2 = m.readlines()[0]
	d = f.readlines()[0]
	actual_m = ast.literal_eval(m2)
	actual_d = ast.literal_eval(d)
	if(sv_type == 4): 
		checker_flag = "<INV>"
	else: 
		checker_flag = '<DEL>'
	all_actual_muts = []
	counter = 0
	for i in range(len(actual_d)):
		for k in range(len(actual_d[i])):
			if(int(actual_m[i][k]) == sv_type):
				tup = (int(actual_d[i][k][0]), int(actual_d[i][k][1]), actual_d[i][k][-1])
				all_actual_muts.append(tup)
	#print(all_actual_muts)
	if(isinstance(vcf_file, str)):
		f2 = open(vcf_file, 'rt')
		d2 = f2.readlines()
		called_muts = []
		for i in d2: 
			z = i.split('\t')
			if(len(z) > 9 and z[0] in clist):
				if(z[4] == checker_flag):
					string_cont = z[7]
					match = re.findall('END=([0-9]*)', string_cont)
					regex_parsed = int(match[0])
					tup = (int(z[1]), regex_parsed, z[0])
					called_muts.append(tup)
	else: 
		called_muts = []
		for ind_file in vcf_file: 
			f2 = open(ind_file, 'rt')
			d2 = f2.readlines()
			for i in d2: 
				z = i.split('\t')
				if(len(z) > 9 and z[0] in clist):
					if(z[4] == checker_flag):
						string_cont = z[7]
						match = re.findall('END=([0-9]*)', string_cont)
						regex_parsed = int(match[0])
						tup = (int(z[1]), regex_parsed, z[0])
						called_muts.append(tup)
	metrika = intervalSetMetric(called_muts, all_actual_muts)
	return metrika

def compareVcftoInfoSNV(vcf_file, info_file, mut_file):
	f = open(info_file, 'r')
	m = open(mut_file, 'r')
	m2 = m.readlines()[0]
	d = f.readlines()[0]
	actual_m = ast.literal_eval(m2)
	actual_d = ast.literal_eval(d)
	all_actual_muts = set()
	counter = 0
	for i in range(len(actual_d)):
		for k in range(len(actual_d[i])):
			if(actual_m[i][k] == 0):
				tup = (actual_d[i][k][-1], actual_d[i][k][-2])
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

def checkPaired(search_dir):
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
    if(list_of_parameters[7] == ' True'): 
        return True
    else: 
         return False


def main():
	
	data_name = sys.argv[1]
	result_file = open('/home/assrivat/simulation_results/results/{}/SNPcallerstats.txt'.format(data_name), 'w')
	sv_resultfile = open('/home/assrivat/simulation_results/results/{}/SVcallerstats.txt'.format(data_name),'w')
	samples = getTumorDirectories(data_name)
	list_of_list_vals = []
	sv_lol_vals  = []
	for i in range(len(samples)): 
		f = '/home/assrivat/simulation_results/{}/tumor_{}/information_list.txt'.format(data_name, i)  
		m = '/home/assrivat/simulation_results/{}/tumor_{}/mutation_list.txt'.format(data_name, i)
		aggregate_vcfs = []
		subtuples_list = []
		aggregate_sv_vcfs = []
		max_sample_num = 0
		for j in samples[i]:
			r1 = re.compile('samplenum_([0-9]*)')
			sample_num = r1.findall(j)[0]
			the_number = int(sample_num)
			search_dir = '/home/assrivat/simulation_results/{}/tumor_{}/samplenum_{}'.format(data_name, i,the_number)
			if(checkPaired(search_dir)): 
				sv_flag = True
			else: 
				sv_flag = False
			if(int(sample_num) > max_sample_num): 
				max_sample_num = int(sample_num)
			if 'singlecell' in j:
				single_cell_flag = True
				regex = re.compile('singlecell_([0-9]*)')
				single_cell_num = regex.findall(j)[0]
				subtuples_list.append((sample_num,single_cell_num))	
				results_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}_singlecell_{}/snv/strelka/results/variants/somatic.snvs.vcf.gz'.format(data_name, i,sample_num, single_cell_num)
				if(sv_flag): 
					sv_results_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}_singlecell_{}/sv/delly/tumorB.vcf'.format(data_name, i,sample_num, single_cell_num)
			else: 
				single_cell_flag = False
				single_cell_num = 0
				subtuples_list.append((sample_num,-1))
				results_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}/snv/strelka/results/variants/somatic.snvs.vcf.gz'.format(data_name, i,sample_num)
				if(sv_flag): 
					sv_results_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}/sv/delly/tumorB.vcf'.format(data_name, i,sample_num)
			aggregate_vcfs.append(results_dir)
			if(sv_flag): 
				aggregate_sv_vcfs.append(sv_results_dir)
				deletion_metric = compareVcftoInfoSV(sv_results_dir,f, m, 2) #TODO: might be a type flag here
				inversion_metric = compareVcftoInfoSV(sv_results_dir,f,m, 4)
				sv_resultfile.write('tumor {}, sample {}, deletion_metric {}, inversion_metric {}\n'.format(i,j,deletion_metric, inversion_metric))
			tp,fp,mm =compareVcftoInfoSNV(results_dir, f,m)
			result_file.write('tumor {}, sample {}, values: {}, {}, {}\n'.format(i,j,tp,fp,mm))
		for k in range(max_sample_num+1): 
			total_run = [item for item in subtuples_list if item[0] == str(k)]
			subaggregate_list = []
			search_dir = '/home/assrivat/simulation_results/{}/tumor_{}/samplenum_{}/'.format(data_name,i,k)
			if(checkPaired(search_dir)): 
				sv_flag = True
			else: 
				sv_flag = False

			#TOFIX!
			list_of_parameters = generateParameters(search_dir)
			for item in total_run: 
				if(item[1] == -1): 	
					the_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}/snv/strelka/results/variants/somatic.snvs.vcf.gz'.format(data_name, i,item[0])
				else: 
					the_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}_singlecell_{}/snv/strelka/results/variants/somatic.snvs.vcf.gz'.format(data_name, i,item[0], item[1])
				subaggregate_list.append(the_dir)
			tp,fp,mm= compareVcftoInfoSNV(subaggregate_list,f,m)
			list_of_parameters.extend([i,k,str(-1),tp,fp,mm])
			list_of_list_vals.append(list_of_parameters)
			result_file.write('tumor {}, sample {}, values: {}, {}, {}\n'.format(i,k,tp,fp,mm))
			if(sv_flag):
				sv_param_list = generateParameters(search_dir)
				sv_subaggregate_list = []
				for item in total_run: 
					if(item[1] == -1): 
						the_sv_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}/sv/delly/tumorB.vcf'.format(data_name, i, item[0])
					else: 
						the_sv_dir = '/home/assrivat/simulation_results/results/{}/tumor_{}/samplenum_{}_singlecell_{}/sv/delly/tumorB.vcf'.format(data_name, i, item[0], item[1])
				sv_subaggregate_list.append(the_sv_dir)
				deletion_metric = compareVcftoInfoSV(sv_subaggregate_list, f,m, 2)
				inversion_metric = compareVcftoInfoSV(sv_subaggregate_list, f,m,4)
				sv_resultfile.write('tumor {}, sample {}, del met {}, inv met {}\n'.format(i,k,deletion_metric, inversion_metric))
				sv_param_list.extend([i,k,deletion_metric, inversion_metric])
				sv_lol_vals.append(sv_param_list)
		tp,fp,mm = compareVcftoInfoSNV(aggregate_vcfs,f,m)
		new_write_list = ['-1']*13
		new_write_list.extend([i,tp,fp,mm])
		list_of_list_vals.append(new_write_list)
		del_tot = compareVcftoInfoSV(aggregate_sv_vcfs, f, m, 2)
		inv_tot = compareVcftoInfoSV(aggregate_sv_vcfs, f, m, 4)
		svwl = ['-1']*13
		svwl.extend([i,del_tot,inv_tot])
		sv_lol_vals.append(svwl)
		sv_resultfile.write('AGGREGATE TUMOR {}, del val {} inv value {}\n'.format(i,del_tot, inv_tot))
		result_file.write('AGGREGATE TUMOR {}, values: {}, {}, {}\n'.format(i,tp,fp,mm))
		aggregate_sv_vcfs.clear()
		aggregate_vcfs.clear()
		subtuples_list.clear()
	final_results = pd.DataFrame(list_of_list_vals, columns=['num_leaves', 'dir_conc','cell_pop', 'coverage', 'num_single_cells', 'read_len', 'frag_len', 'paired', 'exome','poisson_time', 'error_rate','tumor_num', 'sample_num', 'aggregate_tumor','true_positive','false_positive','missed_mutation'])
	sv_final_results = pd.DataFrame(sv_lol_vals, columns=['num_leaves', 'dir_conc','cell_pop', 'coverage', 'num_single_cells', 'read_len', 'frag_len', 'paired', 'exome','poisson_time', 'error_rate','tumor_num', 'sample_num', 'aggregate_tumor','del metric', 'inv metric'])
	final_results.to_csv('/home/assrivat/simulation_results/results/{}/resultsSNP.csv'.format(data_name), sep='\t')
	sv_final_results.to_csv('/home/assrivat/simulation_results/results/{}/resultsSV.csv'.format(data_name), sep='\t')
if __name__ == '__main__':
	main()
 
