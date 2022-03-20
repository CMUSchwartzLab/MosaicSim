#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 09:54:33 2021

@author: leo, arjunsrivatsa
"""

'''
integrate pipeline in to Python?
'''
import os, errno
import subprocess
import sys
import re
import time
import pdb
import glob
import pickle
import collections

'''
1. parse the parameter text, retrieve the read length and Seq data type
2. choose the aligner based on the read length
3. run bash from with Python

folder structure:
random 7 letters for parent folder (e.g kjwcrae)

kjwcrae
    - reference
        parameter_list.txt
        reference fasta
    - tumor_0
        information_list.txt
        mutation_list.txt
        tree_sequence.tree
        - samplenum_0
            parameter_list.txt
            tumor.fasta
    - tumor_1
        ...
        
1. test on bwa, bowtie2, minimap2, make sure the run can be done as previously 
2. test on samtools sort and index
3. test on SNVs caller

4. for each parent folder, analyze all the samples inside it -- think more about this.
'''
def checkpath(path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    return


def get_info(txt):
    '''
    retrieve parameters information and store in pkl file for all
    '''
    dic = {}
    with open(txt, 'rb') as txtfile:
        for line in txtfile:
            key, val = line.decode("utf-8").split("\n")[0].split(":")
            dic[key] = val
    return dic
            
        
def alignment(input_file, output_file, ref, aligner='bwa', threads=4):
    '''
    function for alignment, using the specific aligner
    input_file: a list, single-end or paired-end sequencing data
    output_file: string, output file in SAM format
    ref: string, define the reference genome file
    '''
    prefix = output_file.split('.')[0]
    bam_name = '%s.sorted.bam.bai' % prefix
    '''
    if os.path.exists(bam_name): # check bam index file for less strict checking
        print("Aligned File already exists, quit...")
        exit()
    '''    
    if aligner == 'bwa':
        print("Aligning using BWA...")
        if len(input_file) == 2: #paired-end
            cmd = "bwa mem -t %s %s.fa %s %s > %s" % \
                (threads, ref, input_file[0], input_file[1], output_file)
        else: #single end
            cmd = "bwa mem -t %s %s.fa %s > %s" % \
                (threads, ref, input_file[0], output_file)
        os.system(cmd)
        
    elif aligner == 'bowtie2':
        print("Aligning using BOWTIE2...")
        if len(input_file) == 2: #paired-end
            cmd = 'bowtie2 --threads %s -x %s_index -1 %s -2 %s -S %s' % \
                 (threads, ref, input_file[0], input_file[1], output_file)
        else:
            cmd = 'bowtie2 --thread %s -x %s_index -U %s -S %s' % \
                 (threads, ref, input_file[0], output_file)
        #pdb.set_trace()
        os.system(cmd)
    
    elif aligner == 'minimap2':
        print("Aligning using MINIMAP2...")
        if len(input_file) == 2:
            cmd = "minimap2 -t %s -a %s.fa %s %s > %s" % \
                (threads, ref, input_file[0], input_file[1], output_file)
        else:
            cmd = "minimap2 -t %s -a %s.fa %s > %s" % \
                (threads, ref, input_file[0], output_file) 
        os.system(cmd)
    else:
        print("Only bwa, bowtie2, minimap2 are available, existing...")
        exit()
        
        
    
    
def samtools_sort_index(input_file, threads=4):
    '''
    use samtools for sort and index
    '''
    prefix = input_file.split('.')[0]
    bam_name = '%s.bam' % prefix
    bam_sort_name = '%s.sorted.bam' % prefix
    
    # 1. convert to bam
    if not os.path.exists(bam_name):
        bam_cmd = 'samtools view -@ %s -bS %s -o %s' % (threads, input_file, bam_name)
        os.system(bam_cmd)
    
    # 2. sort the bam
    if not os.path.exists(bam_sort_name):
        sort_cmd = 'samtools sort -@ %s %s -o %s' % (threads, bam_name, bam_sort_name)
        os.system(sort_cmd)
    
    # 3. index the sorted bam
    if not os.path.exists('%s.sorted.bam.bai' % prefix):
        index_cmd = 'samtools index -@ %s %s' % (threads, bam_sort_name)
        os.system(index_cmd)
    


def callSNV(normal_bam, tumor_bam, ref, result_dir, caller='freebayes', threads=4):
    '''
    function to call SNV from sorted bam file
    normal_bam: normal sorted bam file for calling SNV
    tumor_bam: tumor sorted bam file to call SNV
    ref: reference genome
    result_dir: dir to save the results
    caller: SNV caller to use, freebayes or strelka
    '''
    checkpath(result_dir)
    if caller == 'freebayes':
        cmd = 'freebayes -f %s.fa %s > %s/freebayes.vcf' % (ref, tumor_bam, result_dir)
        os.system(cmd)
    elif caller == 'strelka':
        '''
        strelka_install_path/bin/configureStrelkaSomaticWorkflow.py \
        --normalBam HCC1187BL.bam \
        --tumorBam HCC1187C.bam \
        --referenceFasta hg19.fa \
        --runDir ${STRELKA_ANALYSIS_PATH}
        '''
        # 1. configuration for strelka caller
        cmd = '/home/assrivat/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py ' + \
            '--normalBam %s ' % normal_bam + \
            '--tumorBam %s ' % tumor_bam + \
            '--referenceFasta %s.fa ' % ref + \
            '--runDir %s'  % result_dir
        os.system(cmd)
        
        # 2. call strelka
        cmd = '%s/runWorkflow.py -m local -j %s' % (result_dir, threads)
        os.system(cmd)
    else:
        print('Only freebayes and strekla are available...')
        exit()
    
def callCNV(normal_bam, tumor_bam, ref, result_dir, caller='cnvkit', wgs=True, BED=None, threads=4):
    '''
    function to call CNV from sorted bam file
    normal_bam: normal sorted bam for calling CNV
    tumor_bam: timor sorted bam file to call CNV
    ref: reference genome
    result_dir: dir to save results
    caller: CNV caller to use, now we have cnvkit
    wgs: the Seq data is WGS or WES, if WES, BED file needed
    '''
    # 1. configuration for CNVkit
    if caller == 'cnvkit':
        prefix = tumor_bam.split('/')[-1].split('.')[0] + '.' + tumor_bam.split('/')[-1].split('.')[1]
        # work with WGS data, no BED file needed
        if wgs: 
            cmd1 = 'cnvkit.py batch %s ' % tumor_bam + \
                '--normal %s ' % normal_bam + \
                '--seq-method wgs ' + \
                '--segment-method cbs ' + \
                '--processes %s ' % threads + \
                '--output-dir %s ' % result_dir + \
                '--output-reference %s/%s_reference.cnn ' % (result_dir, prefix) + \
                '--fasta %s.fa ' % ref + \
                '--diagram --scatter |& tee -a %s/%s.log' % (result_dir, prefix)
            os.system(cmd1)
            # transfer to vcf file
            cmd2 = 'cnvkit.py export vcf %s/%s.call.cns -o %s/%s.cnv.vcf' % (result_dir, prefix, result_dir, prefix)
            os.system(cmd2)
        # work with WES data, BED file needed    
        else:
            if not BED:
                print('BED file needed for calling WES data...')
                exit()
            else:
                pass

def callSV(normal_bam, tumor_bam, ref, result_dir, caller='delly', threads=4):
    '''
    function to call SV from bam file
    '''
    prefix = tumor_bam.split('/')[-1].split('.')[0]
    if caller == 'delly':
        # call SV with delly
        #FIRST CHECK IF PAIRED, WHOLE GENOME 
        #cmd0 = 'EXPORT OMP_NUMTHREADS='
        cmd1 = '/home/assrivat/delly_v0.9.1_linux_x86_64bit call -o %s/%s.bcf ' % (result_dir, prefix) + \
            '-g %s.fa ' % ref + \
            '%s %s' % (tumor_bam, normal_bam)
        print(cmd1)
        # get vcf output
        cmd2 = '/home/assrivat/bcftool/bin/bcftools view %s/%s.bcf > %s/%s.vcf' % (result_dir, prefix, result_dir, prefix)
        os.system(cmd1)
        os.system(cmd2)
  


def align_normal():
    '''
    align normal fasta
    '''
    pass

def align_tumor():
    '''
    align tumor fasta
    '''
    pass

'''
results folder structure:
random 7 letters for parent folder (e.g kjwcrae)
kjwcrae
- results
  - normal
     normal.bam
     ...
  - tumor_0
    - sample_0:
        tumorB.bam
        - SNV
          freebayes.vcf
          strelka.vcf
    - sample_1
    ...
  - tumor_1
  ...
'''


def run_align_sort_index(data_name, ref, threads=4, normal=1, tumor_num=0, sample_num=0, single_cell_flag = False, single_cell_num = 0):

    '''
    go with semi-automatic first -- define which tumor and cell smaple for run
    '''
    data_dir = "/home/assrivat/simulation_results/%s" % data_name
    result_dir ="/home/assrivat/simulation_results/results/%s" % data_name
    checkpath(result_dir)
    if normal == 1:
        normal_dir = data_dir + '/reference'
        normal_res_dir = result_dir + '/normal'
        checkpath(normal_res_dir)
        # 0. define aligner
        info_dic = get_info(normal_dir + '/parameter_list.txt')
        read_len = int(info_dic['read len'])
        
        if read_len <= 500:
            aligner = 'bowtie2'
        else:
            aligner = 'minimap2'
        
        # 1. define variables
        input_file = glob.glob('%s/*.fasta' % normal_dir)
        output_file = normal_res_dir + '/normal.sam'
   
                
    else:
        info_directory = '{}/tumor_{}/samplenum_{}'.format(data_dir, tumor_num, sample_num)
        if not single_cell_flag: 
            tumor_dir = '%s/tumor_%s/samplenum_%s' % (data_dir, tumor_num, sample_num)
            tumor_res_dir = '%s/tumor_%s/samplenum_%s' % (result_dir, tumor_num, sample_num)
            checkpath(tumor_res_dir)
        else: 
            tumor_dir = '{}/tumor_{}/samplenum_{}_singlecell_{}'.format(data_dir, tumor_num, sample_num, single_cell_num)
            tumor_res_dir = '{}/tumor_{}/samplenum_{}_singlecell_{}'.format(result_dir, tumor_num, sample_num, single_cell_num)
            checkpath(tumor_res_dir)
        # 0. define the aligner
        info_dic = get_info(info_directory + '/parameter_list.txt')
        read_len = int(info_dic['read len'])
        
        if read_len <= 500:
            aligner = 'bowtie2'
        else:
            aligner = 'minimap2'
        
        # 1. define variables
        input_file = glob.glob('%s/*.fasta' % tumor_dir)
        print(input_file)
        output_file = tumor_res_dir + '/tumorB.sam'

    ref = ref
    alignment(input_file, output_file, ref, aligner=aligner, threads=threads)
    samtools_sort_index(output_file, threads=threads) 
        
    
    
def run_variant(data_name, ref, tumor_num, sample_num, \
                snv_caller="None", 
                cnv_caller="None", 
                sv_caller="None", 
                wgs=True,
                BED=None,
                threads=4, single_cell_flag = False, single_cell_num=0):
    '''
    call variants specific caller
    '''
    result_dir = "/home/assrivat/simulation_results/results/%s" % data_name
    # 0. specify the parameters
    normal_bam = '/home/assrivat/simulation_results/results/{}/normal/normal.sorted.bam'.format(data_name)
    if not single_cell_flag:
        tumor_bam = '%s/tumor_%s/samplenum_%s/tumorB.sorted.bam' % (result_dir, tumor_num, sample_num)
    else:
        tumor_bam = '{}/tumor_{}/samplenum_{}_singlecell_{}/tumorB.sorted.bam'.format(result_dir, tumor_num,sample_num, single_cell_num)
    ref = ref
    # 1. call snv
    if snv_caller != 'None':
        if not single_cell_flag:
            snv_result_dir = '%s/tumor_%s/samplenum_%s/snv' % (result_dir, tumor_num, sample_num)
        else: 
            snv_result_dir = '%s/tumor_%s/samplenum_%s_singlecell_%s/snv' % (result_dir, tumor_num, sample_num,single_cell_num)
        #pdb.set_trace()
        if snv_caller == 'freebayes':
            freebayes_dir = snv_result_dir + '/freebayes'
            checkpath(freebayes_dir)
            callSNV(normal_bam, tumor_bam, ref, freebayes_dir, threads=threads)

        elif snv_caller == 'strelka':
            strelka_dir = snv_result_dir + '/strelka'
            callSNV(normal_bam, tumor_bam, ref, strelka_dir, caller=snv_caller, threads=threads)

        else:
            print('Please choose freebayes or strelka...')
            exit()
    
    if cnv_caller != 'None':
        if not single_cell_flag:
            cnv_result_dir = '%s/tumor_%s/samplenum_%s/cnv' % (result_dir, tumor_num, sample_num)
        else: 
            cnv_result_dir = '%s/tumor_%s/samplenum_%s_singlecell_%s/cnv' % (result_dir, tumor_num, sample_num,single_cell_num) 
        if cnv_caller == 'cnvkit':
            cnvkit_dir = cnv_result_dir + '/cnvkit'
            checkpath(cnvkit_dir)
            callCNV(normal_bam, tumor_bam, ref, cnvkit_dir, wgs=wgs, BED=BED, threads=threads)
        else:
            print('Please choose cnvkit...')
            exit()
            
    if sv_caller != 'None':
        if not single_cell_flag:
            sv_result_dir = '%s/tumor_%s/samplenum_%s/sv' % (result_dir, tumor_num, sample_num)
        else: 
            sv_result_dir = '%s/tumor_%s/samplenum_%s_singlecell_%s/sv' % (result_dir, tumor_num, sample_num,single_cell_num) 
        checking_directory = '/home/assrivat/simulation_results/{}/tumor_{}/samplenum_{}'.format(data_name, tumor_num,sample_num)
        if(checkPaired(checking_directory)):
            if sv_caller == 'delly':
                delly_dir = sv_result_dir + '/delly'
                checkpath(delly_dir)
                callSV(normal_bam, tumor_bam, ref, delly_dir)
            else:
                print('Please choose delly...')
                exit()
    
def getTumorDirectories(data_name): 
    data_path = '/home/assrivat/simulation_results/{}/'.format(data_name)
    total_num_tumors = sum([os.path.isdir(data_path+i) for i in os.listdir(data_path)])-1
    list_of_samples = []
    for i in range(total_num_tumors): 
        current_tumor_path = data_path+'tumor_{}/*/'.format(i)
        subsamples = glob.glob(current_tumor_path)
        list_of_samples.append(subsamples)
    return list_of_samples

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
    # 1. data folder, also the parent folder to save results
    data_name = sys.argv[1]
    #parent_dir = "/home/assrivat/haoyun_files/%s" % data_name
    # 2. to align or not, 0 or 1
    align = int(sys.argv[2])
    # 3. if not normal sample, define the tumor sample: 0, 1, 2, ...
    tumor_num = sys.argv[3]
    # 4. under each tumor sample, define the cell sample: 0, 1, 2, ...
    sample_num = sys.argv[4]
    # 5. define the threads to use
    threads = sys.argv[5]
    # 6. define the caller
    snv_caller = sys.argv[6]
    cnv_caller = sys.argv[7]
    sv_caller = sys.argv[8]
    ref_name = sys.argv[9]
    ref = '/home/assrivat/simulation_results/ref/{}/{}'.format(ref_name,ref_name)
    samples = getTumorDirectories(data_name)
    # 1. run alignment
    if(align == 1):
        run_align_sort_index(data_name, ref, normal=1, tumor_num=tumor_num, sample_num=sample_num, threads=threads)
        print('sorted normal')
        for i in range(len(samples)): 
            for j in samples[i]:
                print('sorting:{},{}'.format(i,j))
                r1 = re.compile('samplenum_([0-9]*)')
                sample_num = r1.findall(j)[0]
                if 'singlecell' in j:
                    single_cell_flag = True
                    regex = re.compile('singlecell_([0-9]*)')
                    single_cell_num = regex.findall(j)[0]
                else: 
                    single_cell_flag = False
                    single_cell_num = 0
                run_align_sort_index(data_name, ref, normal = 0, tumor_num=str(i), sample_num=sample_num, threads = threads, single_cell_flag= single_cell_flag, single_cell_num= single_cell_num)
    # 2. run snv calling
    for i in range(len(samples)): 
        for j in samples[i]:
            print(j)
            print('sorting:{},{}'.format(i,j))
            r1 = re.compile('samplenum_([0-9]*)')
            sample_num = r1.findall(j)[0]
            print(sample_num)
            if 'singlecell' in j:
                single_cell_flag = True
                regex = re.compile('singlecell_([0-9]*)')
                single_cell_num = regex.findall(j)[0]
            else: 
                single_cell_flag = False
                single_cell_num = 0
            run_variant(data_name, ref, tumor_num=str(i), sample_num=sample_num, snv_caller=snv_caller, cnv_caller=cnv_caller, sv_caller=sv_caller, threads=threads,single_cell_flag = single_cell_flag,single_cell_num= single_cell_num)
if __name__ == '__main__':
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
