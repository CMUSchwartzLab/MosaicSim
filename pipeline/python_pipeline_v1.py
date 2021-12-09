#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 09:54:33 2021

@author: leo
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
    if os.path.exists(bam_name): # check bam index file for less strict checking
        print("Aligned File already exists, quit...")
        exit()
        
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
        cmd = '/home/haoyunl/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py ' + \
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
        cmd1 = '/home/haoyunl/software/delly/delly_v0.8.7_linux_x86_64bit call -o %s/%s.bcf ' % (result_dir, prefix) + \
            '-g %s.fa ' % ref + \
            '%s %s' % (tumor_bam, normal_bam)
        # get vcf output
        cmd2 = 'bcftools view %s/%s.bcf > %s/%s.vcf' % (result_dir, prefix, result_dir, prefix)
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


def run_align_sort_index(data_name, ref, threads=4, normal=1, tumor_num=0, sample_num=0):

    '''
    go with semi-automatic first -- define which tumor and cell smaple for run
    '''
    data_dir = "/home/assrivat/haoyun_files/%s" % data_name
    result_dir ="/home/assrivat/haoyun_files/results/%s" % data_name
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
        tumor_dir = '%s/tumor_%s/samplenum_%s' % (data_dir, tumor_num, sample_num)
        tumor_res_dir = '%s/tumor_%s/samplenum_%s' % (result_dir, tumor_num, sample_num)
        checkpath(tumor_res_dir)
        
        # 0. define the aligner
        info_dic = get_info(tumor_dir + '/parameter_list.txt')
        read_len = int(info_dic['read len'])
        
        if read_len <= 500:
            aligner = 'bowtie2'
        else:
            aligner = 'minimap2'
        
        # 1. define variables
        input_file = glob.glob('%s/*.fasta' % tumor_dir)
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
                threads=4):
    '''
    call variants specific caller
    '''
    result_dir = "/home/assrivat/haoyun_files/results/%s" % data_name
    # 0. specify the parameters
    normal_bam = '/home/assrivat/haoyun_files/results/normal/normal.sorted.bam'
    tumor_bam = '%s/tumor_%s/samplenum_%s/tumorB.sorted.bam' % (result_dir, tumor_num, sample_num)
    ref = ref
    # 1. call snv
    if snv_caller != 'None':
        snv_result_dir = '%s/tumor_%s/samplenum_%s/snv' % (result_dir, tumor_num, sample_num)
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
        cnv_result_dir = '%s/tumor_%s/samplenum_%s/cnv' % (result_dir, tumor_num, sample_num)
        if cnv_caller == 'cnvkit':
            cnvkit_dir = cnv_result_dir + '/cnvkit'
            checkpath(cnvkit_dir)
            callCNV(normal_bam, tumor_bam, ref, cnvkit_dir, wgs=wgs, BED=BED, threads=threads)
        else:
            print('Please choose cnvkit...')
            exit()
            
    if sv_caller != 'None':
        sv_result_dir = '%s/tumor_%s/samplenum_%s/sv' % (result_dir, tumor_num, sample_num)
        if sv_caller == 'delly':
            delly_dir = sv_result_dir + '/delly'
            checkpath(delly_dir)
            callSV(normal_bam, tumor_bam, ref, delly_dir)
        else:
            print('Please choose delly...')
            exit()
    
    

def main():
    ref = '/home/assrivat/haoyun_files/ref/hg38/hg38'
    # 1. data folder, also the parent folder to save results
    data_name = sys.argv[1]
    #parent_dir = "/home/assrivat/haoyun_files/%s" % data_name
    # 2. if normal or not, 0 or 1
    normal = sys.argv[2]
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
    
    # 1. run alignment
    # run_align_sort_index(data_name, ref, normal=normal, \
    #     tumor_num=tumor_num, sample_num=sample_num, threads=threads)
    # 2. run snv calling
    run_variant(data_name, ref, tumor_num, sample_num, 
                snv_caller=snv_caller, 
                cnv_caller=cnv_caller,
                sv_caller=sv_caller,
                threads=threads)
    

if __name__ == '__main__':
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
