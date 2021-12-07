import ast
import math
f = open('/home/assrivat/haoyun_files/atjiiuk/tumor_0/information_list.txt', 'r')
d = f.readlines()[0]

actual_d = ast.literal_eval(d)
all_actual_muts = set()
counter = 0
for i in actual_d:
        for k in actual_d[i]:
		#snv check
                alternate_tuple = (k[0], int(k[2]/2))
                counter += 1
                all_actual_muts.add(tuple(k))
                all_actual_muts.add(tuple(alternate_tuple))
#print(all_actual_muts)
print(len(all_actual_muts))
f2 = open('/home/assrivat/somatic.snvs.vcf')
d = f2.readlines()
print('loaded')
chrom_dict = {'chr1':0,'chr10':1, 'chr11':2, 'chr12':3, 'chr13':4,'chr14':5, 'chr15':6 ,'chr16':7, 'chr17':8, 'chr18':9, 'chr19':10,'chr2':11, 'chr20':12, 'chr21':13, 'chr22':14, 'chr3':15, 'chr4':16, 'chr5':17, 'chr6':18, 'chr7':19, 'chr8':20,'chr9':21, 'chrX':22, 'chrY':23}
clist = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13','chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19','chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8','chr9', 'chrX', 'chrY']
vcf_lines = 0
match = 0
mismatch = 0
hitset = set()
missset = set()
for i in d: 
        z = i.split('\t')
        #EDIT THIS LINE TO FIT HTE FILE
        if(len(z) == 11 and z[0] in clist):
                chrom = chrom_dict[z[0]]
                pos = int(z[1])-1
                a_tuple = (pos, chrom)
                if(a_tuple in all_actual_muts): 
                        hitset.add(a_tuple)
                        match += 1
                else: 
                        missset.add(a_tuple)
                        mismatch +=1
                vcf_lines += 1
print(vcf_lines)
print(match)
print(mismatch)
#print(hitset)
#print(missset)
for i in missset:
        mindist = 1e8
        for j in all_actual_muts:
                if (i[1] == j[1]):
                        dist = abs(i[0] - j[0])
                        if(dist < mindist):
                                mindist = dist
        #print(mindist)

