import sys
import os
import pickle  ## python 3.7
from scipy import stats
import math
#from itertools import islice

workdir = 'WHERE_GSE_FILES_PUT'
os.chdir(workdir)
os.getcwd()
STRING_HGNC_PATH = "9606.String_ID2HGNC.txt"
NETWORK_PATH = '9606.protein.links.v11.0.txt'
#SAMPLE_POOL_PATH = 'GSE155489_gv_pcos_counts.csv'##'data_nm_protein_log_FPKM.txt.human.5'
SAMPLE_POOL_PATH = 'GSE84958_untreat.txt'
FPKM_state = False
String_conf_level = 900

'''
data transform
reads to log2(FPKM+1)


'''
if(FPKM_state==False):
    fpool = open(SAMPLE_POOL_PATH)
    ## SAMPLE_POOL_PATH:
    ## GENE1	Expressions(numeric)	Expressions(numeric)......
    fo = open(SAMPLE_POOL_PATH+'_Transformed','w')

    for line in fpool:    
        seq=line.rstrip().split(',')    
        tmp=[]
        fo.write(seq[0])

        for one in seq[1:]: 
            one = '%.5f' % math.log(abs(float(one))+1,2)
            fo.write(','+str(one)) 
        fo.write('\n')

    fo.close()
    fpool.close()
	
class Person_PCC_Data:
## 构建PCC文件
    def __init__(self, string2HGNC , HGNC2string , EXP , POOL_LENGTH , PCC_POOL ):
        self.string2HGNC = string2HGNC
        self.HGNC2string = HGNC2string
        self.EXP = EXP
        self.POOL_LENGTH = POOL_LENGTH
        self.PCC_POOL = PCC_POOL
		
string2HGNC={}
HGNC2string={}

fa=open(STRING_HGNC_PATH)
## STRING_id	hgnc_id	hgnc_name
## 9606.ENSP00000228468	100	ACCN2
## 9606.ENSP00000228468	100	ASIC1

line_num = 0
for line in fa:
    seq=line.rstrip().split('\t') 
    string=seq[0]
    HGNC=seq[2]
    string2HGNC[string] = HGNC
    HGNC2string[HGNC] = string
    #print(HGNC)

fa.close()

## process STRING file

EDGE=set()
POINT=set()

fnet=open(NETWORK_PATH)
## NWTWORK_PATH FORMAT:
## protein1 protein2 combined_score
## 9606.ENSP00000000233 9606.ENSP00000272298 490

for line in fnet:
    seq=line.rstrip().split(' ')
    try:
        score=float(seq[2])
        if score >= String_conf_level:
            p1=seq[0]
            p2=seq[1]
            edge=[p1,p2]
            edge.sort()
            edge=':'.join(edge)
            EDGE.add(edge)
            POINT.add(p1)
            POINT.add(p2)
            #print(edge)
    except Exception as e:
        pass
print('At confidence score '+str(String_conf_level)+', there are '+str(len(EDGE))+' links.\n')
fnet.close()

## run batch
samples = ['GSM2254708',
'GSM2254709',
'GSM2254710',
'GSM2254711',
'GSM2254712',
'GSM2254713',
'GSM2254714',
'GSM2254715',
'GSM2254716',
'GSM2254717',
'GSM2254718',
'GSM2254719',
'GSM2254720',
'GSM2254721',
'GSM2254722',
'GSM2254723',
'GSM2254724',
'GSM2254725',
'GSM2254726',
'GSM2254727',
'GSM2254728',
'GSM2254729',
'GSM2254730',


'GSM2254731',
'GSM2254732',
'GSM2254733',
'GSM2254734',
'GSM2254735',
'GSM2254736',
'GSM2254737',
'GSM2254738',
'GSM2254739',
'GSM2254740',
'GSM2254741',
'GSM2254742',
'GSM2254743',
'GSM2254744',
'GSM2254745']

nsam = 37

print(samples[nsam])

EXP={}
Person_EXP={}
POOL_LENGTH=0

if(FPKM_state == False):
    fpool = open(SAMPLE_POOL_PATH+'_Transformed')
else:
    fpool = open(SAMPLE_POOL_PATH)
    
## SAMPLE_POOL_PATH:
## GENE1	Expressions(numeric)	Expressions(numeric)......
for line in fpool:
    
    seq=line.rstrip().split(',')    
    tmp=[]
    
    #for one in seq[7:]:  ## define background samples for building reference network
    for one in seq[1:24]: 
        tmp.append(float(one))
    POOL_LENGTH=len(tmp)
    EXP[seq[0]]=tmp    
    
    #Person_EXP[seq[0]]=float(seq[2]) ## define new individuals for Perturbed network
    Person_EXP[seq[0]]=float(seq[nsam+1])
    #print(EXP[seq[0]],", ",Person_EXP[seq[0]])

fpool.close()
print('Sample size is '+str(POOL_LENGTH)+'.\n')

from scipy import stats
PCC_POOL={}
for edge in EDGE:
    ps=edge.split(':')
    p1=ps[0]
    p2=ps[1]
    try:
        p1_HGNC = string2HGNC[p1]
        p2_HGNC = string2HGNC[p2]
        p1_exp = EXP[p1_HGNC]
        p2_exp = EXP[p2_HGNC]
        #p1_exp = [math.log([a+1,2]) for a in EXP[p1_HGNC]]        
        #p2_exp = [math.log([a+1,2]) for a in EXP[p2_HGNC]]  
        #print(p1_exp, '\n', p2_exp)
        pcc = float(stats.pearsonr(p1_exp, p2_exp)[0])
        #print(pcc)
        if math.isnan(pcc)== False:
            PCC_POOL[edge]= pcc
            #print(pcc)
    except Exception as e:
        pass


import sys
#fo = open(sys.argv[1],'wb')
fo = open('temp.pkl','wb')
data = Person_PCC_Data(string2HGNC , HGNC2string , EXP , POOL_LENGTH , PCC_POOL)
pickle.dump(data,fo)
print('There are '+str(len(PCC_POOL))+' relations(edges).\n')
fo.close()

## Raise PearsonRConstantInputWarning
## 原因：某序列没有变化 ，因此其标准偏差等于0，这将导致spearmanr()函数的除法为零，从而返回NaN 。

## 皮尔森相关系数的bug：
## 1 pearson系数不适用于小的或者非常稀疏的数据集，这样标准差会为0，不能作为分母。
## 2 只是描述了两组数据变化的趋势，对绝对值不敏感。

#print(dir(data))
#print(data.string2HGNC)
#print(type(data.string2HGNC))
#data.HGNC2string
#data.EXP
#data.POOL_LENGTH
#print(data.PCC_POOL)

print( '''
file1: PCCdata
file2: EXP, format: [col1: gene, col2: log2(FPKM+1)]
file3: Output

''')

'''
print('loading...')
#fdata = open(sys.argv[1])
fdata = open('temp.txt','rb')
data = pickle.load(fdata)
fdata.close()
print(data.string2HGNC)
print("loading done !")

'''


output=[]
#out_name = 'GSM4705188'
out_name = samples[nsam]
#fo = open(sys.argv[3],'wb')
fo = open(out_name,'w')
#fo1 = open(sys.argv[3]+'.sig','wb')
fo1 = open(out_name+'.sig','w')

for edge in data.PCC_POOL:
    ps = edge.split(':')
    p1 = ps[0]
    p2 = ps[1]
    #print(ps[0])
    
    try:
        p1_HGNC = data.string2HGNC[p1]
        #print(p1_HGNC)
        p2_HGNC = data.string2HGNC[p2]
        p1_exp_new = data.EXP[p1_HGNC]+[Person_EXP[p1_HGNC]]
        #print(p1_exp_new)
        p2_exp_new = data.EXP[p2_HGNC]+[Person_EXP[p2_HGNC]]
        pcc_new = float(stats.pearsonr(p1_exp_new,p2_exp_new)[0])
        #print(pcc_new)
        pcc = data.PCC_POOL[edge]
        delta_pcc = pcc_new - pcc
        #print(delta_pcc)
        z = delta_pcc /( (1-pcc**2)/(data.POOL_LENGTH-1) )
        p = float(stats.norm.sf(abs(z))*2)
        #print(p)
        line = p1_HGNC +'\t'+p1+'\t'+ p2_HGNC +'\t'+ p2 +'\t'+ str(pcc)+'\t'+str(pcc_new)+'\t'+str(z)+'\t'+str(p)
        if '\tnan' not in line and '\t-inf\t' not in line and '\tinf\t' not in line:
           # print(line)
            output.append([p,-abs(z),p1,p2,line])
    except Exception as e:
        pass #print e

output.sort()



## 挑出最显著的关系

ALL={}
test_time=len(output)
for one in output:
    #print(one)
    pvad = one[0] * test_time
    #print(pvad, type(pvad))
    fo.write(one[-1]+'\t'+ str(pvad) +'\n'  )
    #print(str(pvad) +'\n'  )
    if pvad <0.05:
        fo1.write(one[-1]+'\t'+ str(pvad) +'\n' )
    p1_= one[2]
    p2 = one[3]
    if p1 in ALL:
            ALL[p1].append(pvad)
    else:
            ALL[p1]=[pvad]
    if p2 in ALL:
            ALL[p2].append(pvad)
    else:
            ALL[p2]=[pvad]

fo.close()
fo1.close()



## 挑出最显著的gene

output=[]
for p in ALL:
    num_sigpv=0
    for one in ALL[p]:
        if one < 0.05:
            num_sigpv += 1
    output.append([-num_sigpv,p])
output.sort()
#print(output)


#fo2=open(sys.argv[3]+'.gene','w')
#fo3=open(sys.argv[3]+'.gene.sig','w')
fo2=open(out_name+'.gene','w')
fo3=open(out_name+'.gene.sig','w')
count = 0

for one in output:
    p=one[1]
    human=data.string2HGNC[p]
    #min_pv=str(one[0])
    num_sigpv=str(-one[0])
    
    fo2.write(p+'\t'+human+'\t'+num_sigpv+'\n')
    if float(num_sigpv)>0:
        fo3.write(p+'\t'+human+'\t'+num_sigpv+'\n')
        count+=1
        
print('There are '+str(count)+' significant genes.\n')
fo2.close()
fo3.close()

#compared z-value =  $2 - $1

print('''
file1: step1 output 1
file2: step1 output 2
output: 2-1

''')

PCC_origin={}
PCC1={}
PCC2={}

PAIR1={}
f1 = open('GSM4705187.sig')
#f1=open(sys.argv[1])
 
for line in f1:
    seq=line.rstrip().split('\t')
    tag=':'.join(seq[0:4])
    #print seq
#    try:
    z=float(seq[6])
    PCC1[tag]=seq[5]
    PCC_origin[tag]=seq[4]
 #   except Exception as e:
  #      print seq
    PAIR1[tag]=z
f1.close()

PAIR2={}
    
f2 = open('GSM4705188.sig')
#f2=open(sys.argv[2])

for line in f2:
    seq=line.rstrip().split('\t')
    tag=':'.join(seq[0:4])
    z=float(seq[6])
    PAIR2[tag]=z
    PCC2[tag]=seq[5]
f2.close()


output=[]
for one in PAIR1:
    if one in PAIR2:
        z=PAIR2[one]-PAIR1[one]
        pc0=PCC_origin[one]
        pc1=PCC1[one]
        pc2=PCC2[one]
        
        tag=one.replace(':', '\t')
        output.append([abs(z),z,tag,pc0,pc1,pc2])


output.sort(reverse=True) 

#fo=open(sys.argv[3],'w')
fo=open("compared",'w')

for one in output:
    #if one[0] > 1.97:
        fo.write(one[2]+'\t'+str(one[1])+'\t'+str(one[3])+'\t'+str(one[4])+'\t'+str(one[5])+'\n')

fo.close()
    
print('''
get gene
''')

#fi=open(sys.argv[1])
fi=open("compared")
fo=open("compared"+'.gene','w')

PAIR={}

for line in fi:
    seq=line.split('\t')
    p1=seq[0]
    p2=seq[2]
    tag=[p1,p2]
    tag.sort()
    tag=':'.join(tag)
    PAIR[tag]=float(seq[4])
fi.close()
    
GENE={}

for pair in PAIR:
    p1=pair.split(':')[0]
    p2=pair.split(':')[1]
    z=abs(PAIR[pair])
    if p1 in GENE:
        GENE[p1] += z
    else:
        GENE[p1]=z
    if p2 in GENE:
        GENE[p2] += z
    else:
        GENE[p2]=z

output=[]
for one in GENE:
    output.append([GENE[one],one])
i=1

output.sort(reverse=True)
for one in output:
    fo.write(str(i)+'\t'+one[1]+'\t'+str(one[0])+'\n')
    i+=1
fo.close()

