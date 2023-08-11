
from itertools import count
from multiprocessing import Pool,Process
from random import sample
import pandas as pd 
import sys
import os
import time
import subprocess
import fcntl
import math

def init(sample:str):
    global sam_file
    global raw_fq
    global tr_fq
    global nu_fq
    global temp_dir
    sam_file = os.path.join(dev_dir,sample,sample+'.sam')
    raw_fq = os.path.join(dev_dir,sample,sample)
    temp_dir = os.path.join(dev_dir,sample,'tmp')
    tr_fq = os.path.join(temp_dir,sample+'_tr_out')
    nu_fq = os.path.join(temp_dir,sample+'_nu_out')
    '''
    cmd = "cat /dev/null > "+tr_fq+"_1.fq"
    os.system(cmd)
    cmd = "cat /dev/null > "+tr_fq+"_2.fq"
    os.system(cmd)
    cmd = "cat /dev/null > "+nu_fq+"_1.fq"
    os.system(cmd)
    cmd = "cat /dev/null > "+nu_fq+"_2.fq"
    os.system(cmd)
    '''

def grep_f(seq:str,file:str,rate:float,spli_num:int):

    file1=file+"_1_val_1.fq"
    file2=file+"_2_val_2.fq"
    seq1=seq
    seq2=seq.replace('1:N:0','2:N:0')
    

    try:
        out1 = fq1.index('@'+seq1+'\n')
        out2 = fq2.index('@'+seq2+'\n')
    except:
        return 
    
    f_nu_io1 = open(nu_fq+'_1.fq_'+str(spli_num),'a')
    f_nu_io2 = open(nu_fq+'_2.fq_'+str(spli_num),'a')
    f_tr_io1 = open(tr_fq+'_1.fq_'+str(spli_num),'a')
    f_tr_io2 = open(tr_fq+'_2.fq_'+str(spli_num),'a')


    if (rate >= 0.9):
        for i in range(4):
            print (fq1[out1+i],file=f_tr_io1,end='')
            print (fq2[out2+i],file=f_tr_io2,end='')

    if (rate <= 0.5):
        for i in range(4):
            print (fq1[out1+i],file=f_nu_io1,end='')
            print (fq2[out2+i],file=f_nu_io2,end='')
    f_nu_io1.close()
    f_nu_io2.close()
    f_tr_io1.close()
    f_tr_io2.close()


def split_raw_sam(sample):
    temp_dir = os.path.join(dev_dir,sample,'tmp')
    if not os.path.exists( temp_dir ):
        os.makedirs(temp_dir)
    else :
        for f in os.listdir(temp_dir):
            os.remove(os.path.join(temp_dir,f))
    lines=[]
    last=''
    with open(sam_file, 'r') as f:
        for line in f:
            temp_seq = line.split('\t')[0]
            meth_state=line.split('\t')[-3]
            if (meth_state.startswith("XM:")):
                '''
                # if (temp_seq!=last):
                    # meth_state=line.split('\t')[-3]
                    # lines.append(line)
                # last = temp_seq
                '''
                lines.append(line)
    print (len(lines))
    count = 0
    spli_num = 1
    len_per_sam = math.ceil((len(lines)/2)/split_k)*2
    fw = open(os.path.join(temp_dir,sample+".sam_"+str(spli_num)),'w')
    for line in lines :
        count += 1
        fw.write(line)
        if (count>=len_per_sam):
            count = 0 
            spli_num += 1
            fw = open(os.path.join(temp_dir,sample+".sam_"+str(spli_num)),'w')

def calc_per_split_num(sam_file_split):

    spli_num = sam_file_split.split('_')[-1]
    print ("========= Begin "+str(spli_num)+" ================")
    os.system("date")
    os.system("wc -l "+sam_file_split)
    methy_dict={}
    with open(sam_file_split, 'r') as f:
        for line in f:
            meth_state=line.split('\t')[-3]
            if (not meth_state.startswith("XM:Z:")):
                return
            meth_state=meth_state.split("XM:Z:")[-1]
            seq=line.split('\t')[0]
            seq=seq.replace('_',' ')
            count_x=meth_state.count('x')
            count_z=meth_state.count('z')
            count_h=meth_state.count('h')
            count_u=meth_state.count('u')
            count_X=meth_state.count('X')
            count_Z=meth_state.count('Z')
            count_H=meth_state.count('H')
            count_U=meth_state.count('U')

            non_CG_count_meth = count_X+count_U+count_H
            non_CG_count_all=non_CG_count_meth+count_x+count_u+count_h

            if (non_CG_count_all==0):
                continue
            rate = non_CG_count_meth/non_CG_count_all
            if (methy_dict.__contains__(seq)):
                methy_dict[seq].append(rate)
            else :
                methy_dict[seq] = []
                methy_dict[seq].append(rate)
    
        for seq in methy_dict :
            avr_methy_rate = sum(methy_dict[seq])/len(methy_dict[seq])
            grep_f (seq,raw_fq,avr_methy_rate,spli_num)
        print ("========= Finish "+str(spli_num)+" ================")
        os.system("date")
        os.system("wc -l "+sam_file_split)
# ################
def merge(sample):
    temp_dir = os.path.join(dev_dir,sample,'tmp')
    tr_fq = os.path.join(dev_dir,sample,sample+'_tr_out')
    nu_fq = os.path.join(dev_dir,sample,sample+'_nu_out')
    w_nu_io1 = open(nu_fq+'_1.fq','w')
    w_nu_io2 = open(nu_fq+'_2.fq','w')
    w_tr_io1 = open(tr_fq+'_1.fq','w')
    w_tr_io2 = open(tr_fq+'_2.fq','w')
    for spli_num in range(1,split_k+1,1):
        tmp_tr_fq = os.path.join(temp_dir,sample+'_tr_out')
        tmp_nu_fq = os.path.join(temp_dir,sample+'_nu_out')
        read_nu_io1 = open(tmp_nu_fq+'_1.fq_'+str(spli_num),'r')
        read_nu_io2 = open(tmp_nu_fq+'_2.fq_'+str(spli_num),'r')
        read_tr_io1 = open(tmp_tr_fq+'_1.fq_'+str(spli_num),'r')
        read_tr_io2 = open(tmp_tr_fq+'_2.fq_'+str(spli_num),'r')
        for line in read_nu_io1:
            w_nu_io1.write(line)
        for line in read_nu_io2:
            w_nu_io2.write(line)
        for line in read_tr_io1:
            w_tr_io1.write(line)
        for line in read_tr_io2:
            w_tr_io2.write(line)
    nuf1=nu_fq+'_1.fq'
    nuf2=nu_fq+'_2.fq'
    trf1=tr_fq+'_1.fq'
    trf2=tr_fq+'_2.fq'
    print ("nu")
    cmd = "echo `wc -l "+nuf1+"|awk '{print $1}' `\"/4\" |bc"
    os.system(cmd)
    cmd = "echo `wc -l "+nuf2+"|awk '{print $1}' `\"/4\" |bc"
    os.system(cmd)
    print ("tr")
    cmd = "echo `wc -l "+trf1+"|awk '{print $1}' `\"/4\" |bc"
    os.system(cmd)
    cmd = "echo `wc -l "+trf2+"|awk '{print $1}' `\"/4\" |bc"
    os.system(cmd)
    os.system("date")
    print (sample,"============== All Finish! ================")
# ################

def run(sample:str):
    os.system("date")
    init(sample)
    lines=[]
    last=''
    # 分治
    split_raw_sam(sample)

    file1=raw_fq+"_1_val_1.fq"
    file2=raw_fq+"_2_val_2.fq"
    global fq1
    fq1=[]
    with open(file1, 'r') as f:
        for line in f:
            fq1.append(line)
    global fq2
    fq2=[]
    with open(file2, 'r') as f:
        for line in f:
            fq2.append(line)
    if (len(fq1) !=len(fq2)):
        print ("Error",sample)
        return
    
    # 对于每个文件
    sam_file_split_l = []
    for split_num in range(1,split_k+1,1):
    
        sam_file_split = os.path.join(temp_dir,sample+".sam_"+str(split_num))
        sam_file_split_l.append(sam_file_split)    

    pool=Pool(core_num)
    pool.map(calc_per_split_num,sam_file_split_l)
    pool.close()
    
    print ("==============================")
    print(sample)
    merge(sample)



core_num = 256
split_k = core_num



dev_dir = "20230530_MT_6G/dev_data2"

sample_l = ['MT-54_L4']

my_process=[]
for i in sample_l:
    print (i)
    my_process.append(Process(target=run,args=(i,)))
for i in range(len(my_process)):
    my_process[i].start()
