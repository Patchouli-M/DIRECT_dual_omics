from multiprocessing import Pool
import pandas as pd 
import sys
import os
import time

tool_lib=sys.path[0]

sample_lib="20230530_MT_6G"
reference_lib="reference/xx_reference/hg19"
reference_lamdba_lib='reference/lamdbaDNA/reference_bowtie1'

def get_sample(dir:str=sample_lib):
    data_dir=os.path.join( dir, "data" )
    print ("=========== Check Sample ===========")
    global sample_l
    sample_l=[]
    for root,dir,file in os.walk(data_dir):
        if (not file):
            continue
        file.sort()
        if "MD5.txt" in file :
            file.remove('MD5.txt')
        if (len(file)!=2):
            print (file)
            print ("Error sample, check please. ")
            sys.exit()
        for i in range(len(file)):
            if ( (file[i].split('_')[-1]) != str(i+1)+".fq.gz" ):
                print ("Error sample, check please. ")
                sys.exit()
        f1 = file[0].split('.fq.gz')[0][:-2]
        f2 = file[1].split('.fq.gz')[0][:-2]
        if (f1 != f2 ):
            print ("Error sample, check please. ")
            sys.exit()
        sample_l.append(f1)
        print (file)
    sample_l.sort()

def multi_sys_cmd(cmd:str):
    print(cmd)
    os.system(cmd)

def trim():
    # Input: Raw data
    # Output : trimmed data 
    print ("=========== Trim ===========")
    out_dir=os.path.join(sample_lib,"trimmed")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    cmd_l=[]
    for f in sample_l:
        f1 = f+"_1.fq.gz"
        f2 = f+"_2.fq.gz"
        f1 = os.path.join(sample_lib,"data",f,f1)
        f2 = os.path.join(sample_lib,"data",f,f2)
        cmd="trim_galore --quality 20 --phred33 --stringency 3 --gzip --length 36 --paired --cores 4 --trim1 --output_dir %s %s %s " % (out_dir,f1,f2)
        cmd_l.append(cmd)
    core_num = 16
    pool = Pool(core_num)
    pool.map(multi_sys_cmd, cmd_l)
    pool.close()
    print ("======== Trim Finish ========")

def lamdba_bismark():
    # Input: trimmed data 
    # Output : bismark bam
    print ("=========== lamdba bismark compare ===========")
    out_dir=os.path.join(sample_lib,"lamdba_bismark")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    cmd_l=[]
    for f in sample_l:
        f1 = f+"_1_val_1.fq.gz"
        f2 = f+"_2_val_2.fq.gz"
        f1 = os.path.join(sample_lib,"trimmed",f1)
        f2 = os.path.join(sample_lib,"trimmed",f2)
        cmd="bismark -fastq --output_dir  %s --temp_dir %s/tmp  --multicore 4 --non_directional -bowtie2 %s -1 %s -2 %s" % (out_dir,out_dir,reference_lamdba_lib,f1,f2)
        cmd_l.append(cmd)
    # core_num = 12
    # pool = Pool(core_num)
    # pool.map(multi_sys_cmd, cmd_l)
    # pool.close()
    for cmd in cmd_l:
        multi_sys_cmd(cmd)
    print ("======== lamdba bismark compare Finish ========")

def bismark():
    # Input: trimmed data 
    # Output : bismark bam
    print ("=========== bismark compare ===========")
    out_dir=os.path.join(sample_lib,"bismark")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    cmd_l=[]
    for f in sample_l:
        f1 = f+"_1_val_1.fq.gz"
        f2 = f+"_2_val_2.fq.gz"
        f1 = os.path.join(sample_lib,"trimmed",f1)
        f2 = os.path.join(sample_lib,"trimmed",f2)
        cmd="bismark -fastq --output_dir  %s --temp_dir %s/tmp --multicore 4 --non_directional -bowtie2 %s -1 %s -2 %s" % (out_dir,out_dir,reference_lib,f1,f2)
        cmd_l.append(cmd)
    # core_num = 12
    # pool = Pool(core_num)
    # pool.map(multi_sys_cmd, cmd_l)
    # pool.close()
    for cmd in cmd_l:
        multi_sys_cmd(cmd)
    print ("======== bismark compare Finish ========")

def bismark_dup():
    work_dir=os.path.join(sample_lib,"bismark")
    raw_dir=os.path.join(sample_lib,"bismark_raw")
    if not os.path.exists(raw_dir):
        os.makedirs(raw_dir)
    cmd = 'cp ' + work_dir + '/* ' + raw_dir
    print (cmd)
    os.system(cmd)
    print ("======== cp Raw Finish ========")
    for f in sample_l:
        f1 = f+"_1_val_1.fq.gz_bismark_bt2_pe.bam"
        f1 = os.path.join(work_dir,f1)
        cmd = 'deduplicate_bismark -p --bam '+f1
        print(cmd)
        os.system(cmd)
        d_f = f+"_1_val_1.fq.gz_bismark_bt2_pe.deduplicated.bam"
        d_f = os.path.join(work_dir,d_f)
        cmd = 'mv '+ d_f + ' ' +f1
        print (cmd)
        os.system(cmd)

def prepare_distin():
    out_dir=os.path.join(sample_lib,"dev_data")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for i in sample_l:
        sample_dir=os.path.join(sample_lib,"dev_data",i)
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
        bismark_f = i+"_1_val_1.fq.gz_bismark_bt2_pe.bam"
        bismark_dir=os.path.join(sample_lib,"bismark")

        f1 = i+"_1_val_1.fq.gz"
        f2 = i+"_2_val_2.fq.gz"
        f1 = os.path.join(sample_lib,"trimmed",f1)
        f2 = os.path.join(sample_lib,"trimmed",f2)

        print ("======== copy fq.gz "+i+" ========")
        cmd = "cp " + f1 + " " + os.path.join(sample_dir)
        multi_sys_cmd(cmd)
        cmd = "cp " + f2 + " " + os.path.join(sample_dir)
        multi_sys_cmd(cmd)
        print ("======== unzip "+i+" ========")
        cmd = "gunzip " +  os.path.join(sample_dir,i+"_1_val_1.fq.gz")
        multi_sys_cmd(cmd)
        cmd = "gunzip " +  os.path.join(sample_dir,i+"_2_val_2.fq.gz")
        multi_sys_cmd(cmd)
        print ("======== copy bam "+i+" ========")
        cmd = "cp " + os.path.join(bismark_dir,bismark_f) + " " + os.path.join(sample_dir,i+".bam")
        multi_sys_cmd(cmd)
        print ("======== samtools view -h "+i+" ========")
        cmd = "samtools view -h "+ os.path.join(sample_dir,i+".bam") + " > " + os.path.join(sample_dir,i+".sam")
        multi_sys_cmd(cmd)
        cmd = "rm " + os.path.join(sample_dir,i+".bam")
        multi_sys_cmd(cmd)

get_sample()
trim()
lamdba_bismark()
bismark() 
bismark_dup()
prepare_distin()