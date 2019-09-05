#title         : get_rMATS_res.py
#description      : get rMATS(4.0.2) result 
#author         : Xin Dong
#version        :0.1
#notes         : only worked in rMATS 4.0.2
#python_version     :3.6


import pandas as pd
import os,sys
import argparse

def arg():
    arg = argparse.ArgumentParser()
    arg.add_argument('--path',type=str,help='path to rMATS reslut dir')
    arg.add_argument('--b1',type=str,help='b1')
    arg.add_argument('--b2',type=str,help='b2')
    arg.add_argument('--out',type=str,help='b2')
    parse = arg.parse_args()
    return parse

def get_A3SS(file,A3SS):
    data = pd.read_csv(file,index_col=0,sep='\t')
    A3SS_data = pd.read_csv(A3SS,index_col=0,sep='\t')
    def get_pos(x):
        if x[3] == '-':
            return '\t'.join([x[1], x[2],x[3],str(x[7]),str(x[8])])
        else:
            return '\t'.join([x[1], x[2],x[3],str(x[9]),str(x[6])])
    A3SS_pos = data.apply(get_pos,axis=1)
    res = pd.concat([A3SS_data,A3SS_pos],axis=1).loc[:,['SJC_SAMPLE_1','SJC_SAMPLE_2',0]]
    return res

def get_A5SS(file,A5SS):
    data = pd.read_csv(file,index_col=0,sep='\t')
    A5SS_data = pd.read_csv(A5SS,index_col=0,sep='\t')
    def get_pos(x):
        if x[3] == '+':
            return '\t'.join([x[1], x[2],x[3],str(x[7]),str(x[8])])
        else:
            return '\t'.join([x[1], x[2],x[3],str(x[9]),str(x[6])])
    A5SS_pos = data.apply(get_pos,axis=1)
    res = pd.concat([A5SS_data,A5SS_pos],axis=1).loc[:,['SJC_SAMPLE_1','SJC_SAMPLE_2',0]]
    return res

def get_SE(file,SE):
    data = pd.read_csv(file,index_col=0,sep='\t')
    SE_data = pd.read_csv(SE,index_col=0,sep='\t')
    def get_pos(x):
        return '\t'.join([x[1], x[2],x[3],str(x[7]),str(x[8])])
    SE_pos = data.apply(get_pos,axis=1)
    res = pd.concat([SE_data,SE_pos],axis=1).loc[:,['SJC_SAMPLE_1','SJC_SAMPLE_2',0]]
    return res

def get_RI(file,RI):
    data = pd.read_csv(file,index_col=0,sep='\t')
    RI_data = pd.read_csv(RI,index_col=0,sep='\t')
    def get_pos(x):
        return '\t'.join([x[1], x[2],x[3],str(x[7]),str(x[8])])
    RI_pos = data.apply(get_pos,axis=1)
    res = pd.concat([RI_data,RI_pos],axis=1).loc[:,['IJC_SAMPLE_1','IJC_SAMPLE_2',0]]
    return res

def get_MXE(file,MXE):
    data = pd.read_csv(file,index_col=0,sep='\t')
    MXE_data = pd.read_csv(MXE,index_col=0,sep='\t')
    def get_pos(x):
        return '\t'.join([x[1], x[2],x[3],str(x[9]),str(x[10])])
    MXE_pos = data.apply(get_pos,axis=1)
    res = pd.concat([MXE_data,MXE_pos],axis=1).loc[:,['SJC_SAMPLE_1','SJC_SAMPLE_2',0]]
    return res

def main(arg):
    path, b1, b2 = arg.path, arg.b1, arg.b2
    with open(b1) as f1,open(b2) as f2:
        l1, l2 = f1.readline(),f2.readline()
    A3SS_files = ('%sA3SS.MATS.JC.txt'%path,'%sfromGTF.A3SS.txt'%path)
    A5SS_files = ('%sA5SS.MATS.JC.txt'%path,'%sfromGTF.A5SS.txt'%path)
    SE_files = ('%sSE.MATS.JC.txt'%path,'%sfromGTF.SE.txt'%path)
    RI_files = ('%sRI.MATS.JC.txt'%path,'%sfromGTF.RI.txt'%path)
    MXE_files = ('%sMXE.MATS.JC.txt'%path,'%sfromGTF.MXE.txt'%path)
    A3SS, A5SS, SE, RI, MXE = get_A3SS(A3SS_files[1],A3SS_files[0]), get_A5SS(A5SS_files[1],A5SS_files[0]),get_SE(SE_files[1],SE_files[0]),get_RI(RI_files[1],RI_files[0]), get_MXE(MXE_files[1],MXE_files[0])
    A3SS['AS_type'] = 'A3SS'
    A5SS['AS_type'] = 'A5SS'
    SE['AS_type'] = 'SE'
    RI['AS_type'] = 'RI'
    MXE['AS_type'] = 'MXE'
    A3SS.columns = ['sample1','sample2','id','AS_type']
    A5SS.columns = ['sample1','sample2','id','AS_type']
    SE.columns = ['sample1','sample2','id','AS_type']
    RI.columns = ['sample1','sample2','id','AS_type']
    MXE.columns = ['sample1','sample2','id','AS_type']
    total = pd.concat([A3SS, A5SS, SE, RI, MXE],axis=0)
    total = total.loc[ (total['sample1'].notna() &  total['sample2'].notna()) ,:]
    total['sample1'] = total['sample1'].astype(str)
    total['sample2'] = total['sample2'].astype(str)
    total = pd.concat([total['sample1'].str.split(',',expand=True),total['sample2'].str.split(',',expand=True),total['id'],total['AS_type']],axis=1)
    header = l1.strip().split(',')+l2.strip().split(',')+['gene\tchr\tstrand\tstart\tend']+['AS_type']
    total.columns = header
    return total

def get_sample(x):
    return ','.join(x.loc[x!=0].index)

arg = arg()
total = main(arg)
total = total.melt(id_vars=['gene\tchr\tstrand\tstart\tend','AS_type'])
total['value'] = total['value'].astype(int)
total = total.loc[total['value'] > 0,:]
total.to_csv(arg.out,sep='\t',index=False)
