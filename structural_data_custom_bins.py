#!/usr/bin/env python3
import numpy as np
import os,sys,argparse
import pandas as pd
import time
import glob
import matplotlib as plt
import multiprocessing as mp
from multiprocessing import Pool
import glob
import pandas as pd
from functools import partial
from multiprocessing import Process

def structure(genes, structure, transcript, utr5, cds, utr3, output):
    # list of genes to check
    first = pd.read_csv(genes, delim_whitespace=True)
    # Structural data
    DMS_merged = pd.read_csv(structure, delim_whitespace=True, names=['chr', 'start_map', 'end_map', 'n', 'score', 'strand', 'chr_l', 'exon_start', 'exon_end', 'strand1', 'gene_name', 'transcript_id', 'utr5_length', 'cds_length', 'utr3_length', 'exon_length', 'exon_number'], usecols=['chr', 'start_map', 'end_map','score', 'strand', 'exon_start', 'exon_end', 'gene_name', 'transcript_id', 'utr5_length', 'cds_length', 'utr3_length', 'exon_length', 'exon_number'])
    first_DMS = first.merge(DMS_merged, on=['transcript_id'])
    first_structure = first_DMS[['transcript_id']]
    first_structure = first_structure.drop_duplicates().reset_index()
    del first_structure['index']
    first_structure.to_csv(transcript, index=False, sep='\t')
    
    # 5'UTR preparation
    utr5 = pd.read_csv(utr5, delim_whitespace=True, names=['chr', 'start', 'end', 'transcript_id', 'z', 'strand', 'start1', 'end1', 'n', 'number_of_exons', 'length', 'cordinates'], usecols=['chr', 'start', 'end', 'transcript_id', 'strand', 'number_of_exons', 'length', 'cordinates'])
    utr5['transcript_id'] = utr5['transcript_id'].str.rstrip('_5utr')
    utr5['utr5_length'] = utr5['length'].apply(lambda x: sum(map(int, x.split(',')[:-1])))
    utr5_pos = utr5[utr5['strand'] == '+']
    utr5_neg = utr5[utr5['strand'] == '-']
    
    #Positive data utr5
    results_pos = pd.DataFrame(columns=['transcript_id', 'chr', 'strand', 'start', 'end', 'exons_start', 'exons_end'])
    for index, row in utr5_pos.iterrows():
        coord = row['cordinates'].split(',')
        coord = coord[:-1]
        length = row['length'].split(',')
        length = length[:-1]
        start = row['start']
        end = row['end']
        tr = row['transcript_id']
        ch = row['chr']
        st = row['strand']
        x = []
        y = []
        for i in range(len(length)):
            x.append(int(coord[i]) + int(start))
            y.append(int(coord[i]) + int(start) + int(length[i]))
        results_pos.loc[index] = [tr, ch, st, start, end, x, y] 
    test_first = first_DMS.merge(results_pos, on=['transcript_id', 'chr', 'strand'])
    zobry_utr5 = pd.DataFrame()
    for i in test_first['transcript_id'].unique():
        tmp = test_first[test_first['transcript_id'] == i]
        lol = []
        for m in range(len(tmp.iloc[0]['exons_start'])):
            for j in range(tmp.iloc[0]['exons_start'][m], (tmp.iloc[0]['exons_end'][m])):
                x = pd.DataFrame()
                x = tmp[tmp['start_map'] == j]
                if len(x) > 0:
                    lol.append(x['score'].values[0])
                else:
                    lol.append(0)
        df = pd.DataFrame()
        df[i] = lol
        df1 = df.T
        for n in range(1, 21):
            if n == 1:
                df1[str(n) + 'A'] = (len(df)/20) * n
                df1[str(n) + '_final'] = df1.iloc[:, :round(df1[str(n) + 'A'].values[0])].values.mean()
            else:
                df1[str(n) + 'A'] = (len(df)/20) * n
                df1[str(n) + '_final'] = df1.iloc[:, round(df1[str(n-1) + 'A'].values[0]):round(df1[str(n) + 'A'].values[0])].values.mean()
        names1 = list(range(1, 21))
        mystring = '_final'
        names = ["{}{}".format(l,mystring) for l in names1]
        a7a = pd.DataFrame(df1, columns=names)
        zobry_utr5 = zobry_utr5.append(a7a)   
        
    #Negative data utr5
    results_neg = pd.DataFrame(columns=['transcript_id', 'chr', 'strand', 'start', 'end', 'exons_start', 'exons_end'])
    for index, row in utr5_neg.iterrows():
        coord = row['cordinates'].split(',')
        coord = coord[:-1]
        length = row['length'].split(',')
        length = length[:-1]
        start = row['start']
        end = row['end']
        tr = row['transcript_id']
        ch = row['chr']
        st = row['strand']
        x = []
        y = []
        for i in range(len(length)):
            x.append(int(coord[i]) + int(start))
            y.append(int(coord[i]) + int(start) + int(length[i]))
        results_neg.loc[index] = [tr, ch, st, start, end, x, y] 
        
    test_second = first_DMS.merge(results_neg, on=['transcript_id', 'chr', 'strand'])
    zobry_utr5_neg = pd.DataFrame()
    for i in test_second['transcript_id'].unique():
        tmp = test_second[test_second['transcript_id'] == i]
        lol = []
        for m in range(len(tmp.iloc[0]['exons_start'])):
            for j in range(tmp.iloc[0]['exons_start'][m], (tmp.iloc[0]['exons_end'][m])):
                x = pd.DataFrame()
                x = tmp[tmp['start_map'] == j]
                if len(x) > 0:
                    lol.append(x['score'].values[0])
                else:
                    lol.append(0)
        df = pd.DataFrame()
        df[i] = lol
        df1 = df.T
        for n in range(1, 21):
            x = 21 - n
            if n == 1:
                df1[str(n) + 'A'] = (len(df)/20) * n
                df1[str(x) + '_final'] = df1.iloc[:, :round(df1[str(n) + 'A'].values[0])].values.mean()
            else:
                df1[str(n) + 'A'] = (len(df)/20) * n
                df1[str(x) + '_final'] = df1.iloc[:, round(df1[str(n-1) + 'A'].values[0]):round(df1[str(n) + 'A'].values[0])].values.mean()
        names1 = list(range(1, 21))
        mystring = '_final'
        names = ["{}{}".format(l,mystring) for l in names1]
        a7a = pd.DataFrame(df1, columns=names)
        zobry_utr5_neg = zobry_utr5_neg.append(a7a)    
    frames = [zobry_utr5, zobry_utr5_neg]
    utr5 = pd.concat(frames)
    
    # CDS preparation
    CDS = pd.read_csv(cds, delim_whitespace=True, names=['chr', 'start', 'end', 'transcript_id', 'z', 'strand', 'start1', 'end1', 'n', 'number_of_exons', 'length', 'cordinates'], usecols=['chr', 'start', 'end', 'transcript_id', 'strand', 'number_of_exons', 'length', 'cordinates'])
    CDS['transcript_id'] = CDS['transcript_id'].str.rstrip('_cds')
    CDS['CDS_length'] = CDS['length'].apply(lambda x: sum(map(int, x.split(',')[:-1])))
    CDS_pos = CDS[CDS['strand'] == '+']
    CDS_neg = CDS[CDS['strand'] == '-']
    
    #Positive data CDS
    results_pos = pd.DataFrame(columns=['transcript_id', 'chr', 'strand', 'start', 'end', 'exons_start', 'exons_end'])
    for index, row in CDS_pos.iterrows():
        coord = row['cordinates'].split(',')
        coord = coord[:-1]
        length = row['length'].split(',')
        length = length[:-1]
        start = row['start']
        end = row['end']
        tr = row['transcript_id']
        ch = row['chr']
        st = row['strand']
        x = []
        y = []
        for i in range(len(length)):
            x.append(int(coord[i]) + int(start))
            y.append(int(coord[i]) + int(start) + int(length[i]))
        results_pos.loc[index] = [tr, ch, st, start, end, x, y] 
    test_first = first_DMS.merge(results_pos, on=['transcript_id', 'chr', 'strand'])
    zobry_cds = pd.DataFrame()
    for i in test_first['transcript_id'].unique():
        tmp = test_first[test_first['transcript_id'] == i]
        lol = []
        for m in range(len(tmp.iloc[0]['exons_start'])):
            for j in range(tmp.iloc[0]['exons_start'][m], (tmp.iloc[0]['exons_end'][m])):
                x = pd.DataFrame()
                x = tmp[tmp['start_map'] == j]
                if len(x) > 0:
                    lol.append(x['score'].values[0])
                else:
                    lol.append(0)
        df = pd.DataFrame()
        df[i] = lol
        df1 = df.T
        for n in range(1, 101):
            z = n + 20
            if n == 1:
                df1[str(n) + 'A'] = (len(df)/100) * n
                df1[str(z) + '_final'] = df1.iloc[:, :round(df1[str(n) + 'A'].values[0])].values.mean()
            else:
                df1[str(n) + 'A'] = (len(df)/100) * n
                df1[str(z) + '_final'] = df1.iloc[:, round(df1[str(n-1) + 'A'].values[0]):round(df1[str(n) + 'A'].values[0])].values.mean()
        names1 = list(range(21, 121))
        mystring = '_final'
        names = ["{}{}".format(l,mystring) for l in names1]
        a7a = pd.DataFrame(df1, columns=names)
        zobry_cds = zobry_cds.append(a7a)

    #Negative data CDS
    results_neg = pd.DataFrame(columns=['transcript_id', 'chr', 'strand', 'start', 'end', 'exons_start', 'exons_end'])
    for index, row in CDS_neg.iterrows():
        coord = row['cordinates'].split(',')
        coord = coord[:-1]
        length = row['length'].split(',')
        length = length[:-1]
        start = row['start']
        end = row['end']
        tr = row['transcript_id']
        ch = row['chr']
        st = row['strand']
        x = []
        y = []
        for i in range(len(length)):
            x.append(int(coord[i]) + int(start))
            y.append(int(coord[i]) + int(start) + int(length[i]))
        results_neg.loc[index] = [tr, ch, st, start, end, x, y] 
        
    test_second = first_DMS.merge(results_neg, on=['transcript_id', 'chr', 'strand'])
    zobry_cds_neg = pd.DataFrame()
    for i in test_second['transcript_id'].unique():
        tmp = test_second[test_second['transcript_id'] == i]
        lol = []
        for m in range(len(tmp.iloc[0]['exons_start'])):
            for j in range(tmp.iloc[0]['exons_start'][m], (tmp.iloc[0]['exons_end'][m])):
                x = pd.DataFrame()
                x = tmp[tmp['start_map'] == j]
                if len(x) > 0:
                    lol.append(x['score'].values[0])
                else:
                    lol.append(0)
        df = pd.DataFrame()
        df[i] = lol
        df1 = df.T
        for n in range(1, 101):
            x = 101 - n
            z = x + 20
            if n == 1:
                df1[str(n) + 'A'] = (len(df)/100) * n
                df1[str(z) + '_final'] = df1.iloc[:, :round(df1[str(n) + 'A'].values[0])].values.mean()
            else:
                df1[str(n) + 'A'] = (len(df)/100) * n
                df1[str(z) + '_final'] = df1.iloc[:, round(df1[str(n-1) + 'A'].values[0]):round(df1[str(n) + 'A'].values[0])].values.mean()
        names1 = list(range(21, 121))
        mystring = '_final'
        names = ["{}{}".format(l,mystring) for l in names1]
        a7a = pd.DataFrame(df1, columns=names)
        zobry_cds_neg = zobry_cds_neg.append(a7a)    
    frames = [zobry_cds, zobry_cds_neg]
    cds = pd.concat(frames)    
    
    # 3'UTR preparation
    utr3 = pd.read_csv(utr3, delim_whitespace=True, names=['chr', 'start', 'end', 'transcript_id', 'z', 'strand', 'start1', 'end1', 'n', 'number_of_exons', 'length', 'cordinates'], usecols=['chr', 'start', 'end', 'transcript_id', 'strand', 'number_of_exons', 'length', 'cordinates'])
    utr3['transcript_id'] = utr3['transcript_id'].str.rstrip('_utr3')
    utr3['utr3_length'] = utr3['length'].apply(lambda x: sum(map(int, x.split(',')[:-1])))
    utr3_pos = utr3[utr3['strand'] == '+']
    utr3_neg = utr3[utr3['strand'] == '-']
    #Positive data 3'UTR
    results_pos = pd.DataFrame(columns=['transcript_id', 'chr', 'strand', 'start', 'end', 'exons_start', 'exons_end'])
    for index, row in utr3_pos.iterrows():
        coord = row['cordinates'].split(',')
        coord = coord[:-1]
        length = row['length'].split(',')
        length = length[:-1]
        start = row['start']
        end = row['end']
        tr = row['transcript_id']
        ch = row['chr']
        st = row['strand']
        x = []
        y = []
        for i in range(len(length)):
            x.append(int(coord[i]) + int(start))
            y.append(int(coord[i]) + int(start) + int(length[i]))
        results_pos.loc[index] = [tr, ch, st, start, end, x, y] 
    test_first = first_DMS.merge(results_pos, on=['transcript_id', 'chr', 'strand'])
    zobry_utr3 = pd.DataFrame()
    for i in test_first['transcript_id'].unique():
        tmp = test_first[test_first['transcript_id'] == i]
        lol = []
        for m in range(len(tmp.iloc[0]['exons_start'])):
            for j in range(tmp.iloc[0]['exons_start'][m], (tmp.iloc[0]['exons_end'][m])):
                x = pd.DataFrame()
                x = tmp[tmp['start_map'] == j]
                if len(x) > 0:
                    lol.append(x['score'].values[0])
                else:
                    lol.append(0)
        df = pd.DataFrame()
        df[i] = lol
        df1 = df.T
        for n in range(1, 71):
            z = n + 120
            if n == 1:
                df1[str(n) + 'A'] = (len(df)/70) * n
                df1[str(z) + '_final'] = df1.iloc[:, :round(df1[str(n) + 'A'].values[0])].values.mean()
            else:
                df1[str(n) + 'A'] = (len(df)/70) * n
                df1[str(z) + '_final'] = df1.iloc[:, round(df1[str(n-1) + 'A'].values[0]):round(df1[str(n) + 'A'].values[0])].values.mean()
        names1 = list(range(121, 191))
        mystring = '_final'
        names = ["{}{}".format(l,mystring) for l in names1]
        a7a = pd.DataFrame(df1, columns=names)
        zobry_utr3 = zobry_utr3.append(a7a)
        
    #Negative data utr3
    results_neg = pd.DataFrame(columns=['transcript_id', 'chr', 'strand', 'start', 'end', 'exons_start', 'exons_end'])
    for index, row in utr3_neg.iterrows():
        coord = row['cordinates'].split(',')
        coord = coord[:-1]
        length = row['length'].split(',')
        length = length[:-1]
        start = row['start']
        end = row['end']
        tr = row['transcript_id']
        ch = row['chr']
        st = row['strand']
        x = []
        y = []
        for i in range(len(length)):
            x.append(int(coord[i]) + int(start))
            y.append(int(coord[i]) + int(start) + int(length[i]))
        results_neg.loc[index] = [tr, ch, st, start, end, x, y] 
        
    test_second = first_DMS.merge(results_neg, on=['transcript_id', 'chr', 'strand'])
    zobry_utr3_neg = pd.DataFrame()
    for i in test_second['transcript_id'].unique():
        tmp = test_second[test_second['transcript_id'] == i]
        lol = []
        for m in range(len(tmp.iloc[0]['exons_start'])):
            for j in range(tmp.iloc[0]['exons_start'][m], (tmp.iloc[0]['exons_end'][m])):
                x = pd.DataFrame()
                x = tmp[tmp['start_map'] == j]
                if len(x) > 0:
                    lol.append(x['score'].values[0])
                else:
                    lol.append(0)
        df = pd.DataFrame()
        df[i] = lol
        df1 = df.T
        for n in range(1, 71):
            x = 71 - n
            z = x + 120
            if n == 1:
                df1[str(n) + 'A'] = (len(df)/70) * n
                df1[str(z) + '_final'] = df1.iloc[:, :round(df1[str(n) + 'A'].values[0])].values.mean()
            else:
                df1[str(n) + 'A'] = (len(df)/70) * n
                df1[str(z) + '_final'] = df1.iloc[:, round(df1[str(n-1) + 'A'].values[0]):round(df1[str(n) + 'A'].values[0])].values.mean()
        names1 = list(range(121, 191))
        mystring = '_final'
        names = ["{}{}".format(l,mystring) for l in names1]
        a7a = pd.DataFrame(df1, columns=names)
        zobry_utr3_neg = zobry_utr3_neg.append(a7a)
    
    frames = [zobry_utr3, zobry_utr3_neg]
    utr3 = pd.concat(frames)        
    result = pd.concat([utr5, cds, utr3], axis=1, join="inner")
    result.to_csv(output, index=False, sep='\t')
def getArgs():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-genes', required=True)
    parser.add_argument('-structure', required=True)
    parser.add_argument('-transcript', required=True)
    parser.add_argument('-utr5', required=True)
    parser.add_argument('-cds', required=True)
    parser.add_argument('-utr3', required=True)
    parser.add_argument('-output', required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = getArgs()
    start = time.time()
    structure = structure(args.genes, args.structure, args.transcript, args.utr5, args.cds, args.utr3, args.output)
    end = time.time()
    print ('time elapsed:' + str(end - start))
