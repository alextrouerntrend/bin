#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 14:08:13 2020

@author: wegrzynlab
"""
from Bio import SeqIO
from statistics import mean, variance
import pandas as pd
import os

df = pd.DataFrame(columns=['ClustID', 'seqCount', 'avgLen', 'lenVariance', 'species', 'fraser', 'canaan', 'balsam'])

with open("firCentroids90.fnn") as f:
     centroids = [ rec.id for rec in SeqIO.parse(f, "fasta")]
    
for filename in sorted(os.listdir("fir90")):
    print(filename)
    filename = "fir90/" + filename
    sizes = [len(rec) for rec in SeqIO.parse(filename, "fasta")]
    headers = [rec.id for rec in SeqIO.parse(filename, "fasta")]
    total = len(sizes)
    lenmean = mean(sizes)
    if len(sizes) > 1:
        vari = variance(sizes)
    else:
        vari = "NULL"
    ff = 0
    cf = 0
    bf = 0
    spec = ""
    for i in headers:
        if i in centroids:
            centr = i
        if i.startswith("F"):
            ff += 1
            if ff == 1:
                spec += "F"
        elif i.startswith("C"):
            cf += 1
            if cf == 1:
                spec += "C"
        elif i.startswith("B"):
            bf += 1
            if bf == 1:
                spec += "B"
    spec = ''.join(sorted(spec))
    perF = ff/total      
    perC = cf/total
    perB = bf/total
    df = df.append({'ClustID' : filename , 'seqCount' : total , 'avgLen' : lenmean , 'lenVariance' : vari , 'species' : spec , 'fraser' : perF , 'canaan' : perC , 'balsam' : perB}, ignore_index=True)
df.to_csv("vsearchStats.tsv", sep='\t', encoding='utf-8')
