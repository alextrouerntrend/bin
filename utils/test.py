"""
Created on Wed Jul  3 11:46:29 2019

@author: Alex Trouern-Trend
"""
import argparse
import sys
import numpy as np
import pandas as pd
import os 

parser = argparse.ArgumentParser(
        prog='tandemgenes.py',
        usage='''python tandemgenes.py \
        -g [.gtf containing gene model coordinates. ] \
        -l [list of orthogroups to be assessed] \
        -o [orthofinder outfile (Orthogroups.csv)] \
        -d [output directory]
        -a [type of analysis {1, 2 and 3}, given by their sum]\
        -p [proximity search upstream and downstream of gene model?]''',
        description='''
###############################################################################
# Discovering tandemly arranged paralogs in a genome bears significance towards
#     the genome architecture of the organism and the evolutionary history of
#     that gene family. It might also be significant to the expression
#     regulation of that gene cluster.
###############################################################################
        ''',
        epilog='''
###############################################################################
# FUNCTION #
###############################################################################
# Given a gtf containing the info of each gene model across the species of
#     interest, the results of the orthofinder run, and a list of orthogroups
#     to be investigated, produces the following:
#         4. A file for each orthogroup that gives coordinates of tandemly
#            arranged models ordered into "clusters"
#         2. Another file that considers all of the groups of interest together
#         1. A file that considers other models in proximity to those belonging
#            to the orthogroups of interest
###############################################################################
        ''')

parser.add_argument("-g", "--gtf",
                    type=str,
                    help="The gtf file to be subset",
                    required=True)
parser.add_argument("-l", "--list",
                    type=str,
                    help="One-per-line the orthogroups to-be-assessed.",
                    required=True)
parser.add_argument("-o", "--ortho",
                    type=str,
                    help="The oufile from othofinder analysis.",
                    required=True)
parser.add_argument("-d", "--dir",
                    type=str,
                    help="Output directory to write results of analysis.",
                    required=False)
parser.add_argument("-a", "--anlys",
                    type=int,
                    help="type of analysis to perform.",
                    required=False,
                    default=4)
parser.add_argument("-p",
                    "--prox",
                    type=int,
                    help="Define search distance",
                    required=False,
                    default=10000)

args = parser.parse_args()
gtf = args.gtf  # called
li = args.list  # called
ort = args.ortho  # called
pth = args.dir
ana = args.anlys  # used
proxi = args.prox

proxx = str(proxi)

# Extract genenames of interest from orthofinder results & produce OG-specific\
#     lists of genes.

with open(li) as p:
    orthonames = p.readlines()

orthonames = [x.strip() for x in orthonames]

ogs = pd.read_csv(ort, sep='\t', header=0, index_col=0)

ogpar = ogs.loc[orthonames]

wantgenes = []
oglist = []

for col in ogpar.columns:
    for row in ogpar.index:
        #print(row)
        genelist = ogpar.loc[row, col]
        genelist = str(genelist)
        #print(genelist)
        genelist = genelist.split(', ')
        for i in genelist:
            i = i.strip()
            wantgenes.append(i)
            oglist.append(row)

# Distill gtf file to relevant data, depending on analysis.

orthog = []
scaff = []
genename = []
start = []
stop = []
strand = []

lines = [line.rstrip('\n') for line in open(gtf)]

lgtfp = [1, 3, 5, 7]
sgtfp = [2, 4, 6]
if ana in sgtfp:
    for row in lines:
        row2 = row.split("\t")
        if len(row2) == 9:
            if row2[2] == "gene":
                query = row2[8]
                query = query.strip()
                if query in wantgenes:
                    ogpos = wantgenes.index(query)
                    grabog = oglist[ogpos]
                    scaffold = row2[0]
                    gene = row2[8]
                    stranded = row2[6]
                    if stranded == "+":
                        sta = row2[3]
                        sto = row2[4]
                    elif stranded == "-":
                        sta = row2[4]
                        sto = row2[3]
                    orthog.append(grabog)
                    scaff.append(scaffold)
                    genename.append(gene)
                    start.append(sta)
                    stop.append(sto)
                    strand.append(stranded)
else:
    sys.exit("Analysis argument ('-a', '--anyls') not understood, should be\
             integer 1-7. Type 'python tandemgenes.py -h' for additional\
             details.")

# Subset lists to include only data on scaffolds that contain more
#     than one gene model.

def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

nuorthog = []
nuscaff = []
nugenename = []
nustart = []
nustop = []
nustrand = []

mult = set([x for x in scaff if scaff.count(x) > 1])
for item in mult:
    pos = list_duplicates_of(scaff, item)
    for i in pos:
        nuorthog.append(orthog[i])
        nuscaff.append(scaff[i])
        nugenename.append(genename[i])
        nustart.append(start[i])
        nustop.append(stop[i])
        nustrand.append(strand[i])

# Sanity Check
lenlen = []
lenlen.append(len(orthog))
lenlen.append(len(scaff))
lenlen.append(len(genename))
lenlen.append(len(start))
lenlen.append(len(stop))
lenlen.append(len(strand))

listlen = []
listlen.append(len(nuorthog))
listlen.append(len(nuscaff))
listlen.append(len(nugenename))
listlen.append(len(nustart))
listlen.append(len(nustop))
listlen.append(len(nustrand))

def all_the_same(elements):
   if len(elements) < 1:
       return True
   return len(elements) == elements.count(elements[0])

if all_the_same(lenlen) == False:
    print("Uh, Houston, we've had a problem: lenlen")
    sys.exit()
else:
    print("Sanity check passed!")

if all_the_same(listlen) == False:
    print("Uh, Houston, we've had a problem")
    sys.exit()
else:
    print("Sanity check passed!")

# Small Analysis Dataframe:

df = pd.DataFrame(list(zip(nuscaff, nuorthog, nugenename, nustrand, nustart, nustop)), columns=['Scaffold', 'Orthogroup', 'Gene', 'Strand', 'Start', 'Stop'])

def file_size(fname):
        statinfo = os.stat(fname)
        return statinfo.st_size

# Analyze on per Orthogroup basis.

setogs = df.Orthogroup.unique()

ogoutpath = "orthogroup_tandem"

try:
    os.mkdir(ogoutpath)
except OSError:
    print("Creation of directory %s failed" % ogoutpath)
else:
    print("Creation of directory %s successful" % ogoutpath)

clusthist = []
if ana > 3:    
    for i in setogs:
        clustcount = 0
        places = df.index[df['Orthogroup'] == i]
        subset = df.loc[places]
        setscaff = subset.Scaffold.unique()
        ortoname = str(i)
        outname = ogoutpath+ "/" + ortoname + "_" + proxx + "_tandems.txt"
        f = open(outname, "w")
        f.write("Tandem Arrayed Genes listed by cluster below\n")
        f.close()
        for i in setscaff:
            areas = subset.index[subset['Scaffold'] == i]
            geneset = subset.loc[areas]
            for idx, ro in geneset.iterrows():
                if ro['Strand'] ==  "+":
                    hi = int(ro['Stop'])
                    lo = int(ro['Start'])
                elif ro['Strand'] ==  "-":
                    lo = int(ro['Stop'])
                    hi = int(ro['Start'])
                top = hi + proxi
                bot = lo - proxi
                #print(str(top) + "\t" + str(bot))
                clust = []
                clust.append(ro['Gene'])
                for idx2, ro2 in geneset.iterrows():
                    if ro2['Strand'] ==  "+":
                        if int(ro2['Stop']) <= top and int(ro2['Stop']) > hi:
                            top = int(ro2['Stop']) + proxi
                            clust.append(ro2['Gene'])
                        elif int(ro2['Start']) < lo and  int(ro2['Start']) >= bot:
                            bot = int(ro2['Start']) - proxi
                            clust.append(ro2['Gene'])
                    elif ro2['Strand'] ==  "-":
                        if int(ro2['Stop']) >= bot and int(ro2['Stop']) < lo:
                            bot = int(ro2['Stop']) - proxi
                            clust.append(ro2['Gene'])
                        elif int(ro2['Start']) > hi and  int(ro2['Start']) <= top:
                            top = int(ro2['Start']) + proxi
                            clust.append(ro2['Gene'])
                if len(clust) > 1:
                    if clust[1] not in clusthist:
                        clustcount += 1
                        final = df[df['Gene'].isin(clust)]
                        csv = final.to_csv(index=False, sep = "\t")
                        f = open(ogoutpath+ "/" + ortoname + "_" + proxx + "_tandems.txt", "a")
                        f.write("\n==> " + ortoname + " cluster_" + str(clustcount) + " <==\n")
                        f.write(csv)
                        f.close()
                        clusthist[1:1] = clust
        popul = file_size(outname)
        if int(popul) < 50:
            os.remove(outname)               
# Analyze across orthogroups given by user.

setoutpath = "set_tandem"

try:
    os.mkdir(setoutpath)
except OSError:
    print("Creation of directory %s failed" % setoutpath)
else:
    print("Creation of directory %s successful" % setoutpath)

clustcount = 0
clusthist = []

setscaff = df.Scaffold.unique()

setana = [2, 3, 6, 7]
if ana in setana:
    for i in setscaff:
        clustcount = 0
        scaffname = str(i)
        outname = setoutpath + "/" + scaffname + "_" + proxx + "_tandems.txt"
        f = open(outname, "w")
        f.write("Tandem Arrayed Genes listed by cluster below\n")
        f.close()
        areas = df.index[df['Scaffold'] == i]
        geneset = df.loc[areas]
        for idx, ro in geneset.iterrows():
            if ro['Strand'] ==  "+":
                hi = int(ro['Stop'])
                lo = int(ro['Start'])
            elif ro['Strand'] ==  "-":
                lo = int(ro['Stop'])
                hi = int(ro['Start'])                                                        
            top = hi + proxi
            bot = lo - proxi
            #print(str(top) + "\t" + str(bot))
            clust = []
            clust.append(ro['Gene'])
            for idx2, ro2 in geneset.iterrows():
                if ro2['Strand'] ==  "+":
                    if int(ro2['Stop']) <= top and int(ro2['Stop']) > hi:
                        top = int(ro2['Stop']) + proxi
                        clust.append(ro2['Gene'])
                    elif int(ro2['Start']) < lo and  int(ro2['Start']) >= bot:
                        bot = int(ro2['Start']) - proxi
                        clust.append(ro2['Gene'])
                elif ro2['Strand'] ==  "-":
                    if int(ro2['Stop']) >= bot and int(ro2['Stop']) < lo:
                        bot = int(ro2['Stop']) - proxi
                        clust.append(ro2['Gene'])
                    elif int(ro2['Start']) > hi and  int(ro2['Start']) <= top:
                        top = int(ro2['Start']) + proxi
                        clust.append(ro2['Gene'])
            if len(clust) > 1:
                if clust[1] not in clusthist:
                    clustcount += 1
                    final = df[df['Gene'].isin(clust)]
                    csv = final.to_csv(index=False, sep = "\t")
                    f = open(setoutpath+ "/" + scaffname + "_" + proxx + "_tandems.txt", "a")
                    f.write("\n==>  " + scaffname + " cluster_" + str(clustcount) + "  <==\n")
                    f.write(csv)
                    f.close()
                    clusthist[1:1] = clust
        popul = file_size(outname)
        if int(popul) < 50:
            os.remove(outname)
