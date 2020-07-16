from Bio import SeqIO
from statistics import mean 
from statistics import median
import pylab
import argparse

parser = argparse.ArgumentParser(
	prog='lenDist.py',
	usage='''python lenDist.py -i [input multifasta] -o [output directory]''',
	description='''make some histograms that give info about the length distrubutions of a set of genes.''',
	epilog='''Requires Biopython package''')

parser.add_argument("-i", type=str, help="The sequences to be analyzed in fasta format", required=True)
parser.add_argument("-o", type=str, help="The output directory", required=True)

args = parser.parse_args()
inpu = args.i
outpu = args.o

# Read in the sequence length data
sizes = [len(rec) for rec in SeqIO.parse(inpu, "fasta")]


lenmean = mean(sizes)
lenmean = str(round(lenmean, 2))
lenmed = median(sizes)

ogname=inpu.replace('.','/').split("/")[-2]
print(ogname)
savename = outpu+"/lenDist_"+ogname+".png"
print(savename)

# Make a plot

pylab.hist(sizes, bins=20)
pylab.title(ogname+" has %i sequences\nLengths %i to %i\nMean len = %i" \
            % (len(sizes),min(sizes),max(sizes),round(mean(sizes),2)))
pylab.xlabel("Sequence length (bp)")
pylab.ylabel("Count")
pylab.legend("Mean len = {lenmean}\nMedian len = {lenmed}")
pylab.savefig(savename, format='png')

