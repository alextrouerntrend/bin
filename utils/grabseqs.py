from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
	prog='grabseqs.py',
	usage='''python grabseqs.py -f [Fasta file to subset] -l [list of headers representing desired subset] -o [output subset fasta]''', description='''This script produces a subet of a fasta file based on a text file containing sequence headers of interest.''',
	epilog='''Requirements: Biopython library''')

parser.add_argument("-f", "--fasta", type=str, help="The fasta file to be subset", required=True)
parser.add_argument("-l", "--list", type=str, help="A one-per-line texfile the headers to be taken from the fasta file.", required=True)
parser.add_argument("-o", "--out", type=str, help="The resulting subsetted fasta file.", required=True)

args = parser.parse_args()
fa = args.fasta
li = args.list
ou = args.out

wanted = set()
with open(li) as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)

fasta_sequences = SeqIO.parse(open(fa),'fasta')
with open(ou, "w") as f:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], f, "fasta")
