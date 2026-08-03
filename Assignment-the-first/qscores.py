#!/usr/bin/env python

import bioinfo
import gzip
import argparse
import matplotlib.pyplot as plt

#Set up arguments for command line input.
def get_args():
    parser = argparse.ArgumentParser(description="A program that takes a FQ file and outputs a histogram of the per-base read quality scores.")
    parser.add_argument("-a", "--file1", help="Name/path of input file")
    parser.add_argument("-l", "--read-length", help="Length of read (int)")
    parser.add_argument("-r", "--read", help="Read number (R1, R2, R3, R4), for graph naming")
    return parser.parse_args()
args = get_args()

# #argparse files
file1 = args.file1
length = int(args.read_length)

#test files
# file1 = "../TEST-input_FASTQ/testin_R1.fq.gz"
# length = 101

qscore_dict: dict = {}
#Keys: Nucleotide position 0-100
#Values: Sum of converted qscore values (ints) (later to be turned into avgs)
for i in range(length):
    qscore_dict[i]=0

#fill dictionary with sums per nt position
with gzip.open(file1, "rt") as fh:
    line_index=0
    record_count=0
    while True:
        line = fh.readline().strip()
        if line_index%4==3: #Qscore
            record_count+=1
            for score_index, score in enumerate(line):
                qscore_dict[score_index]+=bioinfo.convert_phred(score)
        line_index+=1

        if line=="":
            break

#turn sums into averages
for position in qscore_dict:
    qscore_dict[position]= qscore_dict[position]/record_count

#Record table and chart graph
x = []
y = []
for pos in qscore_dict:
    x.append(pos)
    y.append(qscore_dict[pos])

with open(f"tables/{args.read}.txt", "wt") as tablefile:
    for a in range(len(x)):
        tablefile.write(f"{x[a]}\t{y[a]}\n")

plt.bar(x, y, color="cornflowerblue")
plt.title(f"Quality score distribution per nucleotide ({args.read})")
plt.xlabel("Nucleotide Position")
plt.ylabel("Average Quality Score")
plt.savefig(f"Distribution_{args.read}.png")