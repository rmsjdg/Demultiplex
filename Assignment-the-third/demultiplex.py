#!/usr/bin/env python

import bioinfo
import argparse
import gzip
import matplotlib.pyplot as plt

def get_args(): #CLI
    parser = argparse.ArgumentParser(description="A program that takes four gzipped Illumina read files (R1-R4), an index file, \
                                     a quailty score cutoff (optional), and an output destination (optional). \
                                     Demultiplexes and organizes outputs data into separate FQ files for each barcode permutation.\
                                     Generates a summary.md file. \
                                     \nNote: Index file structure must include header and barcodes must be in fifth column of each row.")
    parser.add_argument("-a", "--file1", help="Path of R1 input file", required = True, type = str)
    parser.add_argument("-b", "--file2", help="Path of R2 input file", required = True, type = str)
    parser.add_argument("-c", "--file3", help="Path of R3 input file", required = True, type = str)
    parser.add_argument("-d", "--file4", help="Path of R4 input file", required = True, type = str)
    parser.add_argument("-i", "--indexfile", help="Path of index file", required = True, type = str)
    parser.add_argument("-o", "--output", help="Name of pre-existing output folder if desired. \
                        Leave blank to write output files to current dir.", default = "./", required = False, type = str)
    parser.add_argument("-q", "--qcutoff", help="Quality score cutoff level. Leave blank to default to 2.", default = 2, required = False, type = int)
    return parser.parse_args()
args = get_args()

#Reassign argparse to variables
index_file = args.indexfile
file1 = args.file1
file2 = args.file2
file3 = args.file3
file4 = args.file4
qcutoff = args.qcutoff
output = args.output
#ensure output folder includes
if output[-1] != "/":
    output += "/"


# #test files
# index_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"
# file1 = "../TEST-input_FASTQ/testin_R1.fq.gz"
# file2 = "../TEST-input_FASTQ/testin_R2.fq.gz"
# file3 = "../TEST-input_FASTQ/testin_R3.fq.gz"
# file4 = "../TEST-input_FASTQ/testin_R4.fq.gz"
# qcutoff = 20
# output = "output"

# ./demultiplex.py -a ../TEST-input_FASTQ/testin_R1.fq.gz -b ../TEST-input_FASTQ/testin_R2.fq.gz -c ../TEST-input_FASTQ/testin_R3.fq.gz -d ../TEST-input_FASTQ/testin_R4.fq.gz -i /projects/bgmp/shared/2017_sequencing/indexes.txt -o /projects/bgmp/roj/bioinfo/Bi622/Demultiplex/Assignment-the-third/output -q 20

#Initialize barcodes set
barcodes_set: set = set()
with open(index_file) as fh:
    fh.readline()
    for line in fh:
        line = line.strip().split()
        barcode = line[4]
        barcodes_set.add(barcode)

#Sorted list for summary.md
barcodes_sorted_list: list = []
for barcode in barcodes_set:
    barcodes_sorted_list.append(barcode)
barcodes_sorted_list.sort()

#Set up counting dictionary with all combos initalized to 0.
    #Keys - every index combination (1:1, 1:2... 24:23, 24:24) unknown will be counted separaately.
    #Values - Count of each 
counting_dict: dict = {}
for a in barcodes_set:
    for b in barcodes_set:
        counting_dict[a,b] = 0
unknown=0

#Create all output fq files (and leave open), storing file handles in io_dict
io_dict: dict = {}
of_un1 = open(f"{output}unknown_R1.fq", "wt")
of_un2 = open(f"{output}unknown_R2.fq", "wt")
of_hp1 = open(f"{output}hopped_R1.fq", "wt")
of_hp2 = open(f"{output}hopped_R2.fq", "wt")
for barcode in barcodes_set:
    io_dict[barcode]=(open(f"{output}{barcode}_R1.fq", "wt"), open(f"{output}{barcode}_R2.fq", "wt"))

# Main - Open input files and capture a full record for each file. In a try so we close the files regardless.
try:
    with gzip.open(file1, "rt") as r1, \
        gzip.open(file2, "rt") as r2, \
        gzip.open(file3, "rt") as r3, \
        gzip.open(file4, "rt") as r4:
        record_number: int = 0
        print("Reading fq files.", flush=True)
        while True: 
            record_r1 = bioinfo.read_fqrecord(r1)
            record_r2 = bioinfo.read_fqrecord(r2)
            record_r3 = bioinfo.read_fqrecord(r3)
            record_r4 = bioinfo.read_fqrecord(r4) 

            if record_r1[0]=="": #Escape loop if we reach end of file
                break
            record_number += 1
            if (record_number%10000000)==0: #Periodic print statement every 10mil for sanity
                print(f"On record {record_number:,} ({(record_number/363246735)*100:.2f}%).", flush=True)

            record_r3[1] = bioinfo.reverse_complement(record_r3[1]) #turn R3 barcode into revcomp
            bar1=record_r2[1]
            bar2=record_r3[1]

            record_r1[0]+=f" {bar1}-{bar2}" #append barcode combo to headers
            record_r4[0]+=f" {bar1}-{bar2}" #append barcode combo to headers

            #find mis-reads (also identifies N)
            if (bar1 not in barcodes_set) or (bar2 not in barcodes_set):
                unknown += 1
                for i in range(4):
                    of_un1.write(f"{record_r1[i]}\n")
                    of_un2.write(f"{record_r4[i]}\n")

            # #if under cutoff (average for the whole barcode)
            # elif (bioinfo.qual_score(record_r2[3]) <= qcutoff) or (bioinfo.qual_score(record_r3[3]) <= qcutoff):
            #     unknown += 1
            #     for i in range(4):
            #         of_un1.write(f"{record_r1[i]}\n")
            #         of_un2.write(f"{record_r4[i]}\n")

            #alternatate "if under cutoff" but using individual qscores instead of average.
            elif (any(bioinfo.convert_phred(a) <= qcutoff for a in record_r2[3])) or (any(bioinfo.convert_phred(b) <= qcutoff for b in record_r3[3])):
                unknown += 1
                for i in range(4):
                    of_un1.write(f"{record_r1[i]}\n")
                    of_un2.write(f"{record_r4[i]}\n")

            #if matched
            elif bar1 == bar2:
                counting_dict[bar1,bar2] += 1
                for i in range(4):
                    io_dict[bar1][0].write(f"{record_r1[i]}\n")
                    io_dict[bar1][1].write(f"{record_r4[i]}\n")

            #if hopped
            else:
                counting_dict[bar1, bar2] += 1
                for i in range(4):
                    of_hp1.write(f"{record_r1[i]}\n")
                    of_hp2.write(f"{record_r4[i]}\n")

#Systematically close all output fq files
finally:
    for barcode in barcodes_set:
        io_dict[barcode][0].close()
        io_dict[barcode][1].close()
    of_un1.close()
    of_un2.close()
    of_hp1.close()
    of_hp2.close()



print(f"Printing figure and summary file.", flush=True)
#Calculate percent matched and percent hopped
matched = 0
hopped = 0
# go through counting dict an
for a in counting_dict:
    if a[0] == a[1]:
        matched += counting_dict[a]
    elif a[0] != a[1]:
        hopped += counting_dict[a]


y=[]
for a in barcodes_sorted_list:
    y.append((counting_dict[(a,a)]/record_number)*100)
plt.bar(barcodes_sorted_list, y, color="cornflowerblue")
plt.title(f"Percent of matched indexes (Qscore avg cutoff: {qcutoff})")
plt.xlabel("Barcode")
plt.ylabel("Percent (of all records)")
plt.xticks(rotation=60, ha="right")
plt.savefig(f"{output}matched_and_unknown.png", bbox_inches='tight')


#Write out summary file
with open(f"{output}summary.md", "wt") as opf:
    #Basic stats
    ref_num = 0
    opf.write(f"## General Info:\n\
Total number of records: {record_number:,} \n\n\
Total number of matched reads: {matched:,} ({((matched/record_number)*100):.2f}%)\n\n\
Total number of hopped reads: {hopped:,} ({((hopped/record_number)*100):.2f}%)\n\n\
Total number of unknown reads: {unknown:,} ({(unknown/record_number)*100:.2f}%)\n\n")

    #TABLE 1
    opf.write(f"\n\n## Matched Indexes:\n\
![Count of matched/unknown indexes. X is unknown](matched_and_unknown.png)\n\
### Table 1: \n\
| Ref # | Index | Count | Percent |\n\
|----------|----------|----------|----------|\n")
    for a in barcodes_sorted_list:
        ref_num += 1
        opf.write(f"| {ref_num} | {a} | {counting_dict[(a,a)]} | {(counting_dict[(a,a)]/record_number)*100:.4f}% |\n")

    #TABLE 2
    opf.write(f"\n\n## Hopped Indexes:\n\
### Table 2: \n\
| Ref # | Index 1 | Index 2 | Count | Percent |\n\
|----------|----------|----------|----------|----------|\n")
    for a in barcodes_sorted_list:
        for b in barcodes_sorted_list:
            if a != b:
                ref_num +=1
                opf.write(f"| {ref_num} | {a} | {b} | {counting_dict[(a,b)]} | {(counting_dict[(a,b)]/record_number)*100:.4f}% |\n")