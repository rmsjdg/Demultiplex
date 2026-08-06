#!/usr/bin/env python

import bioinfo
import argparse
import gzip

#Set up arguments for command line input
def get_args():
    parser = argparse.ArgumentParser(description="A program that takes four read files (R1-R4), an index file, \
                                     and a quailty score cutoff. Creates 54 fq files with barcode-storted data (fw and rv)\
                                     and a summary .md file summarizing hopping data.")
    parser.add_argument("-a", "--file1", help="Path of R1 input file", required = False, type = str)
    parser.add_argument("-b", "--file2", help="Path of R2 input file", required = False, type = str)
    parser.add_argument("-c", "--file3", help="Path of R3 input file", required = False, type = str)
    parser.add_argument("-d", "--file4", help="Path of R4 input file", required = False, type = str)
    parser.add_argument("-i", "--indexfile", help="Path of index file", required = False, type = str)
    parser.add_argument("-o", "--output", help="Name of output folder. Do not include final slash. \
                        Leave blank to write output files to current dir.", \
                        default = ".", required = False, type = str)
    parser.add_argument("-q", "--qcutoff", help="Quality score cutoff level", required = False, type = int)
    return parser.parse_args()
args = get_args()

# argparse files
index_file = args.indexfile
file1 = args.file1
file2 = args.file2
file3 = args.file3
file4 = args.file4
qcutoff = args.qcutoff
output = args.output

#test files
index_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"
file1 = "../TEST-input_FASTQ/testin_R1.fq.gz"
file2 = "../TEST-input_FASTQ/testin_R2.fq.gz"
file3 = "../TEST-input_FASTQ/testin_R3.fq.gz"
file4 = "../TEST-input_FASTQ/testin_R4.fq.gz"
qcutoff = 20
output = "output"


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
of_un1 = open(f"{output}/unknown_R1.fq", "w")
of_un2 = open(f"{output}/unknown_R2.fq", "w")
of_hp1 = open(f"{output}/hopped_R1.fq", "w")
of_hp2 = open(f"{output}/hopped_R2.fq", "w")
for barcode in barcodes_set:
    io_dict[barcode]=(open(f"{output}/{barcode}_R1.fq", "w"), open(f"{output}/{barcode}_R2.fq", "w"))

# Main - Open input files and capture a full record for each file. In a try so we close the files regardless.
try:
    with gzip.open(file1, "rt") as r1, \
        gzip.open(file2, "rt") as r2, \
        gzip.open(file3, "rt") as r3, \
        gzip.open(file4, "rt") as r4:
        record_number: int = 0
        while True: 
            record_r1 = bioinfo.read_fqrecord(r1)
            record_r2 = bioinfo.read_fqrecord(r2)
            record_r3 = bioinfo.read_fqrecord(r3)
            record_r4 = bioinfo.read_fqrecord(r4) 

            if record_r1[0]=="":
                break
            record_number += 1

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

            #if under cutoff
            elif (bioinfo.qual_score(record_r2[3]) <= qcutoff) or (bioinfo.qual_score(record_r3[3]) <= qcutoff):
                unknown += 1
                for i in range(4):
                    of_un1.write(f"{record_r1[i]}\n")
                    of_un2.write(f"{record_r4[i]}\n")

            # #alternatate "if under cutoff" but using individual qscores instead of average.
            # if (any(bioinfo.convert_phred(a) for a in record_r2[3]) <= qcutoff) or (any(bioinfo.convert_phred(b) for b in record_r3[3]) <= qcutoff):
            #     pass


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



# Probably need to add some extra steps in here to manipulate the data so I can appropriately report all required info.


#Calculate percent matched and percent hopped
matched = 0
hopped = 0
# go through counting dict an
for a in counting_dict:
    if a[0] == a[1]:
        matched += counting_dict[a]
    elif a[0] != a[1]:
        hopped += counting_dict[a]

#Write out summary file
with open(f"{output}/summary.md", "wt") as opf:
    #Basic stats
    opf.write(f"## Basic Stats:\nTotal number of records: {record_number}\n\nTotal number of matched reads: {matched} ({((matched/record_number)*100):.2f}%)\n\nTotal number of hopped reads: {hopped} ({((hopped/record_number)*100):.2f}%)\n\nTotal number of unknown reads: {unknown} ({(unknown/record_number)*100:.2f}%)\n\n")


    #TABLE 1
    opf.write(f"\n\n## Table 1: Matched Indexes:\n\
| Index | Count | Percent |\n\
|----------|----------|----------|\n")

    for a in barcodes_sorted_list:
        opf.write(f"| {a} | {counting_dict[(a,a)]} | {(counting_dict[(a,a)]/record_number)*100:.2f}% |\n")

    #TABLE 2
    opf.write(f"\n\n## Table 2: Hopped Indexes:\n\
| Index | Count | Percent |\n\
|----------|----------|----------|\n")
    for a in barcodes_sorted_list:
        for b in barcodes_sorted_list:
            if a != b:
                opf.write(f"| {a} | {b} | {counting_dict[(a,b)]} | {(counting_dict[(a,b)]/record_number)*100:.2f}% |\n")