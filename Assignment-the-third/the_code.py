#!/usr/bin/env python

# Resulting files (x52):
# output/:
#     [barcode]_R1.fq     x24
#     [barcode]_R2.fq     x24
#     hopped_r1.fq
#     hopped_r2.fq
#     unknown_r1.fq
#     unknown_r2.fq

import bioinfo
import argparse
# import itertools

#Set up arguments for command line input.
def get_args():
    parser = argparse.ArgumentParser(description="A program that takes four read files (R1-R4), an index file, and a quailty score ")
    parser.add_argument("-a", "--file1", help="Name/path of R1 input file")
    parser.add_argument("-b", "--file2", help="Name/path of R2 input file")
    parser.add_argument("-c", "--file3", help="Name/path of R3 input file")
    parser.add_argument("-d", "--file4", help="Name/path of R4 input file")
    parser.add_argument("-i", "--indexfile", help="Name/path of index file")
    parser.add_argument("-q", "--qcutoff", help="Quality score cutoff level")
    return parser.parse_args()
args = get_args()

#argparse files
index_file = args.indexfile
file1 = args.file1
file2 = args.file2
file3 = args.file3
file4 = args.file4
qcutoff = args.qcutoff

#test files
index_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"
file1 = "../TEST-input_FASTQ/test_R1.fq"
file2 = "../TEST-input_FASTQ/test_R2.fq"
file3 = "../TEST-input_FASTQ/test_R3.fq"
file4 = "../TEST-input_FASTQ/test_R4.fq"


#Initialize barcodes set
barcodes_set = set()
with open(index_file) as fh:
    fh.readline()
    for line in fh:
        line = line.strip().split()
        barcode = line[4]
        barcodes_set.add(barcode)

#Set up counting dictionary with all combos initalized to 0.
counting_dict = {}
counting_dict["unknown"]=0
    #Keys - "unknown" and every index combination (1:1, 1:2... 24:23, 24:24)
    #Values - Count of each 
for a in barcodes_set:
    for b in barcodes_set:
        counting_dict[f"{a}-{b}"] = 0


#Create all output fq files (and leave open)
io_dict = {}
of_un1 = open("output/unknown_R1.fq", "w")
of_un2 = open("output/unknown_R2.fq", "w")
of_hp1 = open("output/hopped_R1.fq", "w")
of_hp2 = open("output/hopped_R2.fq", "w")
for barcode in barcodes_set:
    io_dict[barcode]=(open(f"output/{barcode}_R1.fq", "w"), open(f"output/{barcode}_R2.fq", "w"))

# Main - Open input files and capture a full record for each file. In a try so we close the files regardless.
try:
    with open(file1) as r1, open(file2) as r2, open(file3) as r3, open(file4) as r4:
        line_number = 0
        while True: 
            record_r1 = bioinfo.read_fqrecord(r1)
            record_r2 = bioinfo.read_fqrecord(r2)
            record_r3 = bioinfo.read_fqrecord(r3)
            record_r4 = bioinfo.read_fqrecord(r4) 

            #We have a whole record
            if line_number % 4 == 3:
                record_r3[1] = bioinfo.reverse_complement(record_r3[1])

                record_r1[0]+=f" {record_r2[1]}-{record_r3[1]}"
                record_r4[0]+=f" {record_r2[1]}-{record_r3[1]}"

                #find mis-reads (also identifies N)
                if record_r2[1] not in barcodes_set or record_r3[1] not in barcodes_set:
                    counting_dict["unknown"] += 1
                    for i in range(4):
                        of_un1.write(f"{record_r1[i]}\n")
                        of_un2.write(f"{record_r4[i]}\n")

                #if under cutoff
                if bioinfo.qual_score(record_r2[1]) <= qcutoff or bioinfo.qual_score(record_r3[1]) <= qcutoff:
                    counting_dict["unknown"] += 1
                    for i in range(4):
                        of_un1.write(f"{record_r1[i]}\n")
                        of_un2.write(f"{record_r4[i]}\n")


                # Write record to Unknown files (fw and rv), appending "[Barcode1]-[Barcode2]"
                # counting_dict[unknown] += 1    

                # print(record_r1)
                # print(record_r2)
                # print(record_r3)
                # print(record_r4)
                record_r1 = []
                record_r2 = []
                record_r3 = []
                record_r4 = []

            # if line_number == 7:
            #     break
            
            
            # if r1 == []:
            #     break

            line_number += 1

            #Cheating way to end - fix later
            if line_number == 20:
                break



#     #find low quality barcodes (average -- may introduce a loop and do convert_phred instead if I want to use individual scores for cutoff)



#     else:
#         if R2 == reverse(R3):
#             Write record to the appropriate Matched file, appending "[Barcode1]-[Barcode2]" in the header
            
#         else: #they do not match 
#             Write record to the Hopped files (fw and rv), appending "[Barcode1]-[Barcode2]" in the header





#Systematically close all output fq files
finally:
    for barcode in barcodes_set:
        io_dict[barcode][0].close()
        io_dict[barcode][1].close()