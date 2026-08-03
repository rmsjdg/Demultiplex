#!/usr/bin/env python

# Author: <Robbie Johnson> <roj@uoregon.edu>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "1.1"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

from typing import TextIO

DNA_bases = set("ATGCN")
RNA_bases = set("AUGCN")
complement_dictionary = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}

def convert_phred(letter: str) -> int:
    '''Converts a single phred value into a quality score.'''
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """Takes a phredscore string as an input parameter and calulates the average quality score of the entire phred string."""
    tot = 0
    for score in phred_score:
        tot += convert_phred(score)
    return tot/len(phred_score)

def validate_base_seq(seq, RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only Ns, As, Gs, Cs, and (Ts if DNA, Us if RNA). False otherwise. Case insensitive.'''
    seq = seq.upper()
    seq = set(seq)
    return seq <= (RNA_bases if RNAflag else DNA_bases)

def gc_content(DNA):
    '''Calculates the GC content of the input.'''
    assert validate_base_seq(DNA)
    DNA = DNA.upper()
    countg = DNA.count("G")
    countc = DNA.count("C")
    return (countg + countc) / len(DNA)

def calc_median(sorted_list):
    '''Given a sorted list, returns the median value of the list'''
    length = len(sorted_list)
    if length %2 == 0: #even length list    
        return (sorted_list[length // 2 - 1] + sorted_list[length // 2]) / 2
    else: #odd length list
        return sorted_list[length // 2]

def one_line_fasta(infile, outfile):
    '''Takes a fasta file and empty file. Populates empty file with contents from fasta file, with only one line of sequence.'''
    with open(infile, "r") as inf:
        with open(outfile, "w") as outf: #Open both files
            seq = "" #Start a string to hold the sequence
            for line in inf: #for each line in our in file 
                if line.startswith(">"): #If it's a new record
                    if seq != "": #and seq is not blank (it holds the entirety of the last record)
                        outf.write(f"{seq}\n") #Write the sequence and a new line so we can move to new record
                        seq = "" #blank sequence bc we finished the record
                    outf.write(line) #
                else: #If the line is not a header
                    seq += (line.strip()) #add it to our string
            outf.write(seq) #write the sequence to the outfile

def reverse_complement(DNA: str) -> str:
    '''Take a DNA sequence and return its reverse compliment. For unknown nucleotides, lists "N" as reverse compliment.'''
    rev=""
    for i in range(len(DNA)):
        rev+=(complement_dictionary[DNA[-(i+1)]])
    return rev

def read_fqrecord(file_handle: TextIO) -> list:
    '''Takes a fastq record and returns a list of the next record (four lines) in a list.'''
    record = []
    for i in range(4):
        record.append(file_handle.readline().strip())
    return record

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    
    #Convert Phred
    assert convert_phred("!") == 0, "wrong phred score for '!'"
    assert convert_phred("+") == 10, "wrong phred score for '+'"
    assert convert_phred("5") == 20, "wrong phred score for '5'"
    assert convert_phred("?") == 30, "wrong phred score for '?'"
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    print("Your convert_phred function is working! Nice job")

    #Qual_score
    assert qual_score("+5") == 15.0, "wrong average phred score for '+5'"
    assert qual_score("HI!!!") == 15.8, "wrong average phred score for 'HI!!!'"
    assert qual_score("(***)") == 8.4, "wrong average phred score for '(***)'"
    assert qual_score("K=-+") == 23, "wrong average phred score for 'K=-+'"
    assert qual_score("H:87") == 27.25, "wrong average phred score for 'H:87'"
    print("Your qual_score function is working! Nice job")

    #Validate_base_seq
    assert validate_base_seq("CGAAGTC"), "Validate base seq does not work on DNA"
    assert validate_base_seq("GUUCAAG", True), "Validate base seq does not work on RNA"
    assert validate_base_seq("Ayo, river!!")==False, "Not a DNA string"
    assert validate_base_seq("taccga"), "Validate base seq does not work on lowercase DNA"
    assert validate_base_seq("uggaau", True), "Validate base seq does not work on lowercase RNA"
    print("Your validate_base_seq function is working! Nice job")

    #Gc_content
    assert gc_content("GGGGGGGGCCCCCCC") == 1
    assert gc_content("TATATATAAAAAATTTTTT") == 0
    assert gc_content("ATCGATCGATCG") == 0.5
    print("Your gc_content function is working! Nice job")

    #calc_median
    assert calc_median([1, 50, 100, 200, 500]) == 100
    assert calc_median([5, 10]) == 7.5
    assert calc_median([2,8,8,8,8,8,100]) == 8
    assert calc_median([1]) == 1
    assert calc_median([500000, 50000000, 50000000000]) == 50000000
    print("Your calc_median function is working! Nice job")