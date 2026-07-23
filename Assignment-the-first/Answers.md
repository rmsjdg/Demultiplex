# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz |  |  |  |
| 1294_S1_L008_R2_001.fastq.gz |  |  |  |
| 1294_S1_L008_R3_001.fastq.gz |  |  |  |
| 1294_S1_L008_R4_001.fastq.gz |  |  |  |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**
    
## Part 2
1. Define the problem
    We have four separate files, a read file and index file for read one, and a read file and index file for read two. We need to assess whether the indices match, have hopped, or are not in our list of known indices. We also aim to separate the data into individual files, separated by index. 

2. Describe output
    This code will output 52 files. 48 files containing matched reads (24 from read 1, 24 from read 2), 2 files listing which reads demonstrated index hopping (1 for read 1 one for read 2), and 2 files listing which reads were unable to be matched to an index (read 1 and read 2).

    Additionally, the code will calculate a few statistics on the data, including how many read-pairs had properly matched indices for each potential index pair, how many read-pairs exhibited index hopping, and how many had unknown indices. 

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).

    Using the barcodes CTTA and TGAC, my four test files will result in:
        two records in the 1:1 matching index file (CTTA-TAAG)
        one record in the 2:2 matching index file (TGAC-GTCA)
        one hopped read in the 1:2 index file (CTTA-GTCA)
        and one read in unknown. (NTTA-NAAG)

4. Pseudocode

    [pseudocode.txt](pseudocode.txt)

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
