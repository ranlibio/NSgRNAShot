About the code
Title: NSgRNAShot
Version: 1
Date: 2020-11-02
Description: This method (NSgRNAShot) will efficiently identify the gRNAs that leads to negative selections from the input read count file. 
Author: Chao Wu <wuchao1984@zju.edu.cn>, Ran Li <ranli1993@zju.edu.cn>
Depends: R (>= 3.6.0).

#####################################################################################################
Running NSgRNAShot is extremely easy and convenient. The demo folder contains one example to go through all steps in NSgRNAShot.

Step 1: prepare the read count file

Remember the following rules of a read count file:
The read count file not include a header line of condition labels;
The first column is the gRNA name;
The second column is pre-selection gRNA count for per gRNA;
The third column is post-selection gRNA count for per gRNA;

As example:

AAAAACACCAGCTGCGAACC	168	219

Step 2: run the NSgRNAShot

NSgRNAShot.R gRNA_count_file prefix_of_ouput

Step 3: output results

If successful, you should see files ".result.txt" and ".log.txt". The top lines of ".result.txt" are as follows:

ACTCAACGTGACCACGCGCC	840	0	1034.92119129577	10	-6.69337710121368	-6.32082705069821	1 demo

Each of the column is:
1. gRNA name
2. P1.count
2. P10.count
3. P1.normalize.count	
4. P10.normalize.count
6. Log2FC
7. Log2FC.centered
8. Label: 1 denote the FDR<0.1 of gRNA, 0 denote the gRNA not cause negative selection.
9. The prefix of output file

For ".log.txt" file, it records the total gRNA number, the threshold value that can be used to distinguish if gRNA causes negative selection, 
the number of gRNA of false positive (FDR<0 and log2FC>0), and the number of gRNA that cause negative selection (FDR<0.1 and log2FC<0).

#######################################################################################################################
How to contact?

If you have any question or suggestion about NSgRNAShot, please contact with Ran Li (ranli1993@zju.edu.cn) and Chao Wu (wuchao1984@zju.edu.cn).
