#!/usr/bin/python

import sys, os, gzip, io, datetime, numpy

## honestly, could simplify whole thing
## 1 - take VCF, get header
## 2 - for each line, write as 0, 1, 2, if bi-allelic SNPs

if len(sys.argv) < 3:
    print("plink2treemix.py [gzipped input file] [gzipped output file]")
    print("ERROR: improper command line")
    exit(1)

this_keep_fraction=0.1
outfile = gzip.open(sys.argv[2], "w")

## genotypes are 0/0 -> 2,0
##               0/1 -> 1,1
##               1/1 -> 0,2
def convert_geno(g):
    if (g == "0/0"):
        return "2,0"
    if (g == "0/1"):
        return "1,1"
    if (g == "1/1"):
        return "0,2"
    return "./."


i_line = 0
prev_chr = -1
with io.TextIOWrapper(io.BufferedReader(gzip.open(sys.argv[1], "rb"))) as infile:
    for line in infile:
        i_line = i_line + 1 ## 1-based
        ##line = infile.decode('utf-8')
        if (str(line[0]) == "#") & (str(line[1]) == "C"):
            header_line=line.split("\t")
            header=header_line[9]
            for i_col in range(10, len(header_line)):
                header = header + " " + str(header_line[i_col])
            outfile.write(header.encode())
        line=line.split("\t")            
        if (str(line[0][0]) != "#"):
            chr=line[0]
            if chr != prev_chr:
                print(chr + ", " + str((datetime.datetime.now())))
                prev_chr = chr
            if (len(line[3]) == 1) & (len(line[4]) == 1) & (line[6] == "PASS") & (chr != "X") & (numpy.random.uniform() < this_keep_fraction):
                ## OK - can now strip out!
                has_missing=False
                i_col=9
                g=convert_geno(line[i_col][0:3])                
                if g == "./.":
                    has_missing=True
                out_line=g
                for i_col in range(10, len(line)):
                    g = convert_geno(line[i_col][0:3])
                    out_line = out_line + " " + g
                    if g == "./.":
                        has_missing=True
                out_line = out_line + "\n"
                if has_missing == False:
                    outfile.write(out_line.encode())
        if i_line > 1000:
            done = 1


infile.close()
outfile.close()
exit()

## testing
## python ~/personal/proj/primates/vcf2treemix.py temp.input.vcf.gz temp.vcf.gz && zcat temp.vcf.gz | head
