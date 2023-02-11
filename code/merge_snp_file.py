#!/usr/bin/python
import optparse
import sys 
import pandas as pd
import re
import numpy as np

#Parse options
parser = optparse.OptionParser()
parser.add_option("-s", "--snpfiles", action="store",dest="snpfiles") #Comma separated list of SNP files
parser.add_option("-o", "--outfile", action="store",dest="outfilename") #Output file name 

options, args = parser.parse_args()

#Try to catch potential errors
if not options.snpfiles:   # if filename is not given
    parser.error('A list of SNP files is needed is needed')
    
if options.snpfiles is not None:
    if len(options.snpfiles.split(','))<2:
        parser.error('At least two SNP files are needed')
        
snpfilelist=options.snpfiles.split(',')
snptable={}

for f in range(0,len(snpfilelist),1):
    tempfile=open(snpfilelist[f], "r")
    for line_s1 in tempfile:
        line_s=line_s1.rstrip().split('\t')
        if line_s[0][0]!="#" and len(line_s)>0:
            gene=line_s[0]
            snps=line_s[1:]
            if gene in snptable:
                snptable[gene][1]=snptable[gene][1]+snps
            else:
                snptable[gene]=[gene, snps]
    tempfile.close()
     
#Write Output
outfile=open(options.outfilename, "w")
#outfile.write("#GENE\tSNPS\n")
for x in snptable:
    if len(x)>0:
    #Read through hash table and print out variants
        #snp_out=','.join(snptable[x][1])
        df=pd.DataFrame()
        for ele in snptable[x][1]:
            ser=ele.replace("chr", "")
            ser=ser.replace("X", "23")
            ser=pd.Series(re.split('_|:',ser))
            ser=ser.to_frame().T
            ser[0]=np.array(ser[0].astype(int))
            ser[1]=np.array(ser[1].astype(int))
            df=df.append(ser,ignore_index=True)
        df.columns=["chr", "pos", "vars"]
        df = df.sort_values(['chr','pos'], ascending=True)
        c = [ snptable[x][1][i] for i in df.index]
        #snp_out = '\t'.join(snptable[x][1])
        snp_out = '\t'.join(c)
        outfile.write(str(x)+"\t"+snp_out+"\n")
outfile.close()
                    
                
