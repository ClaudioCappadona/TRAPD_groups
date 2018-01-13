#!/usr/bin/python
import optparse
import sys
import gzip 

#Parse options
parser = optparse.OptionParser()
parser.add_option("-s", "--snpfile", action="store",dest="snpfilename") #File matching SNPs to genes
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename") #Path to vcf file
parser.add_option("-o", "--outfile", action="store",dest="outfilename", default="out.txt") #Output file name 

parser.add_option("--snpformat", action="store",dest="snpformat", default="VCFID") #Field in which to get SNP names. If not VCF ID, then CHR:POS:REF:ALT is used
##parser.add_option("--snpcolname", action="store",dest="snpcolname", default="NA")
parser.add_option("--samplefile", action="store",dest="samplefilename", default="ALL")
##parser.add_option("--recessive", action="store_true",dest="recessive")

#Optional Filters
parser.add_option("--pass", action="store_true", dest="passfilter")
parser.add_option("--maxAF", action="store",dest="maxAF", default=1)
parser.add_option("--maxAC", action="store",dest="maxAC", default=99999)
parser.add_option("--GTfield", action="store",dest="gtfield", default="GT")

options, args = parser.parse_args()
if not options.snpfilename:   # if filename is not given
    parser.error('A file containing a list of SNPs is needed')

if not options.vcffilename:   # if vcf filename is not given
    parser.error('A vcf file is needed')

def findcarriers(vcfline, gtname, snpformat, samplelist):
	#Find the column in the genotype field corresponding to the genotypes
	gtcol=vcfline.split('\t')[8].split(":").index(gtname)

	if snpformat=="VCFID":
		snpid=vcfline.split('\t')[2]
	else:
		snpid=str(vcfline.split('\t')[0]).lstrip("chr")+":"+str(vcfline.split('\t')[1])+":"+str(vcfline.split('\t')[3])+":"+str(vcfline.split('\t')[4])
	
	#Extract genotypes 
	gt=[i.split(':')[gtcol] for i in vcfline.rstrip().split('\t')[9:]]

	#Find carriers
	hetcarriers=[i for i,val in enumerate(gt) if str(val) in ["0/1", "1/0", "0|1", "1|0"]]
	hetcarriers=list(set(hetcarriers) & set(samplelist))
	homcarriers=[i for i,val in enumerate(gt) if str(val) in ["1/1", "1|1"]]
	homcarriers=list(set(homcarriers) & set(samplelist))
	return [hetcarriers, homcarriers]

def findsampleindex(vcfline, samplefilename):
	#This takes the vcf header line and finds the indices corresponding to the individuals present in the sample file
	samplenames=vcfline.rstrip().split('\t')[9:]

	#If User doesn't provide sample list, assume all samples in vcf
	if samplefilename=="ALL":
		sampleindex=range(0, len(samplenames),1)

	#else, find the indices corresponding to the samples in the user-provided list
	else:
		#Generate sample list
		sample_list=[]
		sample_file=open(samplefilename, "r")
		for line_s1 in sample_file:
        			sample_list.append(line_s1.rstrip())
		sample_file.close()
		sampleindex=[i for i,val in enumerate(samplenames) if str(val) in samplelist]
	return sampleindex

def makesnplist(snpfile, snpcolname):
	#Makes a list of SNPs present in the snpfile
	snplist=[]
	#Read in snpfile
	snp_file=open(snpfile, "r")
	
	for line_snp1 in snp_file:
		line_snp=line_snp1.rstrip('\n').split('\t')

		#Find column corresponding to desired snps
		if line_snp[0]!="GENE":
			snplist=snplist+line_snp[snpcol].split(",")
##		if line_snp[0]=="GENE":
##			if snpcolname!="NA":
##				snpcol=line_snp.index(snpcolname)
##			else:
##				snpcol=1
##		elif len(line_snp[snpcol])>1:
##			#Add SNPs to list
##			snplist=snplist+line_snp[snpcol].split(",")

	return set(snplist)
	snp_file.close()


def calculatecount(genesnps, snptable):
	#This will generate an aggregate count for a given gene.
        het_index=[]
	hom_index=[]
        for s in range(0, len(genesnps), 1):
                if genesnps[s] in snptable:
                        tempsnp=genesnps[s]
                        het_index=het_index+snptable[tempsnp][1]
	                hom_index=hom_index+snptable[tempsnp][2]

	
	#Generate number of individuals carrying one variant
        het_ac=len(set([x for x in het_index if het_index.count(x) == 1]))
	ch_ac=len(set([x for x in het_index if het_index.count(x) > 1]))
        hom_ac=len(list(set(hom_index)))
	return [het_ac, ch_ac, hom_ac]

#Make list of all SNPs across all genes present in snpfile
allsnplist=makesnplist(options.snpfilename, options.snpcolname)

#Make a hashtable with keys as each SNP, and stores a list of indices of carriers for that SNP
count_table={} 

#Open vcf file
if str(options.vcffilename)[-3:]==".gz":
	vcffile=gzip.open(options.vcffilename, "rb")
else:
	vcffile=open(options.vcffilename, "r")


for line_vcf1 in vcffile:
	line_vcf=line_vcf1.rstrip().split('\t')
	if line_vcf[0][0]!="#":
		if options.snpformat=="VCF":
			snpid=str(line_vcf[2])
		else: 
			snpid=str(line_vcf[0].lstrip("chr"))+":"+str(line_vcf[1])+":"+str(line_vcf[3])+":"+str(line_vcf[4])
		if snpid in allsnplist:
			count_table[snpid]=[snpid, list(set(sampleindices) & set(templist))]
			counts=findcarriers(line_vcf1, options.gtfield, options.snpcolname, sampleindices)
			count_table[snpid]=[snpid, counts[0], counts[1]]
		
	#Find indices of samples in the sample file
	elif line_vcf[0]=="#CHROM":
		sampleindices=findsampleindex(line_vcf1, options.samplefilename)

vcffile.close()


#Generate output counts
outfile=open(options.outfilename, "w")
snpfile=open(options.snpfilename, "r")
for line_s1 in snpfile:
	line_s=line_s1.rstrip('\n').split('\t')
	if line_s[0]!="GENE":
		genesnplist=list(set(line_s[1].split(',')))
		counts=calculatecount(genesnplist, count_table)
		outfile.write(line_s[0]+"\t"+str(counts[0])+"\t"+str(counts[1])+"\t"+str(counts[2])+'\n')
outfile.close()
snpfile.close()

