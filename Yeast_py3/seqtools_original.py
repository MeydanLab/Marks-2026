# This is a collection of programs for processing fasta, fastq, SAM, and GFF files.
# The more critical functions here are writegene2 and givegene, and some GFF processing tools.

from Bio import SeqIO
from Bio import Seq
import builddense
from Bio import SeqRecord
from BCBio import GFF
import csv
import gentools


# Simple program that reports read sizes.
#recordsin is a generator.
def countreadsizes(recordsin,lowest,highest):
	sizes=[0 for x in range(highest-lowest+1)]
	for record in recordsin:
		reclen=len(record)
		if reclen>highest or reclen<lowest:
			print("ERROR - wrong lengths")
			exit()
		else:
			sizes[reclen-lowest]+=1
	return sizes
	
# Counts up reads as one whole for every suffix of a particular file.
def countreadsizes_wf(filesin,suffixes,fileout,lowest,highest):
	inlists=[]
	if suffixes==[]:		# Default option.
		for filecount in filesin:
			suffixes.append(".fastq")
	if len(suffixes)!=len(filesin):
		print("ERROR IN SUFFIX LIST LENGTH.")
		exit()
		
	for i in range(len(filesin)):
		filein=filesin[i]
		suffix=suffixes[i]
		inlist=[0 for x in range(highest-lowest+1)]
	# Go through each group of files and get sizes, then combine together.
		
		f_in=open(filein+suffix)
		recin=SeqIO.parse(f_in,"fastq")
		# Get counts for each and combine to master
		
		retlist=countreadsizes(recin,lowest,highest)
		for i in range(len(inlist)):
			inlist[i]+=retlist[i]
		
		inlist.reverse()
		fileinname=filein[filein.rfind("/")+1:]
		inlist.append(fileinname+"_"+str(lowest)+"-"+str(highest))
		inlist.reverse()
		inlists.append(inlist)
		f_in.close()
	gentools.writelisttoexcel(inlists,fileout)


# Simple program that writes out reads with all A's on the 3' end past a particular point.
# input and output paths are fastq files.
# bases is number of bases on 5' end that match to the genome.
# A _run script calls this for several fastq files and combines them to a unified file for bowtie.
# Returns only the nonA reads.
def outAreads_wf(fastapathin,fastapathout,bases):
	
	f_in=open(fastapathin)
	recin=SeqIO.parse(f_in,"fastq")
	f_out=open(fastapathout,"w")
	
	reads=outAreads(recin,bases)
	
	
	SeqIO.write(reads, f_out, "fastq")
	f_in.close()
	f_out.close()

def outAreads(recordsin,bases):
	for record in recordsin:
		reclen=len(record)
		if reclen>bases:	#Make sure read has at least on extra base on end.
			trimmedread=record[bases:]
			reclen-=bases
			Acount=0
			while(trimmedread[Acount]=='A'):
				Acount+=1
				if(Acount==reclen):
					break
			if(Acount==reclen):
				yield record[0:bases]		# This can be used to output trimmed or untrimmed reads.

# Just like outAreads, but now it outputs reads with consecutive As starting at least so many bases from 3' end.
def outAreads3(recordsin,bases):
	for record in recordsin:
		reclen=len(record)
		Acount=1
		while(record[-Acount]=='A'):
			Acount+=1
			if(Acount>reclen):
				break
		Acount-=1	# Because we started with 1.
		if(Acount>=bases):
			yield record[0:-Acount]		# Can have -Acount as 3' end or can have just full record yielded.

def outAreads3_wf(fastapathin,fastapathout,bases):
	
	f_in=open(fastapathin)
	recin=SeqIO.parse(f_in,"fastq")
	f_out=open(fastapathout,"w")
	
	reads=outAreads3(recin,bases)
	
	SeqIO.write(reads, f_out, "fastq")
	f_in.close()
	f_out.close()


# This function counts the number of reads with > than some number of contiguous A's (from the right) in a given seed region.
# f_out is a fastq file of reads selected for some number of polyA.
# Put in 0 for seelen if you want the whole read.
# Set fout to 0 to not write out anything.
# Set numA to -1 to count reads that are all A (as well as seedlen=0).
def countA(fastafile,numA,seedlen,f_out):
	count = 0
	totcount=0
	totAcount = 0
	reads = SeqIO.parse(fastafile,"fastq")
	
	for read in reads:
		if seedlen == 0:
			tr=read
		else:	
			tr = read[0:seedlen]
		trimmedread=tr[::-1]
		seedlen=len(tr)		
		totAcount += trimmedread.seq.count('A')
		
		Acount=0
		while(trimmedread[Acount]=='A'):
			Acount+=1
			if(Acount==seedlen):
				break
		
		if numA==-1:
			numA=seedlen-1
		
		if Acount > numA:
			count += 1
			if f_out!=0:
				SeqIO.write(read, f_out, "fastq")
				# write out read	
		totcount+=1
	
	percentA = float(totAcount)/(totcount*seedlen)
	
	print("Number of reads with enough consecutive A's at 3' end = "+str(count))
	print("  in percent terms = "+str(float(count)/totcount))
	print("Percent A in full seed = "+str(percentA))
	print("Total readcount = "+str(totcount))


# This function will go through a fastq and determine the number of reads with exclusively A's 3' of each position in the read.
# bases is the maximum number of bases to count in from the 3' end.
# This is newer than countA.
# New feature. If bases is negative, works just like it does when positive, except that bases is turned positive but only reads of that length and shorter are considered.
#New feature. Now we count from 3 or 5 depending on param whichend being 1 or -1. -1 is for 3' end.
def countallA(fastqfilename,bases,whichend):
	f=open(fastqfilename)
	reads=SeqIO.parse(f,"fastq")
	lencheck=0
	if bases<0:
		bases=-bases
		lencheck=1
	
	
	Aarray=[0 for x in range(bases+1)]

	for read in reads:
		readlen=len(read)
		if lencheck==1 and readlen>bases:
			continue
		
		Acount=1
		
		while(read[whichend*Acount]=='A'):
				Acount+=1
				if(Acount==bases+1):
					break
				if(Acount>readlen):
					break
		Aarray[Acount-1]+=1
	
	f.close()
	return Aarray # Aarray 1st position is 0 As on end, 2nd is 1, etc.


# Workflow for CountA
def countallA_wf(filelist,bases,outfile,whichend):
	writer = csv.writer(open(outfile+".csv", "wb"))
	for fastqfile in filelist:
		print(fastqfile)
		arr=countallA(fastqfile,bases,whichend)
		#writer.writerow(fastqfile)
		arr.insert(0,fastqfile)
		writer.writerow(arr)

			
				
# This function outputs a fasta file that can be used by bowtie to map reads against. 
#seqid is the short sequence id.
# sequence is the sequence, descrip is the description, and f_fsafilepath is an open text file where it will be written.
def createfasta_adhocmsg(sequence,descrip,seqid,f_fsafilepath):

    seqrec=SeqRecord.SeqRecord(Seq(sequence))    
    seqrec.id=seqid
    seqrec.description=descrip
    
    f_fsafilepath.write(seqrec.format("fasta"))
                        
    f_fsafilepath.close()



# Tool for loading the entire yeast genome into memory.
def makeGFFlist(GFFgen):
	GFFlist={}
	for chr in GFFgen:
		GFFlist[chr.id]=chr
	return GFFlist

# Table generator for UTR start and end points (like a GFFlist for UTRs).
# utrgffgen is the gffgen for 3'UTR from Nagalkshmi for R64.
# For pombe, it takes the normal GFF file that has been modified to only have UTR3 or UTR5 entries (since this script looks for either).
def makeutrtable(utrgffgen):
	table={}
	for chr in utrgffgen:
		if chr.id!="I" and chr.id!="II" and chr.id!="III" and chr.id!="AB325691" and chr.id!="MTR" and chr.id!="MT":
			for feature in chr.features:
				stoppoint=(feature.id).find("_")
				if stoppoint==-1:
					print("ERROR in reading UTR GFF file - no _ character")
				else:
					table[feature.id[0:stoppoint]]=[feature.location.start.position,feature.location.end.position]	# Not sure what this was about in Scer.
		else:
			for feature in chr.features:
				table[feature.id]=[feature.location.start.position,feature.location.end.position]
	return table

# This function creates a dictionary for looking up the GFF feature and chrom for gene names. Also includes other gene info about each gene.
def makeidtable(GFFgen):
	table={}
	for chrom in GFFgen:
		num_feat=0
		for feature in chrom.features:
			name_feat=feature.id
			table[name_feat]=[feature,num_feat,chrom.id]	
			num_feat+=1
	return table

#same thing for aliases
def makeidtable2(GFFgen):
	table={}
	for chrom in GFFgen:
		if chrom.id=="gi|49175990|ref|NC_000913.2|":
			keyname="Name"
		# for Scer
		elif chrom.id!="I" and chrom.id!="II" and chrom.id!="III" and chrom.id!="AB325691" and chrom.id!="MTR" and chrom.id!="MT":
			keyname="Alias"
		# For pombe
		else:
			keyname="external_name"
		num_feat=0
		for feature in chrom.features:
			name_feat=feature.id
			for item in feature.qualifiers:
				if keyname in feature.qualifiers:
					alias_feat = feature.qualifiers[keyname][0]
					table[alias_feat]=[feature,num_feat,chrom.id]
			num_feat+=1
	return table
	

# wrapper func
def convertmrnatogenomic(position,chrom,feature,GFFlist):
	if chrom!="I" and chrom!="II" and chrom!="III" and chrom!="AB325691" and chrom!="MTR" and chrom!="MT":
		return convertmrnatogenomicscec(position,chrom,feature,GFFlist)
	else:
		import seqtools_extra
		return seqtools_extra.convertmrnatogenomicpombe(position,chrom,feature,GFFlist)

# This function has been tested and works well. It works only on CDS positions or 3'UTR if position beyond the end. 5'UTR will be correct if value is negative since all position values are with respect to START codon.
# Findseqinstance calls it. I think that's the only thing.
def convertmrnatogenomicscec(position,chrom,feature,GFFlist):
	strand=GFFlist[chrom].features[feature].strand
	# First create 2 lists
	exonlengths=[]    # lengths of exons in gene.
	intronlengths=[]	# lengths of all introns prior to exon. First in list is by definition therefore 0.
	genestart=int(GFFlist[chrom].features[feature].location.start.position)
	geneend=int(GFFlist[chrom].features[feature].location.end.position)
	intronlengths.append(0)
	
	# Check for coli - assume no introns - just add 1st position.. Legacy feature.
	if "gi|49175990|ref|NC_000913.2|" in GFFlist:
		if strand==-1:
			return geneend-position-1		# -1 needed because end position is 1 past actual end of gene as always with GFF parse info.
		else:
			return position+genestart
	
	
	fixswitch=0
	for item in GFFlist[chrom].features[feature].sub_features:
		if item.type == 'intron':		# This 5' UTR intron added for genes starting with intron, previously these were refjected, 9/22/12.
			start = int(item.location.start.position)
			end = int(item.location.end.position)
			intronlengths.append(end-start)
		if item.type == 'CDS':
			start = int(item.location.start.position)
			end = int(item.location.end.position)
			exonlengths.append(end-start)
		if item.type == "five_prime_UTR_intron":		# genes that start with intron. This was fixed on 9/22/12. This is not in pombe.
			fixswitch=1
			startutr5intron = int(item.location.start.position)
			endutr5intron = int(item.location.end.position)
			UTR5intronlen=endutr5intron-startutr5intron
			
	if fixswitch==1:
		for item in GFFlist[chrom].features[feature].sub_features:
			if item.type == 'CDS':
				if strand==1:
					genestart=int(item.location.start.position)
					exon2len=genestart-endutr5intron
				else:
					geneend=int(item.location.end.position)
					exon2len=startutr5intron-geneend

	# Deal with genes that have a break between exons, such as a "region."
	
	weirdexonlength=0
	if len(intronlengths)!=len(exonlengths):
		for item in GFFlist[chrom].features[feature].sub_features:
			start = int(item.location.start.position)
			end = int(item.location.end.position)
			if item.type == 'intron':
				print("Got an intron when looking for a region. Didn't expect that.")
			if item.type == 'CDS' or item.type == 'region':
				weirdexonlength+=(end-start)
		exonlengths=[weirdexonlength]
		intronlengths=[0]	
		
	# This weird gene has two exons that adjoin and no intron. Just a "region" in between. Strange.
	#if GFFlist[chrom].features[feature].id=="YIL009C-A":
	#	intronlengths=[0]
	#	exonlengths=[547]
		
	# This too.
	#if GFFlist[chrom].features[feature].id=="YOR239W":
	#	intronlengths=[0]
	#	exonlengths=[1888]

	if position<0:		# Case of 5'UTR
		if fixswitch==1:
			if abs(position)>exon2len:
				position-=UTR5intronlen		
	
	if len(intronlengths)!=len(exonlengths):
		print("error in intron/exon count.")
		print(GFFlist[chrom].features[feature].id)
	
	if strand==-1:
		del intronlengths[0]
		intronlengths.append(0)
		intronlengths.reverse()
		exonlengths.reverse()
	
	# Now convert position.	
	i=0
	runninglength=0
	while i<(len(intronlengths)):
		runninglength+=exonlengths[i]
		if position<runninglength:	
			if strand==1:
				return (genestart+position+sum(intronlengths[0:i+1]))
			if strand==-1:
				return (geneend-position-sum(intronlengths[0:i+1])-1)
		i+=1
	# If we get here, it's because we're asking for a position beyond the end of the gene.
	i-=1
	if strand==1:
		return (genestart+position+sum(intronlengths[0:i+1]))
	if strand==-1:
		return (geneend-position-sum(intronlengths[0:i+1])-1)
	
	
	#print "loop failure"	# Should never get here because we should return something.
	#print GFFlist[chrom].features[feature].id	
	#print position
		
		
		
		
# Wrapper function to call givegene, outputs a csv file of a given gene's density file.
# See note below on using this with formal names as opposed to aliases.
# Options are:
# Shift values tell you how much to put on the 5' and 3' end of the gene extra.
# densityfiles are the input density file for + and - strand.
# features are the names of the genes (can be alias or formal).
# 3 gff files are inputs also.
# outfile is where to write output as a csv.
# intron=0; standard way with spliced transcript returned. (Sense direction, 5' to 3')
# intron=1; Unspliced transcript is returned. (Sense direction, 5' to 3')
# intron=2; Unspliced transcript returned on sense strand, 5' to 3' direction. And, also returned, unspliced transcript returned on antisense strand, 5' to 3' direction.
# intron=3; same as intron = 2 except the orientation for both sense and antisense will be as they appear on the chromosome. Nothing is flipped.

def writegene2_wf(shift5,shift3,intron,densityfiles,features,gfffile,utrgfffilename5,utrgfffilename3,outfile):	
	writerfile=open(outfile+".csv", "wb")
	writer = csv.writer(writerfile)
	GFFgen=GFF.parse(gfffile)
	tablealias=makeidtable2(GFFgen)
	GFFgen=GFF.parse(gfffile)
	tableformal=makeidtable(GFFgen)
	GFFgen=GFF.parse(gfffile)
	GFFlist=makeGFFlist(GFFgen)
	utrgffgen5=GFF.parse(utrgfffilename5)
	utrtable5=makeutrtable(utrgffgen5)
	utrgffgen3=GFF.parse(utrgfffilename3)
	utrtable3=makeutrtable(utrgffgen3)
	GFFs=[GFFlist,utrtable5,utrtable3]
	
	for densityfile in densityfiles:
		counts1=builddense.readcountsf(densityfile+"_plus_")
		counts2=builddense.readcountsf(densityfile+"_minus_")
		counts=[counts1,counts2]
		for feature in features:
			if len(feature)>=7:
				table=tableformal
			else:
				table=tablealias
			chromosome=table[feature][2]
			featurenum=table[feature][1]
			longfeature=table[feature][0].id
			
			if intron==0:
				retval=givegene(chromosome,featurenum,GFFs,counts,[shift5,shift3,0],0)
			else:
				if(GFFlist[chromosome].features[featurenum].strand==1):
					strand=0
				else:
					strand=1			 
			
				if shift5=="none":
					shift5=0
				if shift3=="none":
					shift3=0
			
				# Start and stop come from UTR annotation if it's there. If not, get start/end of "CDS" (exon) region using a loop.
				if longfeature in utrtable5:
					if strand==0:
						start=int(utrtable5[longfeature][0])
					else:
						end=int(utrtable5[longfeature][0])
				else:	# Get ends of the exon since we don't want the annotated UTR5 intron, if present.
					start=-1
					for item in GFFlist[chromosome].features[featurenum].sub_features:
						if item.type == 'CDS':
							start_feat = int(item.location.start.position)
							if start==-1:
								start=start_feat
												
				if longfeature in utrtable3:
					if strand==0:
						end=int(utrtable3[longfeature][1])
					else:
						start=int(utrtable3[longfeature][1])
				else:	# Get ends of the exon since we don't want the annotated UTR5 intron, if present.
					for item in GFFlist[chromosome].features[featurenum].sub_features:
						if item.type == 'CDS':
							end_feat = int(item.location.end.position)
					end=end_feat
			
				# Flip 3' and 5' extensions if opposite strand. Absolute value since manual UTRs come in negative.
				if strand==1:
					localshift3=abs(shift5)
					localshift5=abs(shift3)
				else:
					localshift3=abs(shift3)
					localshift5=abs(shift5)
					
				unsplicedcounts=counts[strand][chromosome][start-localshift5:end+localshift3]
				antisense_unsplicedcounts=counts[abs(strand-1)][chromosome][start-localshift3:end+localshift5] # Line for getting out antisense reads if requested.
				
				if intron == 1 or intron == 2:
					antisense_unsplicedcounts.reverse()
				
				if strand==1 and (intron == 1 or intron == 2):
					unsplicedcounts.reverse()
					antisense_unsplicedcounts.reverse()
				retval2=[antisense_unsplicedcounts]
				retval=[unsplicedcounts]
			
			writer.writerow([feature+"_"+densityfile.split("/")[-1]]+retval[0])
			if intron==2 or intron==3:
				writer.writerow([feature+"_"+densityfile.split("/")[-1]]+retval2[0])
	writerfile.close()	
	gentools.transposecsv(outfile)	


# This is a wrapper for givegene for distinguishing between pombe and cerevisiae. Coli was supported at one time, but not at this time.
# pombe genome signalled by name of chromosome: for pombe, these are either I, II, or III.

def givegene(chromosome,feature,GFFs,counts,bp,goodgenes):

	if chromosome!="I" and chromosome!="II" and chromosome!="III" and chromosome!="AB325691" and chromosome!="MTR" and chromosome!="MT":
		return(givegene_ScEc(chromosome,feature,GFFs,counts,bp,goodgenes))
	else:
		import seqtools_extra
		return(seqtools_extra.givegene_Sp(chromosome,feature,GFFs,counts,bp,goodgenes))



# This is a function that gives spliced sequences and read counts for full genes. 
# This function takes as input a GFFlist and a dictionary of counts (list with plus and minus tracks) for a genome and a chromosome and feature id.
# The switch goodgenes when set to 1 will output [-1,-1] when a gene is requested that is not: 
# on one of the main chromosomes, if extra region causes it to overlap with another gene, non genes, non dubious genes.
# Output is a list where first element is the counts and 2nd element is the sequence. 
# Third and fourth element are counts/sequence of first intron. Fifth is frame of the intron wrt to main ORF. Last 2 elements are the UTR lengths.
# If goodgenes=2, then overlapping genes are allowed.
# This program will also not return sequences for genes that lack the CDS field. 
# The input list called GFFs has element 0 the GFFlist and 1 and 2 are the utr tables for 5' and 3' utrs, respectively. 
# givegene will compute UTR lengths based on annotation (bp is 0) and extend to value if >0.
# UTR lengths in the output. If no annotation, then length of UTR will be reported as 0 (or the extended value). 
# bp values are the desired extensions to annotated UTRs. Put negative values if you want fully manual UTRs. Put "none" for no UTR at all.
# Note that genes with annotated frameshifts are spliced as well. I have counted 3.
def givegene_ScEc(chromosome,feature,GFFs,counts,bp,goodgenes):
	
	# For no UTR option.
	if bp[0]=="none":
		bp[0]=0
		manualutr5=1
	elif bp[0]<0:
		bp[0]*=-1
		manualutr5=1
	else:
		manualutr5=0
		
	if bp[1]=="none":
		bp[1]=0
		manualutr3=1
	elif bp[1]<0:
		bp[1]*=-1
		manualutr3=1
	else:
		manualutr3=0	
		
	# Option here to handle cases from old functions where shift not given.
	if type(bp)!=int:
		if len(bp)==2:
			bp.append(0)	# Shift not given. Set it to 0.
	else:
		print("Old function needs upgrading.")
		exit()
	
	if type(GFFs)==list:		
		GFFlist=GFFs[0]
		utrtable5=GFFs[1]
		utrtable3=GFFs[2]
	else:
		print("No longer allowed to have no UTR information. Required for overlap checks. July 2020.")
		exit()
				
	# Put in dummy counts if none given - function is then just being used for sequence information.
	if counts==-1:
		outputdata1={}
		outputdata2={}
		for key in list(GFFlist.keys()):
			outputdata1[key]=[0 for x in range(len(GFFlist[key]))]
			outputdata2[key]=[0 for x in range(len(GFFlist[key]))]
		counts=[outputdata1,outputdata2]
	
	if feature>=len(GFFlist[chromosome].features):
		print("not a feature")
		return [-1,-1,-1,-1,-1,-1,-1]
	
	if(GFFlist[chromosome].features[feature].strand==1):
		strand=0
	else:
		strand=1
	
	# July 2020, givegene will not work for mito, 2um, or non-gene entities.
	if (GFFlist[chromosome].id == 'chrMito' or GFFlist[chromosome].id == '2-micron'):
		return [-1,-1,-1,-1,-1,-1,-1]
	# Note that the GFF annotation file as of 2/14/13 has a number (a dozen or more) genes that are flagged here because they are called "region" rather than "gene".
	if (GFFlist[chromosome].features[feature].type!='gene'):
		return [-1,-1,-1,-1,-1,-1,-1]
		
	if goodgenes > 0:		
		if ("orf_classification" in GFFlist[chromosome].features[feature].qualifiers):
			if (GFFlist[chromosome].features[feature].qualifiers["orf_classification"][0]=="Dubious"):
				return [-1,-1,-1,-1,-1,-1,-1]
				
		### ADDED 12/15/16:
		### THESE OBVIOUSLY BAD ANNOTATIONS WILL BE REJECTED NOW 
		#badgenelist=['YOR186W','YGL033W','YJR120W','YHR215W']
		badgenelist=['YGL033W','YJR120W']
		for badgene in badgenelist:
			if GFFlist[chromosome].features[feature].id==badgene:
				return [-1,-1,-1,-1,-1,-1,-1]

	# Do the splicing.
	splicedseq=Seq.Seq('')
	splicedcounts=[]
	start=-1
	#Assemble gene by splicing out introns here. And get start and end of gene (start/stop codon). Must do it here since main gene start/stop include 5'UTR intron.
	for item in GFFlist[chromosome].features[feature].sub_features:
		if item.type == 'CDS':
			start_feat = int(item.location.start.position)
			end_feat = int(item.location.end.position)
			splicedseq+=(GFFlist[chromosome][start_feat:end_feat]).seq
			splicedcounts+=counts[strand][chromosome][start_feat:end_feat]	
			if start==-1:
				start=start_feat		# Set the start of CDS
	end=end_feat	# Set the end of the CDS
			
	# Check to see if there were no CDS entries:
	if splicedcounts==[] or str(splicedseq)=='':
		return [-1,-1,-1,-1,-1,-1,-1]
			
	CDSlen=len(splicedseq)		
	
	# Put on annotated UTRs. Record UTR lengths
	featid=GFFlist[chromosome].features[feature].id
	
	# Check if 5'UTR intron. If yes, handle it by adding in 2nd exon. If 2nd exon is 0 length, nothing changes.
	exon2len=-1	# Set to -1 as indicator of no 5'UTR intron. OTher intron values should never be accessed if this condition true.

	# We are adding on exon2.
	for item in GFFlist[chromosome].features[feature].sub_features:
		if item.type == "five_prime_UTR_intron":
			intronstart = int(item.location.start.position)
			intronend = int(item.location.end.position)
			intronlen=intronend-intronstart
			if strand==0:
				exon2len=start-intronend
			else:
				exon2len=intronstart-end

			if strand==0:
				splicedcounts=counts[strand][chromosome][intronend:start]+splicedcounts
				splicedseq=GFFlist[chromosome][intronend:start].seq+splicedseq
			else:
				splicedseq+=(GFFlist[chromosome][end:intronstart]).seq
				splicedcounts+=counts[strand][chromosome][end:intronstart]
		

	# Get 5'UTR length
	if manualutr5==1:		# We are setting a fixed UTR length
		UTR5len=bp[0]
	else:	# We are adding on annotated UTR and further extending if needed.
		if featid in utrtable5:
			UTR5len=utrtable5[featid][1]-utrtable5[featid][0]+bp[0]
		else:
			UTR5len=0+bp[0]	# No annotation, no UTR appended. 
			
	if strand==0:
		if exon2len==-1:	# NO UTR5 intron
			splicedseq=(GFFlist[chromosome][start-UTR5len:start]).seq+splicedseq
			splicedcounts=counts[strand][chromosome][start-UTR5len:start]+splicedcounts
		else:	# DEal with intron - extend off of it.
			if UTR5len>=exon2len:
				if manualutr5==0 and UTR5len>bp[0]:	# If annotation there, get rid of intron for now since it is in. But don't if manual or if no annotation exists.
					UTR5len-=(intronlen)		# Have to do this for non-manual since input length from Cherry lab GFF includes introns.
				splicedseq=(GFFlist[chromosome][intronstart-(UTR5len-exon2len):intronstart]).seq+splicedseq
				splicedcounts=counts[strand][chromosome][intronstart-(UTR5len-exon2len):intronstart]+splicedcounts
			#else:	# else case is either Manual case or 5'UTR intron without annotation - handled below after shifting.

	else:
		if exon2len==-1:	# NO UTR5 intron
			splicedseq+=(GFFlist[chromosome][end:UTR5len+end]).seq
			splicedcounts+=counts[strand][chromosome][end:UTR5len+end]
		else:	# DEal with intron - extend off of it.
			if UTR5len>=exon2len:
				if manualutr5==0 and UTR5len>bp[0]:				# If annotation there, get rid of intron for now since it is in. But don't if manual or if no annotation exists.
					UTR5len-=(intronlen)		# Have to do this for non-manual since input length from Cherry lab GFF includes introns.
				splicedseq+=(GFFlist[chromosome][intronend:intronend+(UTR5len-exon2len)]).seq
				splicedcounts+=counts[strand][chromosome][intronend:intronend+(UTR5len-exon2len)]
		#	else:	# else case is either Manual case or 5'UTR intron without annotation - handled below after shifting.

		
	
	# Get the actual ends of the feature for neighborchecking below.
	if strand==0:
		if UTR5len>exon2len and exon2len!=-1:	# Use simple greater than since if equal, then it's just exon2.
			revisedstartoffeature=intronstart-(UTR5len-exon2len)	# With 5'UTR intron or manual extended more than exon2
		else:
			revisedstartoffeature=start-UTR5len			# Case of no 5'UTR intron or manual set to be less than or equal to exon 2 (or no annotation).
			if exon2len!=-1:
				revisedstartoffeature-=(exon2len-UTR5len+intronlen)	# Put this at start of intron so shift below can work.
	else:
		if UTR5len>exon2len and exon2len!=-1:	# Use simple greater than since if equal, then it's just exon2.
			revisedendoffeature=intronend+(UTR5len-exon2len)	# With 5'UTR intron or manual extended more than exon2
		else:
			revisedendoffeature=end+UTR5len				# Case of no 5'UTR intron or manual set to be less than or equal to exon 2 (or no annotation).
			if exon2len!=-1:
				revisedendoffeature+=(exon2len-UTR5len+intronlen)	# Put this at start of intron so shift below can work.
			
	
	# Add 3'UTRs on and revise ends of feature for neighborchecking below.
	if manualutr3==1:
		UTR3len=bp[1]
	else:
		# Add 3'UTR length.					
		if featid in utrtable3:
			UTR3len=utrtable3[featid][1]-utrtable3[featid][0]+bp[1]
		else:
			UTR3len=0+bp[1]
			
	if strand==0:
		splicedcounts+=counts[strand][chromosome][end:end+UTR3len]
		splicedseq+=GFFlist[chromosome][end:end+UTR3len].seq
		revisedendoffeature=end+UTR3len
	else:
		splicedcounts=counts[strand][chromosome][start-UTR3len:start]+splicedcounts
		splicedseq=GFFlist[chromosome][start-UTR3len:start].seq+splicedseq
		revisedstartoffeature=start-UTR3len
	
	
	
	# check overlap. Overlapping genes going in different directions will not be considered overlapping. Shift doesn't figure in since all shifted the same.
	# It will get the UTR information for current feature and adjust as needed, and then compare to annotated UTR of neighbors. If no UTR, use gene ends.
	
	if goodgenes == 1:
		feats=neighbors(GFFlist,chromosome,feature)	# This now gives the next feature on the same strand. Update July 2020 (used to be either strand).
		prevfeat=feats[0]
		nextfeat=feats[1]
		prevend=GFFlist[chromosome].features[prevfeat].location.end.position
		nextstart=GFFlist[chromosome].features[nextfeat].location.start.position
		
		prevdir=GFFlist[chromosome].features[prevfeat].strand==GFFlist[chromosome].features[feature].strand
		nextdir=GFFlist[chromosome].features[nextfeat].strand==GFFlist[chromosome].features[feature].strand

		# July 2020, this is updated to be the absolute position, not the distance, so as to deal with 5'UTR introns.
		# UTR length of interest is used and next neighbor on strand considered - so nearest on strand neighbor always considered, not just nearest neighbor on any strand.
		
		prevfeatid=GFFlist[chromosome].features[prevfeat].id
		nextfeatid=GFFlist[chromosome].features[nextfeat].id
		
		# Deal with case where the feature of interest is first or last on the chromosome.
		if prevfeatid[0:3]=="chr":
			prevdir=True
			prevend=0
		if nextfeatid[0:3]=="chr":
			nextdir=True
			nextstart=10000000000000000000000000000000000
		
		
		if (prevdir == False or nextdir == False):	# chromosomes okay.		
			print("Error, neighboring features on wrong strands.") 
			print("feature "+featid+" prevdir "+str(prevfeatid)+" "+str(prevdir)+" nextdir "+str(nextfeatid)+" "+str(nextdir))
			exit()
		
		if strand==0:
			if prevfeatid in utrtable3:
				prevtrueend=utrtable3[prevfeatid][1]
			else:
				prevtrueend=prevend
			if nextfeatid in utrtable5:
				nexttrueend=utrtable5[nextfeatid][0]
			else:
				nexttrueend=nextstart
		else:
			if prevfeatid in utrtable5:
				prevtrueend=utrtable5[prevfeatid][1]
			else:
				prevtrueend=prevend
			if nextfeatid in utrtable3:
				nexttrueend=utrtable3[nextfeatid][0]
			else:
				nexttrueend=nextstart
		
		if revisedstartoffeature<prevtrueend or revisedendoffeature>nexttrueend:	
			return[-2,-2,-2,-2,-2,-2,-2]
		
		
	# Perform shift for ribosome counts. Can assume all genes have already included any 5'UTR introns. For manual with UTR5 shorter than exon2, exception is made.
	if bp[2]!=0: # Check if there is a shift. Then do it.
		if strand==0:
			if bp[2]>0:
				splicedcounts=counts[strand][chromosome][revisedstartoffeature-bp[2]:revisedstartoffeature]+splicedcounts[:-bp[2]]
			else:
				splicedcounts=splicedcounts[-bp[2]:]+counts[strand][chromosome][revisedendoffeature:revisedendoffeature-bp[2]]
		else:
			if bp[2]>0:		
				splicedcounts=splicedcounts[bp[2]:]+counts[strand][chromosome][revisedendoffeature:revisedendoffeature+bp[2]]
			else:
				splicedcounts=counts[strand][chromosome][revisedstartoffeature+bp[2]:revisedstartoffeature]+splicedcounts[:bp[2]]

	
	# RC for -1 strand.
	if strand==1:
		splicedcounts.reverse()
		splicedseq=splicedseq.reverse_complement()
		
	if UTR5len<exon2len:		# Case of manual utr5 being less than exon2 or non-manual without any annotation - shrinkage needed.
		# Shrink to requested length of 5'UTR - manual case only allows this.
		splicedseq=splicedseq[exon2len-bp[0]:]
		splicedcounts=splicedcounts[exon2len-bp[0]:]		
		
# Final check for other errors, or running off the end of a chromosome.
	if len(splicedcounts)!=(CDSlen+UTR5len+UTR3len) or len(splicedseq)!=(CDSlen+UTR5len+UTR3len):
		print("Somewhere givegene ran off the end of a chromosome or shift was too big. Output length not equal to full transcript.")
		print("Feature = "+featid+" CDS="+str(CDSlen)+" UTR5="+str(UTR5len)+" UTR3="+str(UTR3len)+" nt")
		return [-1,-1,-1,-1,-1,-1,-1]
	
	
	# Get 1st intron that is fully in CDS for special output mode.
	gotintron=0
	for item in GFFlist[chromosome].features[feature].sub_features:
		if item.type == 'intron':
			gotintron=1
			start_feat = int(item.location.start.position)
			end_feat = int(item.location.end.position)
			intron1seq=(GFFlist[chromosome][start_feat:end_feat]).seq
			intron1counts=counts[strand][chromosome][start_feat:end_feat]
			if strand==0:
				break		# Should always give first coding intron in gene.
	if gotintron==0:
		return [splicedcounts,splicedseq,-1,-1,-1,UTR5len,UTR3len]
			
	if bp[2]!=0:
		if strand==0:
			if bp[2]>0:		# Note we need start_feat/end_feat in all these.
				intron1counts=counts[strand][chromosome][start_feat-bp[0]-bp[2]:start_feat-bp[0]]+intron1counts[:-bp[2]]
			else:
				intron1counts=intron1counts[-bp[2]:]+counts[strand][chromosome][end_feat+bp[1]:end_feat+bp[1]-bp[2]]
		else:
			if bp[2]>0:		
				intron1counts=intron1counts[bp[2]:]+counts[strand][chromosome][end_feat+bp[1]:end_feat+bp[1]+bp[2]]
			else:
				intron1counts=counts[strand][chromosome][start_feat-bp[0]+bp[2]:start_feat-bp[0]]+intron1counts[:bp[2]]
	# RC for -1 strand.
	if strand==1:
		intron1counts.reverse()
		intron1seq=intron1seq.reverse_complement()
	if strand==0:
		intronframe=(start_feat-start)%3
	else:
		intronframe=(end-end_feat)%3

	return [splicedcounts,splicedseq,intron1counts,intron1seq,intronframe,UTR5len,UTR3len]	


# Returns feature numbers of nearest genes that are nondubious.
# If no neighbor (i.e. the current feature is at the end of gene), returns feature 0, which is the whole chromosome.
# OK to return whole chromosome because chromosomes don't have a strand to strand check fails and first/end features are still reported.
# This is used by the program givegene.
# July 2020, now there is a check to make sure next feature is on same strand.
# Failure to find a neighhbor due to being near the end of a chromosome returns the entire chromosome.
def neighbors(GFFlist,chromosome,feature):
	if feature<0 or feature >=len(GFFlist[chromosome].features):
		print("Illegal feature given to neigbors program.")
		return[-1,-1]
	
	strand=GFFlist[chromosome].features[feature].strand
	
	i=1
	checkvar1=True
	checkvar2=True
	while(checkvar1==True or checkvar2==True):
		if((feature-i)<=0 and checkvar1==True):
			checkvar1=False
			prevfeat=0
		if((feature+i)>=len(GFFlist[chromosome].features) and checkvar2==True):
			checkvar2=False
			nextfeat=0
		if(checkvar1==True):
			# Need to do this has_key check here and below in case we're using a coli genome.
			if "orf_classification" in GFFlist[chromosome].features[feature-i].qualifiers:
				if (GFFlist[chromosome].features[feature-i].qualifiers["orf_classification"][0]!='Dubious'):
					
					dub=False
				else:
					dub=True
			else:
				dub=False
			if (GFFlist[chromosome].features[feature-i].type=='gene' and dub==False and strand==GFFlist[chromosome].features[feature-i].strand):
				prevfeat=feature-i
				checkvar1=False
		
		if(checkvar2==True):
			if "orf_classification" in GFFlist[chromosome].features[feature+i].qualifiers:
				if (GFFlist[chromosome].features[feature+i].qualifiers["orf_classification"][0]!='Dubious'):
					
					dub=False
				else:
					dub=True
			else:
				dub=False			
			if (GFFlist[chromosome].features[feature+i].type=='gene' and dub==False and strand==GFFlist[chromosome].features[feature+i].strand):
				nextfeat=feature+i
				checkvar2=False
		i+=1
		
	return [prevfeat,nextfeat]	   


	


# Converts GFFlist to GFFgen
def GFFlisttogen(GFFlist):
	for chrom in GFFlist:
		yield GFFlist[chrom]


