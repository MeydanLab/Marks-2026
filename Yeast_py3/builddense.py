# This file holds the programs used in the part of the workflow involved in putting ribosome density into data structures.
# Done after bowtie.
# Not every program here is routinely used.

import re
from BCBio import GFF
from Bio import Seq
from Bio import SeqIO
import csv
import seqtools
import gentools
import copy

########################## Workflow for setdense ################## 

# GFFgen_filename is a string of the filename with full path
# splices_filename is the fasta file string with full path used to make the library for bowtie to align splices, 
# juncrange is the size used to make splice junctions (generally readlength -1). This is good to use even if mismatches are more than 0.
# pathandfilestring is the path and first letters where the chromosome binary arrays will be written.
# smallestsize is the smallest readlength accepted, and largestsize the biggest.
# samgenmain_filename is the mapped reads for the genome, samgensplices_filename is for splice junctions.
# The totreads parameter was added on 2/3/12 so that one could easily use the script to generate density and wigfiles from the splice only data. Also to other datasets.
# In any alignment, simply put in the value of totreads computed previously (reads mapping to genome - ie sum of main and jcns). To not do an alignment, put in -1 for the particular samgen_filename.
# If you want to go back to having it compute the # of reads and using that for normalizing, put -1 in for totreads.

#### Starting Aug 2020. This will be called setdense_wf2 to avoid confusion with the old function.
#### No longer doing splice clean. 

# assignmode can be: end5,end3,cov,covflip
# species can be: scer,pombe,coli

def setdense_wf2(GFFgen_filename,samgenmain_filename,samgensplices_filename,juncrange,pathandfilestring,wigpathandfile,smallestsize,largestsize,species,assignmode,totreads):
	
	#First, just a bunch of checks for illegal input and print statements to double verify the inputs.
	if juncrange<0:
		print("Error in junc range. Cannot be negative.")
		exit()
		
	if species!="coli" and (smallestsize<0 or largestsize<=0):
		print("Error in sizes - negative was inputted.")
		exit()
	
	if (assignmode=="cov" or assignmode=="covflip") and species!="scer":
		print("Coverage only works with S. cerevisiae.")
		exit()
	
	if species=="coli":
		print("Species = coli. Assign mode is custom and assign mode input ignored.")
		if smallestsize==0 or largestsize==0:
			print("At this time, coli cannot handle 0 for smallest or largest size.")	# This is a special coli mapping strategy from Li et al.
		smallestsize*=-1
		largestsize*=-1
	elif species=="pombe":
		print("Species = pombe.")
	elif species=="scer":
		print("Species = S cerevisiae.")
	else:
		print("Illegal species.")
		exit()
	
	if species!="coli":
		if assignmode=="end5":
			print("Assignment will be 5' end.")
		elif assignmode=="end3":
			print("Assignment will be 3' end.")
		elif assignmode=="cov":
			print("Assignment will be coverage.")
		elif assignmode=="covflip":
			print("Assignment will be coverage. Strand will be flipped.")
		else:
			print("Illegal assign mode.")
			exit()		

	dens=0
	dens2=0
	outputdata1={}
	outputdata2={}
	GFFgen=GFF.parse(GFFgen_filename)
	
	if samgenmain_filename!=-1:
		f_samgenmain=open(samgenmain_filename)
		samgenmain=csv.reader(f_samgenmain,delimiter='	')
		
	if samgensplices_filename!=-1:
		f_samgensplices=open(samgensplices_filename)
		samgensplices=csv.reader(f_samgensplices,delimiter='	')
	 
	 # Create BioSeq files from the GFF file for modification of the letter annotations field.
	for sequence in GFFgen:
		outputdata1[sequence.id]=[0 for x in range(len(sequence))]
		outputdata2[sequence.id]=[0 for x in range(len(sequence))]
	
	if samgenmain_filename!=-1:
		dens=setdense3(outputdata1,outputdata2,samgenmain,smallestsize,largestsize,assignmode)	
		#print str(dens)+" chromosome reads mapped."
		
	if samgensplices_filename!=-1:
		if (species=="scer"):
			# For cerevisiae
			dens2=setdensesplice3(samgensplices,outputdata1,outputdata2,juncrange,smallestsize,largestsize,assignmode)
			#print str(dens2)+" splice reads mapped."
		elif species=="pombe":
			#For pombe
			GFFgen=GFF.parse(GFFgen_filename)
			GFFtable=seqtools.makeidtable(GFFgen)
			dens2=setdensesplice4(GFFtable,samgensplices,outputdata1,outputdata2,juncrange,smallestsize,largestsize,assignmode)
		else:
			print("Problem with species.")
			exit()

	if totreads==-1:
		totreads=dens+dens2
		print("total reads mapped: "+str(totreads))
		print("Using sum of chromosome and splice mappings for normalization.")
	else:
		print(str(totreads)+" total reads entered. Using this for normalization.")
	
	#Now we normalize data.
	norm_m(outputdata1,(totreads))
	norm_m(outputdata2,(totreads))
	
	# Write to file as float. Both density and wig.
	writecountsf(outputdata1,pathandfilestring+"_plus_")
	countstowig(outputdata1,wigpathandfile+"_plus")
	writecountsf(outputdata2,pathandfilestring+"_minus_")
	countstowig(outputdata2,wigpathandfile+"_minus")
	
	if samgenmain_filename!=-1:
		f_samgenmain.close()
	if samgensplices_filename!=-1:
		f_samgensplices.close()

# Divide every read value by number of mapped reads. Values then in rpm.
# No need to return readcounts since it is a mutable type being a list.
def norm_m(readcounts,reads):
    for chrom in list(readcounts.keys()):
        for position in range(len(readcounts[chrom])):
            readcounts[chrom][position]/=float(reads)
            readcounts[chrom][position]*=1E6            
            
# Function to display individual gene counts.
def countstowig(countsfile,filestring):
	import random
	f=open(filestring+".wig","w")
	filestring=filestring.partition("_")[0][-3:]
	random.seed(filestring)
	c1=random.randint(0,255)
	random.seed(filestring+"xxx")
	c2=random.randint(0,255)
	random.seed(filestring+"000")
	c3=random.randint(0,255)
	f.write("track name=tracklabel viewLimits=-5:5 color="+str(c1)+','+str(c2)+','+str(c3)+"\n")
	for chrom in list(countsfile.keys()):
		if chrom[0:3]=='chr':		# Case for S cer
			f.write("fixedStep  chrom="+chrom+"  start=1  step=1\n")
		elif chrom[0:1]=='I' or chrom[0:1]=='M' or chrom[0:1]=='A':	# Case for pombe
			f.write("fixedStep  chrom="+chrom+"  start=1  step=1\n")
		else:	# Case for e coli
			f.write("fixedStep  chrom=\""+chrom+"\"  start=1  step=1\n")	# Note the extra '' around chrom added for E coli. Should be fine for yeast. Actually it's not. So new code here to detect that.
		for position in countsfile[chrom]:
			f.write(str(position)+"\n")
	f.close()


# Function to convert wigfiles back to a countsfile.
# This is pretty basic and handles 2 types of wigs (fixed and variable) that I have created in the past. 
# Fixedstep assumed to be stepsize 1 and start value to be 1. Can be changed later by adding more code.
# chrlen is a dictionary of chromosomes that is provided for variablestep wigs. It specifies the length of each chrom.
def wigtocounts(wigfile,chrlen):
	resultlist={}
	f=open(wigfile)
	chrom=""
	steptype=-1			# Fixed is 1, Variable is 0, and not assigned is -1.
	for line in f:		# Should be okay to be checking past ends of lines. Sloppy but no problems yet.
		if line[0:9]=='fixedStep':	
			steptype=1
			if line[11:17]=="chrom=":
				chrom=gentools.parsenext(line[17:],'  ')
				chrom=cleanchrom(chrom)  # Remove any end of lines or quotations.
				resultlist[chrom]=[]
			else:
				print("Error - no chrom.")
				return -1
			
		elif line[0:12]=='variableStep':
			steptype=0
			if line[14:20]=="chrom=":
				chrom=gentools.parsenext(line[20:],'  ')
				chrom=cleanchrom(chrom)	# Remove any end of lines or quotations.
				resultlist[chrom]=[float(0) for x in range(chrlen)]
			else:
				print("Error - no chrom.")
				return -1
			
		else:
			if steptype==1:
				resultlist[chrom].append(float(line))
			if steptype==0:
				
				resultlist[chrom][int(gentools.parsenext(line,'	'))-1]=float(gentools.parselast(line,'	'))
			else:
				continue	# This happens on early comment lines.
				#print "Step type not specified. Error."
				#return -1

	f.close()
	return resultlist

# Support function for converting wig to countslist function above.
def cleanchrom(chrom):
	if chrom[-1]=="\n":
		chrom=chrom.strip()
	if chrom[0]=='\"':
		chrom=chrom[1:-1]
	return chrom

# Function to write counts files to disk - Use for float. For int, switch "f" to "i".
def writecountsf(resultlist,filestring):  #Resultlist it is the list of counts, filestring is the file prefix for each chr to be written.
    import struct
    f2=open(filestring+"keys","w")
    for chrom in list(resultlist.keys()):
        f=open(filestring+chrom,"wb")
        for position in resultlist[chrom]:
            f.write(struct.pack("f",position))
        f.close()   
        f2.write(chrom+"\n")
    f2.close()    

# Function to read in counts files. Use for floats. If your wig is integers, then you can change the "f" below to "i."
def readcountsf(filestring):
    import struct
    keys=[]
    resultlist={}
    f2=open(filestring+"keys","r")
    for line in f2:
        keys.append(line.rstrip('\n'))  # rstrip takes off the end of line character.
    for chrom in keys:
        resultlist[chrom]=[]
        with open(filestring+chrom,"rb") as f:
            nextval = f.read(4)
            while nextval != b'': 
                resultlist[chrom].append(struct.unpack("f",nextval)[0])    # This returns a tuple, not an int, so need to take 1st item. If your wig is integer, you can change this to "i." Old function called readcounts() had this.
                nextval= f.read(4)
    f2.close()
    return resultlist
    # Note that there is no file closing, probably not ideal, but everything seems to work fine.
    

# What this function does is to create a big binary array for each chromosome and returns it. 
# This array can then be written to file or converted to a wig file.
# It requires a SAM generator.
# setdensesplice is a sister function that puts in splice jcn reads.
# f_splice is a file handle to a fasta of splice junctions. splicerange is bp on sides of splice, typically 99. 
# smallest is smallest allowed readlength, largest the biggest.

# The inputput outputdata is the dict of density values. It is not returned since it is a mutable type.
# In the splice version, this is called readcounts.
# The value returned is the total number of mapped reads.

def setdense3(outputdata1,outputdata2,samgen,smallest,largest,assignmode):
	if assignmode=="cov":
		coverage=1
		flip=0
	elif assignmode=="covflip":
		coverage=1
		flip=1
	elif assignmode=="end3":
		end3=1
		coverage=0
	elif assignmode=="end5":
		end3=0
		coverage=0
	else:
		print("Illegal input")
		exit()
	
	# Do for ecoli
	if smallest<0 and largest<0:
		smallest*=-1
		largest*=-1
		trimlen=int(round(smallest/float(2)))-1
		print("Density assigned for all positions >"+str(trimlen)+"bp from end.")
		coli=True
		print("Doing e coli actually.")
	else:
		coli=False

	mappedreads=0
	mappedreadslocal=0
# Loop through the samfile.
	for read in samgen:
		if read[0][0] == '@':   # Ignore header lines.
			continue
        
		if read[1] == '4':      # A bowtie mismatch. These are the same as those previously saved in the files bowtie output as reads without matches and reads exceeding the -m limit. No need to write as fastq.
			continue
        
		chrom = read[2]             # chromosome identified for read in bowtie
		readid = read[0]            # read id
		startp = int(read[3]) - 1    # start position. Need to subtract 1 since genomic sequence starts at 1, but list structures for read summing start at 0.
		seq = Seq.Seq(read[9])      # sequence of the read
		length = len(seq)           # length of read

		mappedreads+=1
		
		# Note that Bowtie reverse complements any sequence aligning to the reverse strand.  
		
		# Filter to get rid of reads of particular length. Or a particular strand.
		if (length<smallest or length > largest):
			continue
		mappedreadslocal+=1
		
		if coli==False:

			if coverage==1:
			# Handle Coverage
				for posinread in range(startp,startp+length):
					if (read[1] == '0' and flip==0) or (read[1] == '16' and flip==1):
						outputdata1[chrom][posinread]+=(float(1)/float(length))
					else:
						outputdata2[chrom][posinread]+=(float(1)/float(length))

			else:
			
				if end3==0:
					if (read[1] == '16'):
						start=startp+length-1   
						outputdata2[chrom][start]+=1 
					if (read[1] == '0'):
						outputdata1[chrom][startp]+=1 
				else:
					if (read[1] == '16'):
						outputdata2[chrom][startp]+=1 
					if (read[1] == '0'):
						start=startp+length-1   
						outputdata1[chrom][start]+=1 
		
		else:
			efflength=length-2*trimlen
			if (read[1] == '16'):
				start=startp+length-1  
				for trimpos in range(efflength):
					outputdata2[chrom][start-trimlen-trimpos]+=(float(1)/efflength)
				
			if (read[1] == '0'):
				for trimpos in range(efflength):
					outputdata1[chrom][startp+trimlen+trimpos]+=(float(1)/efflength) 
		               
	print(str(mappedreadslocal)+" reads mapped to genome.")
	return(mappedreads)


def setdensesplice3(samgen,readcounts1,readcounts2,splicerange,smallest,largest,assignmode):
	if assignmode=="cov":
		coverage=1
		flip=0
		end3=-1
	elif assignmode=="covflip":
		coverage=1
		flip=1
		end3=-1
	elif assignmode=="end3":
		end3=1
		coverage=0
	elif assignmode=="end5":
		end3=0
		coverage=0
	else:
		print("Illegal input")
		exit()
		
	mappedreads=0
	mappedreadslocal=0
           # loop through each read from sam file.
	for read in samgen:
		if read[0][0] == '@':   # Ignore header lines.
			continue

		if read[1] == '4':      # A bowtie mismatch or greater then -m limit. 
			continue

		junction = read[2]             # chromosome identified for read in bowtie
		readid = read[0]            # read id
		startj = int(read[3]) - 1    # start position. Need to subtract 1 since genomic sequence starts at 1, but list structures for read summing start at 0.
		seq = Seq.Seq(read[9])      # sequence of the read
		length = len(seq)           # length of read

		mappedreads+=1
		
		# Filter to get rid of reads of particular length. Or a particular strand.
		if (length<smallest or length > largest):
			continue
		mappedreadslocal+=1
        
		chrom=re.split("_",junction)[0]
		intronpos=int(re.split("_",junction)[1])
		intronlength=int(re.split("_",junction)[2])
        # intronpos is the genomic position of the start of the splicerange*2 bp splice junction. Recall this is a 0-based GFF parser coords. Intron start -junction range.
		# intronlength is length of actual intron in bp.
		
		if len(re.split("_",junction))>3:
			intron2pos=int(re.split("_",junction)[3])		# This is start of 2nd intron. No offset.
			intron2length=int(re.split("_",junction)[4])
			interveningexonlen=intron2pos-intronpos-intronlength-splicerange
		else:
			intron2pos=0			# Not used
			intron2length=0
			interveningexonlen=0

		
		if (read[1] == '16' and end3==0) or (read[1]=='0' and end3==1):
			start = startj+length-1
			modstart=start+intronlength+intronpos
			if start>=(splicerange+interveningexonlen):
				modstart+=intron2length		# Should never happen unless 2 introns.
		else:
			start = startj
			modstart=start+intronpos 
			if start>=splicerange:			# Should never happen unless 2 introns.
				modstart+=intronlength		
     
     	# Handle Coverage
		if coverage==0:	
			if read[1] == '0':
				readcounts1[chrom][modstart]+=1   
			if read[1] == '16':
				readcounts2[chrom][modstart]+=1
		else:
			for posinread in range(startj,startj+length):
				modstart=posinread+intronpos
				if posinread>=splicerange:			
					modstart+=intronlength
				if posinread>=(splicerange+interveningexonlen):
					modstart+=intron2length
				if (read[1] == '0' and flip==0) or (read[1] == '16' and flip==1):
					readcounts1[chrom][modstart]+=(float(1)/float(length))
				else:
					readcounts2[chrom][modstart]+=(float(1)/float(length))
	
	if (length-1)>splicerange:
		print("***************")
		print("WARNING. Read length is "+str(length))
		print("However, juncrange is "+str(splicerange))
		print("***************")
		exit()

	print(str(mappedreadslocal)+" reads mapped to junctions.")
	return mappedreads
           
# For pombe. This works in an entirely different way from S.cer. It takes the reads mapped to the cDNA fasta (plus artificial UTR length extensions), converts them to chromosomal, and fills in the density.
def setdensesplice4(GFFtable,samgen,readcounts1,readcounts2,splicerange,smallest,largest,assignmode):
	
	if assignmode=="end3":
		end3=1
	elif assignmode=="end5":
		end3=0
	else:
		print("Illegal input")
		exit()
		
	mappedreads=0
	mappedreadslocal=0
	genecountdict={} 
	genomecounts=[readcounts1,readcounts2]
           # loop through each read from sam file.
	for read in samgen:
		if read[0][0] == '@':   # Ignore header lines.
			continue

		if read[1] == '4':      # A bowtie mismatch or greater then -m limit. 
			continue

		geneheader = read[2]        # The header info pulled by bowtie.
		readid = read[0]            # read id
		startj = int(read[3]) - 1    # start position. Need to subtract 1 since genomic sequence starts at 1, but list structures for read summing start at 0.
		seq = Seq.Seq(read[9])      # sequence of the read
		length = len(seq)           # length of read

		mappedreads+=1
		
		gene=re.split("_",geneheader)[0]
		lengthofgene=int(re.split("_",geneheader)[1])
		if read[1] == '16':
			strand=1
		elif read[1] == '0':
			strand=0
		
		# Filter to get rid of reads of particular length. Or a particular strand.
		if (length<smallest or length > largest):
			continue
		mappedreadslocal+=1
        
		if end3==1:
			if read[1]=='0':
				startj+=(length-1)
			if read[1]=='16':
				startj-=(length-1)
					
		# Create invidual gene lists in a dictionary
		if gene in genecountdict:
			genecountdict[gene][strand][startj]+=1
		else:
			genecountdict[gene]=[[0 for x in range(lengthofgene)],[0 for x in range(lengthofgene)]]
			genecountdict[gene][strand][startj]+=1
			
				
	# Now go through list of genes.
	for genename in list(genecountdict.keys()):
		import seqtools_extra
		seqtools_extra.convertmrnatogenomic3(GFFtable,genename,genecountdict[genename],genomecounts,splicerange)
   
	print(str(mappedreadslocal)+" reads mapped to junctions.")
	return mappedreads


# This function returns a list of counts of reads mapping to a sequence of interest, and their sizes.
def setdense_adhocmsg(samgen,seqlength,smallest,largest):
	mappedreads=0
	mappedreadslocal=0
	# seqlength is the length of the sequence of interest
	countsp=[0 for x in range(seqlength)]
	countsm=[0 for x in range(seqlength)]
	# For now, we will make an assumption that readlengths are under 100. This can be changed in later code.
	lengthlim=100
	lengthsp=[0 for x in range(100)]
	lengthsm=[0 for x in range(100)]

	for read in samgen:
		if read[0][0] == '@':   # Ignore header lines.
			continue

		if read[1] == '4':      # A bowtie mismatch or greater then -m limit. 
			continue

		#junction = read[2]             # chromosome identified for read in bowtie
		readid = read[0]            # read id
		startj = int(read[3]) - 1    # start position. Need to subtract 1 since genomic sequence starts at 1, but list structures for read summing start at 0.
		seq = Seq.Seq(read[9])      # sequence of the read
		length = len(seq)           # length of read
		if length>=lengthlim or length<0:
			print("Error: readlength exceeds hardcoded limit of "+str(lengthlim))
			return [0,0,0,0]		
		
		mappedreads+=1
		
		# Filter to get rid of reads of particular length. Or a particular strand.
		if (length<smallest or length > largest):
			continue
		mappedreadslocal+=1
		
		if read[1] == '16':
			start = startj+length-1
			countsm[start]+=1
			lengthsm[length]+=1
		else:
			start = startj
			countsp[start]+=1
			lengthsp[length]+=1
		
	print(str(mappedreadslocal)+" reads mapped to sequence.")
	return [countsp,countsm,lengthsp,lengthsm]
	

# This function will read in samfiles for a number of different viral sequences.
# It then outputs the densities, which are then written to file.
# Lengths of positive and negative mapped reads are written together as one file (not interlaced).
def viraldense_wf(samfilelist,lengthlist,normalization,outfilelist,smallest,largest):

	if len(samfilelist)!=len(lengthlist) or len(samfilelist)!=len(outfilelist):
		print("Error on lengths of input lists.")
	
	for i in range(len(samfilelist)):
		# Make samgen
		samfile=open(samfilelist[i])
		samgen=csv.reader(samfile,delimiter='	')
		# Create density
		density=setdense_adhocmsg(samgen,lengthlist[i],smallest,largest)
		mappedcounts=sum(density[0])+sum(density[1])
		
		if normalization!=-1:
		
			# Normalize
			for position in range(lengthlist[i]):
				density[0][position]/=float(normalization+mappedcounts)
				density[0][position]*=1E6
				density[1][position]/=float(normalization+mappedcounts)
				density[1][position]*=1E6
		# Do nothing if normalzation is -1.
		
		# Write out density
		density[1]=[-1*x for x in density[1]]	# Negate neg strand values.
		
		print(outfilelist[i])
		fp=open(outfilelist[i]+"_plus","wb")
		fm=open(outfilelist[i]+"_minus","wb")
		gentools.writelistbinint(density[0],fp)
		gentools.writelistbinint(density[1],fm)
		fsizes=open(outfilelist[i]+"_sizes","wb")
		gentools.writelistbinint(density[2]+density[3],fsizes)
		fsizes.close()
		fp.close()
		fm.close()


# This function is generally not used anymore.
# This function creates an averaged density file from a list of input files. 
# Files are opened one by one to save space and the average is returned. 
# There is a wf function below that calls it and a script to call it.
# Also, a switch for averaging vs combining is added. Set combine to 1 to truly combine, not average.
# Added on 7/25/14: the ability to combine according to the total number of reads in each file so that depth relevance is preserved.
# To do this, a new list called totcountslist is added. To simply combine all equally, just make all the values the same.
# Also included a shift value to shift a given density with respect to the original. Shift values must be positive. This isn't useful yet because shifts will be different for +/- genes.

def combinedensities(filelist,combine,totcountslist,shiftlist):
	listlen=len(filelist)
	
	# Open first file.
	counts1p=readcountsf(filelist[0]+"_plus_")
	counts1m=readcountsf(filelist[0]+"_minus_")
	for chrom in list(counts1p.keys()):
		chromlen=len(counts1p[chrom])
		for i in range(chromlen):
			counts1p[chrom][i]=0
			counts1m[chrom][i]=0
	
	counts1=[counts1p,counts1m]
	
	totcounts=float(sum(totcountslist))
	
	# Make average.
	for (ff,countsval,shift) in zip(filelist,totcountslist,shiftlist):
		counts2p=readcountsf(ff+"_plus_")
		counts2m=readcountsf(ff+"_minus_")
		for chrom in list(counts1[0].keys()):
			chromlen=len(counts1[0][chrom])
			
			for i in range(chromlen-shift):
				counts1[0][chrom][i+shift]+=(counts2p[chrom][i]*(countsval/totcounts))
				counts1[1][chrom][i+shift]+=(counts2m[chrom][i]*(countsval/totcounts))
				
	if combine==1:
		for chrom in list(counts1[0].keys()):
			chromlen=len(counts1[0][chrom])
			for i in range(chromlen):
				counts1[0][chrom][i]*=listlen
				counts1[1][chrom][i]*=listlen
	
	return counts1


def combinedensities_wf(combine,filelist,totcountslist,shiftlist,outpathandfilestring,wigpathandfile):
	dens=combinedensities(filelist,combine,totcountslist,shiftlist)
	# Write to file as float. Both density and wig.
	writecountsf(dens[0],outpathandfilestring+"_plus_")
	countstowig(dens[0],wigpathandfile+"_plus")
	writecountsf(dens[1],outpathandfilestring+"_minus_")
	countstowig(dens[1],wigpathandfile+"_minus")

    
# Very simple counter
# Takes a bowtie-created samfile and adds up all reads aligning to each entry
# Use for counting which tRNAs are getting reads. 
# Developed late 2019 for use with 40S profiling. Use wit quant_trna_wf_run.py
def counttrna(samfile,trnafasta,outputfileprefix):

	f_samgenmain=open(samfile)
	samgen=csv.reader(f_samgenmain,delimiter='	')

	trnadictionary={}
	
	fastagen=SeqIO.parse(trnafasta,"fasta")
	
	for geneentry in fastagen:
		aa=geneentry.id[1]
		anticodon=geneentry.id[3:6]
		startSGID=str(geneentry.description).find(":")+1
		endSGID=str(geneentry.description).find(",")
		SGDID=(geneentry.description)[startSGID:endSGID]
		trnadictionary[geneentry.id]=[anticodon,aa,SGDID,0,geneentry.description]
		
	
	for read in samgen:
		if read[0][0] == '@':   # Ignore header lines.
			continue
        
		if read[1] == '4':      # A bowtie mismatch. 
			continue
        
		feature = read[2]             

		if feature in trnadictionary:
			trnadictionary[feature][3]+=1
		else:
			print("Bad fasta - missing genes")
			exit()
		
	# Write Out
	writer = csv.writer(open(outputfileprefix+".csv", "wb"),delimiter=',')
	writer.writerow(["anticodon","aa","SGDID","count","fullinfo"])
	for row in list(trnadictionary.keys()):
		writer.writerow(trnadictionary[row])
	
	
# This was originally developed in late 2019 as a way to determine where 40S reads were mapping in the transcriptome.
# It is used with the quant_regions_wf_run.py script.			
# Counter for reads mapping to regions of the transcriptome.
# Takes many samfiles as inputs
# Adds up reads around starts/stops and UTRs/ORF.
# position settings determines where read ends have to map to be counted. It is distance from the first base of the P site. 4,4 would correspond to the Preiss lab approach. And 
# a value of 4, means they have to be at least 4 nt past the P site (ie the way Preiss lab did it). -3,-5 would correspond to the way we defined it at one point.
# A value of -3 means at least 1 base in an active (E,P,or A site).
def countregions(position5,position3,sizerangelower,sizerangeupper,samfiles,outputfileprefix):
	regiondictionary={}
	regiondictionary["headers"]=["UTR5","start","ORF","stop","UTR3"]
	regionlengths={}
	regionlengths["headers"]=list(range(sizerangelower,sizerangeupper+1))

	for samfile in samfiles:
		f_samgenmain=open(samfile)
		samgen=csv.reader(f_samgenmain,delimiter='	')
		
		sampleend=samfile.rfind(".")
		samplestart=samfile.rfind("/")
		sample=samfile[samplestart+1:sampleend]
		regiondictionary[sample]=[0,0,0,0,0]

		regionlengths["UTR5_"+sample]=[0 for x in range(sizerangeupper-sizerangelower+1)]
		regionlengths["start_"+sample]=[0 for x in range(sizerangeupper-sizerangelower+1)]
		regionlengths["ORF_"+sample]=[0 for x in range(sizerangeupper-sizerangelower+1)]
		regionlengths["stop_"+sample]=[0 for x in range(sizerangeupper-sizerangelower+1)]
		regionlengths["UTR3_"+sample]=[0 for x in range(sizerangeupper-sizerangelower+1)]
	
		for read in samgen:
			if read[0][0] == '@':   # Ignore header lines.
				continue
        
			if read[1] == '4':      # A bowtie mismatch. 
				continue
        
			feature = read[2]            
			end5 = int(read[3]) - 1    # start position. 
			seq = Seq.Seq(read[9])
			length = len(seq)           # length of read
			end3 = end5+length-1
   
			first_=feature.find("_")
			second_=feature.rfind("_")
   
			ORFstart=int(feature[first_+1:second_])
			ORFend=int(feature[second_+1:]) #This is start of 3'UTR

			length-=15
			
			# Decide where read maps:
			if end3<(ORFstart+position5):
				regiondictionary[sample][0]+=1	#UTR5
				regionlengths["UTR5_"+sample][length]+=1
			elif end5>(ORFend-6-position3):
				regiondictionary[sample][4]+=1	#UTR3
				regionlengths["UTR3_"+sample][length]+=1
			elif end5>(ORFstart-position3) and end3<(ORFend-6+position5):
				regiondictionary[sample][2]+=1	#ORF
				regionlengths["ORF_"+sample][length]+=1
			elif end5<=(ORFstart-position3) and end3>=(ORFstart+position5):
				regiondictionary[sample][1]+=1	#start
				regionlengths["start_"+sample][length]+=1
			elif end5<=(ORFend-6-position3) and end3>=(ORFend-6+position5):
				regiondictionary[sample][3]+=1	#stop
				regionlengths["stop_"+sample][length]+=1
			else:
				print("Should never get here")
				print(end5)
				print(end3)
				print(readid)				

		
	regionlengthsraw = copy.deepcopy(regionlengths)		
	
	for samfile in samfiles:	
		sampleend=samfile.rfind(".")
		samplestart=samfile.rfind("/")
		sample=samfile[samplestart+1:sampleend]
		
		# Normalize lengths
		sumlengths=sum(regionlengths["start_"+sample])
		for j in range(sizerangeupper-sizerangelower+1):
			regionlengths["start_"+sample][j]/=float(sumlengths)
		sumlengths=sum(regionlengths["stop_"+sample])
		for j in range(sizerangeupper-sizerangelower+1):
			regionlengths["stop_"+sample][j]/=float(sumlengths)
		sumlengths=sum(regionlengths["UTR5_"+sample])
		for j in range(sizerangeupper-sizerangelower+1):
			regionlengths["UTR5_"+sample][j]/=float(sumlengths)
		sumlengths=sum(regionlengths["UTR3_"+sample])
		for j in range(sizerangeupper-sizerangelower+1):
			regionlengths["UTR3_"+sample][j]/=float(sumlengths)
		sumlengths=sum(regionlengths["ORF_"+sample])
		for j in range(sizerangeupper-sizerangelower+1):
			regionlengths["ORF_"+sample][j]/=float(sumlengths)

	# Write Out
	gentools.writedicttoexcel(regiondictionary,outputfileprefix+"_counts_")
	gentools.writedicttoexcel(regionlengths,outputfileprefix+"_lengths_")
	gentools.writedicttoexcel(regionlengthsraw,outputfileprefix+"_lengths_raw_")
	gentools.transposecsv(outputfileprefix+"_lengths_")	
	gentools.transposecsv(outputfileprefix+"_counts_")	
	gentools.transposecsv(outputfileprefix+"_lengths_raw_")	
