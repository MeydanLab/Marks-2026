import gentools
import seqtools
from BCBio import GFF
from Bio import SeqRecord
from Bio import Seq
from Bio import SeqIO

# This is a new updated genometools as of August 2020. It eliminates all the old genometools except the finseqinstance scripts, which have been updated to work with the new givegene.
# It also has the makeyeasttranscriptome script for region analysis used to make piecharts in 40S paper analysis.
# It also retains scripts to do mapability analysis for repeats and low complexity regions.


# This function finds sequences in the yeast transcriptome
# This function is an update of the old version in the "extra" genometools. This handles all regions but eliminates the searches for charge, structure, and tAI. Use old functions if you want that.
# region - region of interest, 0 for 5'UTR, 1 for ORF, 2 for 3'UTR. Start of these regions is frame of reference.
# inframe is whether to look in frame for the sequence, 1 or 0. Set to 2 if motif is an aa sequence.
# inframe=2.1 means doing out of frame by +1 and 2.2 is out of frame by +2.
# motif is the motif of interest to look for
# If motif is "PENULT," then we get all triplets in 3'UTR that match the penultimate codon exactly.
# mismatches is number of allowed mismatches in motif.
# GFFlists are the main chromosome information and UTR information.
# findstop - Set this to 0 if not doing it. Otherwise, it will report the next stop codon downstream of the search motif that is <= this value. If negative, it looks > than the absolute value downstream.

def findseqinstance2(GFFlists,motif,inframe,region,mismatches,findstop):

	findstoprev=0
	if findstop<0:		# Case of filtering out short ORFs.
		findstop*=-1
		findstoprev=1
		
	if motif=="PENULT":
		penultmotif=1
	else:
		penultmotif=0
	GFFlist=GFFlists[0]
	hits=[]
	dummy={}
	hits.append(["headers","alias","chromosome","featurenum","gen_position","mrna_pos","Note"])
	
	frameadj=0
	if inframe>2 and inframe<3:
		if inframe==2.1:
			frameadj=1
		elif inframe==2.2:
			frameadj=2
		else:
			print("Frame error")
		inframe=2

	for sequence in GFFlist:
		dummy[sequence]=[0 for x in range(len(GFFlist[sequence]))]
	
	for chrom in GFFlist:
		feat_num=0
		for feature in GFFlist[chrom].features:
			localstops=[]			
			gg=seqtools.givegene(chrom,feat_num,GFFlists,[dummy,dummy],[0,0,0],2)	
			genesequence=gg[1]
			UTR5len=gg[5]
			UTR3len=gg[6]
			
			if genesequence == -1:		# For genes that are considered not allowed because they are dubious, nongene, wrong chrom, etc. Overlap ok. 
				feat_num+=1
				continue

			if region==0:
				if UTR5len>0:
					genesequence=gg[1][0:UTR5len]
				else:
					feat_num+=1
					continue
			elif region==1:
				genesequence=gg[1][UTR5len:len(gg[1])-UTR3len]
			
			elif region==2:
				if UTR3len>0:
					genesequence=gg[1][-UTR3len:]
				else:
					feat_num+=1
					continue
					
			if penultmotif==1:
				motif=str(gg[1][len(gg[1])-UTR3len-6:len(gg[1])-UTR3len-3])

			if inframe==2:
				genesequence=genesequence[frameadj:]
			
			if inframe > 1 and inframe<5:
				genesequence=genesequence.translate()
				
			genesequence=str(genesequence)	
			genelength=len(genesequence)
			if inframe==2 or inframe==1 or inframe==0:
				motiflength=len(motif)
			else:
				print("Error - inframe wrong.")
			
			
			i=0
				
			while(i<(genelength-motiflength+1)):
				# Determine if motif matches.
				if inframe<3 and inframe>-1:
					# Deal with any "_" present in test sequence.
					if(motif.count("_")>0):
						if mismatches>0:
							print("Can't have mismatches with wildcards.")
							exit()
						motifbreaks=motif.split("_")
						partialmotiflen=-1
						for partialmotif in motifbreaks:
							if partialmotiflen==-1:
								prevmotiflen=0		# First time through.
							else:
								prevmotiflen+=(partialmotiflen+1)
							partialmotiflen=len(partialmotif)
							if partialmotif!=genesequence[i+prevmotiflen:i+prevmotiflen+partialmotiflen]:
								testseq=0
								break
						else:		# Loop fell through without finding a mismatch.
							testseq=1
					
					elif motif == genesequence[i:i+motiflength]:					
						testseq=1
									
					elif mismatches>0:
						mm=gentools.mismatchcount(motif,genesequence[i:i+motiflength])
						if mm>mismatches:
							testseq=0
						else:
							testseq=1
					else:
						testseq=0
				else:
					print("Error with inframe.")			

				if findstop!=0 and testseq==1:			# Find the next stop in frame.
					motiflengthlocal=3	
					j=i
					j+=3
					testseq=0
					while(j<(genelength-motiflengthlocal+1)):
						if "TGA"==genesequence[j:j+motiflengthlocal] or "TAA"==genesequence[j:j+motiflengthlocal] or "TAG"==genesequence[j:j+motiflengthlocal]:
							testseq=1
							break
						j+=3
					if (testseq==1): # Keep a local list of stops found so not to duplicate and of certain length.
						for stopcodon in localstops:
							if stopcodon==j:
								testseq=0
						localstops.append(j)		# Note, this will duplicate (okay).
						
					if findstoprev==0:				# Eliminate long small ORFs. Short only.
						if ((j-i)>findstop):		
							testseq=0							
					else:							# Eliminate short small ORFs. Long only.
						if ((j-i)<=findstop):		
							testseq=0				
		# Add hit to hit list.
				if testseq==1:
					if inframe>1 and inframe<5:	# Convert back to bp code for aa search.
						pos=(i)*3+frameadj		# frameadj only will be nonzero if inframe is 2.1 or 2.2.
					else:
						if findstop==0:
							pos=i
						else:
							pos=j				
					
					if('Alias' in feature.qualifiers):
						alias=feature.qualifiers["Alias"]
					elif "Name" in feature.qualifiers:				
						alias=feature.qualifiers["Name"][0]
					else:
						alias="NoAlias"	
						
					if "Note" in feature.qualifiers:
						note=feature.qualifiers["Note"]
					else:
						note = "NA"				
					
					if region==0:
						posconv=pos-UTR5len
					elif region==2:
						posconv=pos+len(gg[1])-UTR5len-UTR3len
					else:
						posconv=pos
					genomic_pos=seqtools.convertmrnatogenomic(posconv,chrom,feat_num,GFFlist)
					hits.append([feature.id,alias,chrom,feat_num,genomic_pos,pos,note])
				
				if inframe==1 or inframe==5:	
					i+=3
				elif inframe == 0 or inframe==-1 or inframe == 2 or inframe==3 or inframe==4 or inframe==3.5:
					i+=1
				else:
					print("error")
					i+=3
			feat_num+=1
	print("Number of hits before filters = "+str(len(hits)-1))
	return hits


# Calling script that also processes the output.
# gffpathandfile - GFF file input
# utrgfffile5 - UTR5 annotation file
# utrgfffile3 - UTR3 annotation file
# motif - The motif to be searched. # If motif starts with a "/", the string is opened and used as the hit dict for the filtering tools in this wrapper. findseqinstance never actually called.
# inframe - see above, can be 0, 1, 2, 2.1, or 2.2
# outfilestring - output csv
# mismatches - Number of mismatches to tolerate in motif.
# noconsec - Sets how far apart consecutive motifs are allowed to be. 1 for nothing consecutive. Negative 1 for 1 per gene cases only.
# stopcodon - Throws out genes without this main orf stop codon. 0 for all stops.
# region - 0 for 5'UTR, 1 for main ORF, 2 for 3'UTR.
# regionwindow - Regions to exclude. Set to -1 to ignore.
# findstop - Set to nonzero value to find next downstream stop of motif and report that is <= nt away. Set to negative for > nt away.
def findseqinstance2_wf(gffpathandfile,utrgfffile5,utrgfffile3,motif,inframe,outfilestring,mismatches,noconsec,stopcodon,region,regionwindow,findstop):
	GFFgen=GFF.parse(gffpathandfile)
	GFFlist=seqtools.makeGFFlist(GFFgen)
	utrgff=GFF.parse(utrgfffile5)
	utrtable5=seqtools.makeutrtable(utrgff)
	utrgff=GFF.parse(utrgfffile3)
	utrtable3=seqtools.makeutrtable(utrgff)
	
	GFFlists=[GFFlist,utrtable5,utrtable3]
	
	# Pull in an existing list if required.
	if motif[0]=="/":
		hits=gentools.readexceltolist(motif)
	else:
		hits=findseqinstance2(GFFlists,motif,inframe,region,mismatches,findstop)

	
	
	if noconsec>0:
		goodhits=[]
		prevpos=-10
		prevgene="xxx"
		hitlen=list(range(len(hits)))
		for hit in hitlen:
			if hits[hit][0]=="headers":	# Take care of headers.
				goodhits.append(hits[hit])
				continue
			currenthit=hits[hit][0]
			currentpos=int(hits[hit][5])
			
			if currenthit!=prevgene or abs(currentpos-prevpos)>noconsec:
				goodhits.append(hits[hit])
			prevpos=currentpos
			prevgene=currenthit
		hits=goodhits
		print("Actual hits after consecutive removed = "+str(len(hits)-1)) 		# To get rid of headers. This wasn't included before so all hit #s are off by 1 prior to aam 3/3/13.
	
	# This is the case where you eliminate all genes with >1 instances.
	if noconsec<0:
		goodhits=[]
		hadhits={}
		hitlen=list(range(len(hits)))
		for hit in hitlen:
			if (hits[hit][0] in hadhits)==False:
				hadhits[hits[hit][0]]=1
			else:
				hadhits[hits[hit][0]]+=1
		for hit in hitlen:	
			if hadhits[hits[hit][0]]==1:
				goodhits.append(hits[hit])
		hits=goodhits
		print("Actual hits after consecutive removed = "+str(len(hits)-1)) 		# To get rid of headers. This wasn't included before so all hit #s are off by 1 prior to aam 3/3/13.
	
	
	if stopcodon!=0 and stopcodon!='0':
		dummy={}
		for sequence in GFFlist:
			dummy[sequence]=[0 for x in range(len(GFFlist[sequence]))]
		goodhits=[]
		hitlen=list(range(len(hits)))
		for hit in hitlen:
			if hits[hit][0]=="headers":	# Take care of headers.
				goodhits.append(hits[hit])
				continue
			chrom=hits[hit][2]
			feat_num=int(hits[hit][3])
			currentstop=givestopcodon(chrom,feat_num,GFFlists,dummy)
			if stopcodon==currentstop:
				goodhits.append(hits[hit])
		hits=goodhits
		print("Actual hits after bad stop codons removed = "+str(len(hits)-1))
	
	
	# A new bit of code 6/8/14 to remove hits that fall outside the desired regionwindow of gene.
	if regionwindow!=-1:
	
		dummy={}
		for sequence in GFFlist:
			dummy[sequence]=[0 for x in range(len(GFFlist[sequence]))]
		goodhits=[]
		hitlen=list(range(len(hits)))
		for hit in hitlen:
			if hits[hit][0]=="headers":	# Take care of headers.
				goodhits.append(hits[hit])
				continue
			chrom=hits[hit][2]
			feat_num=int(hits[hit][3])
			gg=seqtools.givegene(chrom,feat_num,GFFlists,[dummy,dummy],[0,0,0],2)	
			genesequence=gg[1]
			ORFlength=len(genesequence)-(gg[6]+gg[5])
			UTR5len=gg[5]
			UTR3len=gg[6]
					
			regioncheck=[0,0]
			regioncheck[0]=regionwindow[0]
			regioncheck[1]=regionwindow[1]
			
			# Check hit position
			mrnapos=int(hits[hit][5])
			
			if region==0:
				regionlen=UTR5len
			elif region==1:
				regionlen=ORFlength
			elif region==2:
				regionlen=UTR3len	
			# Convert regioncheck if it's from 3' end to 5'end referenced.	
			if regionwindow[0]<0:
				regioncheck[0]=regionlen+regionwindow[0]
			if regionwindow[1]<0:
				regioncheck[1]=regionlen+regionwindow[1]
					
			if mrnapos>=regioncheck[0] and mrnapos<=regioncheck[1]:
				goodhits.append(hits[hit])
		
		
		hits=goodhits
		print("Actual hits after targets out of regionwindow removed = "+str(len(hits)-1))

	gentools.writelisttoexcel(hits,outfilestring)

# Helper function for findseqinstance2.
def givestopcodon(chrom,feat_num,GFFlists,dummy):
	if dummy=={}:	# In case need to make dummy on fly.
		for sequence in GFFlist:
			dummy[sequence]=[0 for x in range(len(GFFlist[sequence]))]
	
	gg=seqtools.givegene(chrom,feat_num,GFFlists,[dummy,dummy],["none","none",0],0)	
	codon=str(gg[1][-3:])

	if codon!="TAA" and codon!="TAG" and codon!="TGA":
		return -1
	
	return str(codon)
	
# Outputs a yeast transcriptome fasta with the first number being the first base of ORF and the 2nd number being the first base of the 3'UTR.
# This was created in 2019 as a way to know what regions 40S data map to (quantregions script). 
# This was updated Aug 2020 for new givegene. 
def makeyeasttranscriptome(GFFgen_filename,utrgfffilename5,utrgfffilename3,outputfile):
	GFFgen=GFF.parse(GFFgen_filename)
	GFFlist=seqtools.makeGFFlist(GFFgen)
	utrgffgen3=GFF.parse(utrgfffilename3)
	utrtable3=seqtools.makeutrtable(utrgffgen3)
	utrgffgen5=GFF.parse(utrgfffilename5)
	utrtable5=seqtools.makeutrtable(utrgffgen5)
	f=open(outputfile,"w")
	
	shift=0

	#dummy counts
	outputdata1={}
	outputdata2={}
	for key in list(GFFlist.keys()):
		outputdata1[key]=[0 for x in range(len(GFFlist[key]))]
		outputdata2[key]=[0 for x in range(len(GFFlist[key]))]
	counts=[outputdata1,outputdata2]
	
	writtenfeat_num=0
	for chrom in GFFlist:
		feat_num=0
		
		for feature in GFFlist[chrom].features:

			gg=seqtools.givegene(chrom,feat_num,[GFFlist,utrtable5,utrtable3],counts,[0,0,shift],1)
			if gg[0]==-1 or gg[0]==-2 or gg[1]==-1 or gg[1]==-2:
				feat_num+=1
				continue
				
			else:
				genesequence=str(gg[1])
				bp5=gg[5]
				bp3=gg[6]
				if bp5==0:
					noutr5=True
				else:
					noutr5=False
				if bp3==0:
					noutr3=True
				else:
					noutr3=False
			
			descrip=feature.qualifiers["Name"][0]
			seqid=feature.id
			
			end=len(genesequence)-bp3
			start=bp5 
			
			if noutr3==False and noutr5==False:
				seqrec=SeqRecord.SeqRecord(Seq(genesequence))    
				seqrec.id=seqid+"_"+str(start)+"_"+str(end)
				seqrec.description=descrip
				f.write(seqrec.format("fasta"))
			
				writtenfeat_num+=1
			feat_num+=1

	print("Transcripts written ="+str(writtenfeat_num))			
	f.close()




##### Next 3 functions are for mapability analysis.

# Function that makes fastq file of reads of certain length (usually 28nt) for all possible positions in an input sequence.
# Note that no shifting is done so remember to do that.
def makegenefastq(inputsequence,seqname,outputfastqfilebase,readsize):
	f_ORF=open(outputfastqfilebase+"_ORF.fastq","w")
	orfs=[]
	for ntpos in range(len(inputsequence)-readsize+1):
		seqrec=SeqRecord.SeqRecord(Seq.Seq(inputsequence[ntpos:ntpos+readsize]))   
		seqrec.id=seqname
		descrip=seqname
		seqrec.description=descrip+"_"+str(ntpos)
		seqrec.letter_annotations["phred_quality"]=[30 for x in range(readsize)]
		orfs.append(seqrec)
           	
	SeqIO.write(orfs,f_ORF,"fastq")

# Caller of function above with read in fasta file.
def makefastafastq(inputfasta,seqname,outputfastqfilebase,readsize):
	f=SeqIO.parse(inputfasta,"fasta")
	fseq=str(f.next().seq)
	makegenefastq(fseq,seqname,outputfastqfilebase,readsize)



#Function that makes a fastq file of reads of a certain length for all possible positions in genome.
# readsize is size of fake reads, usually 28 or 16 nt.
# Updated Aug 2020 for new givegene.
def makegenomefastq(shift,GFFlist,utrtable5,utrtable3,outputfastqfilebase,readsize):
	f_ORF=open(outputfastqfilebase+"_ORF.fastq","w")
	f_UTR=open(outputfastqfilebase+"_UTR3.fastq","w")


	
	#dummy counts
	outputdata1={}
	outputdata2={}
	for key in list(GFFlist.keys()):
		outputdata1[key]=[0 for x in range(len(GFFlist[key]))]
		outputdata2[key]=[0 for x in range(len(GFFlist[key]))]
	counts=[outputdata1,outputdata2]
	
	for chrom in GFFlist:
		print(chrom)
		feat_num=0
		for feature in GFFlist[chrom].features:
			utrs=[]
			orfs=[]
		
			# FOR NOW, we are disallowing overlaps, since we rejecting those anyway!
			gg=seqtools.givegene(chrom,feat_num,[GFFlist,utrtable5,utrtable3],counts,[0,0,shift],2)
			if gg[0]==-1 or gg[0]==-2 or gg[1]==-1 or gg[1]==-2:
				feat_num+=1
				continue
				
			else:
				genesequence=str(gg[1])
				genecounts=gg[0]
				bp3=gg[6]
				if bp3==0:
					noutr=True
				else:
					noutr=False
			
			descrip=feature.qualifiers["Name"][0]
			seqid=feature.id
			
			end=len(genesequence)-bp3
			
			st=0
			
			for ntpos in range(st,end-readsize+1):
				seqrec=SeqRecord.SeqRecord(Seq.Seq(genesequence[ntpos:ntpos+readsize]))   
				seqrec.id=seqid
				seqrec.description=descrip+"_"+str(ntpos)
				seqrec.letter_annotations["phred_quality"]=[30 for x in range(readsize)]
				orfs.append(seqrec)
                        
			st=end
			if noutr==False:
				for ntpos in range(st,len(genesequence)-readsize+1):
					seqrecutr=SeqRecord.SeqRecord(Seq.Seq(genesequence[ntpos:ntpos+readsize]))    
					seqrecutr.id=seqid
					seqrecutr.description=descrip+"_"+str(ntpos)
					seqrecutr.letter_annotations["phred_quality"]=[30 for x in range(readsize)]
					utrs.append(seqrecutr)
			
			feat_num+=1
	
			SeqIO.write(orfs,f_ORF,"fastq")
			SeqIO.write(utrs,f_UTR,"fastq")

	f_ORF.close()
	f_UTR.close()
















