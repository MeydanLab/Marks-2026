from BCBio import GFF
import gentools
import csv
import seqtools
import builddense

# There are 3 major functions here: posavg, posstats, and motifavg.
# They make averages of density at known positions (posavg) or compute individual pause scores for those positions (posstats).
# The function motifavg computes averages and pause scores from those averages for every motif of a particular length.

# Simple workflow for makeposavg. It averages ribosome profiling data around sites of interest.
# genelist - sites to average, a csv with gene names in col 0, chrom id in 2, feature number in 3, and chromosome position in 4. *** NOTE ***> In 5 is the mRNA position that can be used by commenting that line out below. For 5'UTR it is required.
# GFFfile - path and name of GFF file
# utrgfffile5/3 - this is either the GFF for the UTRs.
# seqwin - this is a list of the number of nt upstream and downstream of your feature of interest to include in the average. If seqwin is negative, it allows "spillover" into other areas. For example, a motif just passed the stop codon would no longer be ignored, but it would be contaminated with ORF reads on the upstream side.
# densityfile - This is the path and fileroot to the ribosome profiling data.
# outfilestring - root of the csv file you will have generated
# riboshift - Amount to shift density before analysis.
# thresh is the number of rpm counts required in the seqwin window of interest.
# ORFthresh is minimal threshold of the main ORF required, in rpkm units.
# normmethod has 3 options: 0 is for normalizing to the window (so every motif is equally weighted). 1 is for normalizing to the main ORF counts (so weighted by ORF expression), and 2 is to unequally weight (no normalization) the motifs.
# Region is 0, 1, or 2 (5'UTR, ORF, 3'UTR). Legacy option for -2.
# Input sites for makeposavg should be in coordinates that start with the region of interest (region variable): 5'UTR (0), ORF (1), or 3'UTR (2). These coordinates are defined by running givegene using the Cherry Lab's published UTRs annotation based on the Steinmetz lab data.

def makeposavg_wf0(genelist,GFFfile,utrgfffile5,utrgfffile3,seqwin,densityfiles,outfilestring,riboshift,thresh,ORFthresh,normmethod,region,bp5,bp3):
	gffgen=GFF.parse(GFFfile)
	GFFlist=seqtools.makeGFFlist(gffgen)
	
	utrgff5=GFF.parse(utrgfffile5)
	utrtable5=seqtools.makeutrtable(utrgff5)
	utrgff3=GFF.parse(utrgfffile3)
	utrtable3=seqtools.makeutrtable(utrgff3)
	GFFlists=[GFFlist,utrtable5,utrtable3]
	
	writerfile=open(outfilestring+".csv", "wb")
	writer = csv.writer(writerfile)
		
	for densityfile in densityfiles:
		counts1p=builddense.readcountsf(densityfile+"_plus_")
		counts1m=builddense.readcountsf(densityfile+"_minus_")
		readcounts=[counts1p,counts1m]
		seqwin_in=list(seqwin)		# Have to make a separate list since it changes in main func.
		genesinavg=makeposavg(genelist,GFFlists,seqwin_in,readcounts,riboshift,thresh,ORFthresh,normmethod,region,bp5,bp3)
		samplename=densityfile.split("/")[-1]
		writer.writerow([samplename]+genesinavg[6])
		print(densityfile)
		print("positions in avg = "+str(genesinavg[0]))
		print("positions not in avg because zero count in window = "+str(genesinavg[1]))
		print("positions not in avg due to zero ORF count = "+str(genesinavg[2]))
		print("positions exempted because too close to end of gene = "+str(genesinavg[3]))
		print("positions exempted because bad gene in some way = "+str(genesinavg[4]))
		print("positions not in avg because normalization value was 0 = "+str(genesinavg[5]))
	
	writerfile.close()	
	gentools.transposecsv(outfilestring)

# This function will take a list with gene names and positions and a window over which to average and then output the average as a binary.
# input genelist is a csv gene names in col 1, first column alias, second chrom id, third featurenum, and 5th position. Fourth is ignored at this point.
def makeposavg(genelist,GFFlists,seqwin,readcounts,riboshift,thresh,ORFthresh,normmethod,region,bp5,bp3):
	
	zerocountnormcount=0
	zeroORFcount=0
	count=0
	zerocount=0
	tooclosecount=0
	badgene=0
	
	seqcheck=[0,0]
	if seqwin[0]<0:
		seqwin[0]*=-1
		seqcheck[0]=1
	 	
	if seqwin[1]<0:
		seqwin[1]*=-1
		seqcheck[1]=1
	
	f_csv=open(genelist)
	pausedict=gentools.readindict(f_csv)
	avgcounts=[0 for x in range(seqwin[0]+seqwin[1])]
	
	for genename in pausedict:
		if genename=="headers":
			continue
		chrom=pausedict[genename][1]
		feat_num=int(pausedict[genename][2])	
		mrnaposition=int(pausedict[genename][4]) 
		# Note the genome position could also be used if converted by convertmRNAtogenomic. That's not implemented because that increases chance of errors.	
	
		######### On 10/15/16, goodgenes was changed to 1 because the ORF normalization made it clear that goodgenes=2 includes a handful of overlapping genes that dominate the averages.				
		gg=seqtools.givegene(chrom,feat_num,GFFlists,readcounts,[bp5,bp3,riboshift],1)			
		if gg[0]==-1 or gg[0]==-2:
			badgene+=1
			continue
		
		UTR5len=gg[5]
		UTR3len=gg[6]
		genecounts=gg[0]
		genelength=len(gg[0])	

		if region==0:		
			tooclosetostart=0
			tooclosetostop=UTR5len
		elif region==1:		# ORF
			mrnaposition+=(UTR5len)
			tooclosetostop=(genelength-UTR3len)	
			tooclosetostart=UTR5len
		elif region==2:
			mrnaposition+=(genelength-UTR3len)
			tooclosetostop=genelength
			tooclosetostart=(genelength-UTR3len)
		elif region==-2:
			mrnaposition+=(UTR5len)
			tooclosetostop=genelength
			tooclosetostart=(genelength-UTR3len)
		else:
			exit()
			
		windowscheck=0		# Checking to make sure no average outside region of interest go into the average (unless the seqcheck got set to 1).
		if ((seqcheck[0]==1 or (mrnaposition-tooclosetostart)>=seqwin[0]) and (seqcheck[1]==1 or (tooclosetostop-mrnaposition)>=seqwin[1])):
			windowscheck=1
			
		if windowscheck==1:
			if (mrnaposition-seqwin[0] < 0) or (mrnaposition+seqwin[1] > len(genecounts)):
				tooclosecount+=1
				continue
			else:
				loccounts=genecounts[mrnaposition-seqwin[0]:mrnaposition+seqwin[1]]
			ORFcounts=genecounts[UTR5len:genelength-UTR3len]
			
			if sum(loccounts)>=thresh and (1000*sum(ORFcounts)/len(ORFcounts))>=ORFthresh:	# ORFthresh is in rpkm
				count+=1	
			elif sum(loccounts)>=thresh:
				zeroORFcount+=1	# If both thresholds are missed, gets recorded as zeroORFcount since more important to know that. And matters for normalization below.
				continue
			elif (1000*sum(ORFcounts)/len(ORFcounts))>=ORFthresh:
				zerocount+=1
				continue
			else:
				 zeroORFcount+=1
				 continue                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
		else:
			tooclosecount+=1
			continue

		# Make average
		if normmethod==2:
			normval=1	
		elif normmethod==1:
			normval=sum(ORFcounts)/len(ORFcounts)	# This is by nt, not kb.
		elif normmethod==0:
			normval=sum(loccounts)			
		else:
			exit()
			
		if normval==0:
			zerocountnormcount+=1
			continue	
			
		for i in range(len(avgcounts)):
			avgcounts[i]+=loccounts[i]/float(normval)

	if count==0:
		print("No genes to average.")
		return [count,zerocount,zeroORFcount,tooclosecount,zerocountnormcount,badgene,avgcounts]
			
	for i in range(len(avgcounts)):
		if normmethod==0:
			avgcounts[i]/=count					# Doesn't make sense to count loccounts 0s in something that is self-normalized, equally weighted.
		elif normmethod==1:
			avgcounts[i]/=count		# Now we can't include 0 ORFcounts because they are likely artifacts (would weight to infinity). 
		elif normmethod==2:
			avgcounts[i]/=(count+zerocountnormcount)	# Now all are included since unequally weighted, so all should go in.

	return [count,zerocount,zeroORFcount,tooclosecount,zerocountnormcount,badgene,avgcounts]
	
	
# Simple workflow for makeposstats
def makeposstats_wf(filenames,genelist,GFFfile,utrgfffile5,utrgfffile3,seqwin,pausewin,outfilestring,riboshift,siteshift,region):
	
	gffgen=GFF.parse(GFFfile)
	GFFlist=seqtools.makeGFFlist(gffgen)
	utrgff5=GFF.parse(utrgfffile5)
	utrtable5=seqtools.makeutrtable(utrgff5)
	utrgff3=GFF.parse(utrgfffile3)
	utrtable3=seqtools.makeutrtable(utrgff3)
	GFFlists=[GFFlist,utrtable5,utrtable3]
	
	outputs=[]
	firsttime=0
	for fname in filenames:
		counts1p=builddense.readcountsf(fname+"_plus_")
		counts1m=builddense.readcountsf(fname+"_minus_")
		readcounts=[counts1p,counts1m]
		samplename=fname.split("/")[-1]
		localoutput=makeposstats(genelist, GFFlists, seqwin, pausewin, readcounts, samplename,riboshift, siteshift, region)
		
		if firsttime==0:
			numoutputs=len(localoutput)-1
			genedictionarygeneral=localoutput[0]
			for i in range(numoutputs):
				outputs.append(localoutput[i+1])	
			firsttime=1
		else:
			for i in range(numoutputs):
				fillins=len(localoutput[i+1]["headers"])
				outputs[i]=gentools.mergedictionaries(outputs[i],localoutput[i+1],fillins,-10)
			
	for i in range(numoutputs):
		fillins=len(outputs[i]["headers"])
		outputs[i]=gentools.mergedictionaries(genedictionarygeneral,outputs[i],fillins,-10)	
		gentools.writedicttoexcel(outputs[i],outfilestring+"_output"+str(i+1))	

# This script computes a pause score at the input sites of interest by taking the reads locally around the peak and dividing by the average reads in the gene.	
# GFFfile - path and name of GFF file
# utrgfffile5/3 - this is either the GFF for the UTRs.
# seqwin - this is a list of the number of nt upstream and downstream of your feature of interest to include in the reported sequence in the output. 
# pausewin is a list of 2 numbers for nt upstream/downstream of the peak to include in the numerator of the pause score.
# readcounts - This is the ribosome profiling density.
# outfilestring - root of the csv file you will have generated
# riboshift - Amount to shift density before analysis.
# put 0 for siteshift if you want to keep it at the 5' end of the motif. Put in motif length if you want 3' end. And subtract 3 if you want 3 in (ie A site since shift shifts for P site).
# Region is 0, 1, or 2 (5'UTR, ORF, 3'UTR). Legacy option for -2.
# Input sites should be in coordinates that start with the region of interest (region variable): 5'UTR (0), ORF (1), or 3'UTR (2). These coordinates are defined by running givegene using the Cherry Lab's published UTRs annotation based on the Steinmetz lab data.
# genelist - sites to average, a csv with gene names in col 0, chrom id in 2, feature number in 3, and chromosome position in 4. In 5 is the mRNA position that is also required.
# samplename - is name to identify the input file.
def makeposstats(genelist, GFFlists, seqwin, pausewin, readcounts, samplename,riboshift, siteshift, region):
	badgene=0
	tooclosecount=0
	count=0
	outdict={}
	outdict["headers"]=[]
	
	GFFlist=GFFlists[0]
	f_csv=open(genelist)
	
	pausedict=gentools.readindict(f_csv)
	
	pausedict["headers"].append("LocalSeq_"+samplename)
	pausedict["headers"].append("LocalSeq_trans_"+samplename)
	
	outdict["headers"].append("PauseScore_"+samplename)
	outdict["headers"].append("num_rpkm_"+samplename)
	outdict["headers"].append("denom_rpkm_"+samplename)

	for genename in pausedict:
		if genename=="headers":
			continue
		
		outdict[genename]=[]
		chrom=pausedict[genename][1]
		feat_num=int(pausedict[genename][2])
		mrnaposition=int(pausedict[genename][4])
		
		mrnaposition=mrnaposition+siteshift	# This is here so anywhere in motif can be put in the proper site of ribosome and shift corrected. 
		gg=seqtools.givegene(chrom,feat_num,GFFlists,readcounts,[0,0,riboshift],2)	# Allowing overlap. May wish to change this since it mattered for posavg.		
		if gg[0]==-1 or gg[0]==-2:
			badgene+=1
			continue
		UTR5len=gg[5]
		UTR3len=gg[6]
		genecounts=gg[0]
		genelength=len(gg[0])	
		genesequence=gg[1]

		if region==0:		
			tooclosetostart=0
			tooclosetostop=UTR5len
		elif region==1:		# ORF
			mrnaposition+=(UTR5len)
			tooclosetostop=(genelength-UTR3len)	
			tooclosetostart=UTR5len
		elif region==2:
			mrnaposition+=(genelength-UTR3len)
			tooclosetostop=genelength
			tooclosetostart=(genelength-UTR3len)
		elif region==-2:
			mrnaposition+=(UTR5len)
			tooclosetostop=genelength
			tooclosetostart=(genelength-UTR3len)
		else:
			exit()
			
		# Checking to make sure no average outside region of interest go into the average.
		if ((mrnaposition-tooclosetostart)>=pausewin[0] and (tooclosetostop-mrnaposition)>=pausewin[1]):
			count+=1
		else:
			tooclosecount+=1
			continue
		
		# Get sequence
		if ((mrnaposition-tooclosetostart)>=seqwin[0] and (tooclosetostop-mrnaposition)>=seqwin[1]):
			locseq=genesequence[mrnaposition-seqwin[0]:mrnaposition+seqwin[1]]
			adjustval0=(seqwin[0])%3
			adjustval1=(seqwin[1])%3
			transeq=genesequence[mrnaposition+adjustval0-seqwin[0]:mrnaposition+seqwin[1]-adjustval1].translate()
		else:
			locseq="Outsideofgene"
			transeq="Outsideofgene"

		# Compute pause score by calling module. Can compute other things (ie frame info, or use different denom) via modules.
		pauseinfo=pausemodule(genecounts,pausewin,region,UTR5len,UTR3len,mrnaposition)
		pausedict[genename].append(locseq)
		pausedict[genename].append(transeq)		
		outdict[genename].append(pauseinfo[2])
		outdict[genename].append(pauseinfo[0])	# in rpkm
		outdict[genename].append(pauseinfo[1])

	return [pausedict,outdict]		#future could add other outputs
	

def pausemodule(genecounts,pausewin,region,UTR5len,UTR3len,mrnaposition):
	if region==0:
		pause_denom=sum(genecounts[0:UTR5len])/UTR5len
	elif region==1:
		pause_denom=sum(genecounts[UTR5len:len(genecounts)-UTR3len])/(len(genecounts)-UTR3len-UTR5len)
	elif region==2 or region==-2:
		pause_denom=sum(genecounts[len(genecounts)-UTR3len:])/UTR3len
	pause_numerator=sum(genecounts[mrnaposition-pausewin[0]:mrnaposition+pausewin[1]+1]) # Note includes middle base (+1).
	pause_numerator/=(pausewin[1]+pausewin[0]+1)			
	if pause_denom>0:
		pausescore=pause_numerator/pause_denom
	else:
		pausescore=-1
				
	# To convert to rpkm
	pause_numerator*=1000
	pause_denom*=1000			 
	
	return([pause_numerator,pause_denom,pausescore])
	
	
	
# CODING NOTES
# motifavg2 is a highly experimental function that is under continued development.
# In the current code, it computes the average plot for every possible motif of a given length (either aa or nt).
# This is output as a binary file.
# It also computes the pause scores from those plots. 
# Originally, this code did more (pause scores and drop off scores) at each position, but this wasn't very robust. (See versions in seqtools_extra for this).
# Code now is same as that developed for 1st Github release, but made to work with new givegene in Aug 2020.
# In the future, the code here could be merged with (or work on outputs of) posavg, and maybe posstats.
# The current difference between this and posavg/posstats is that those functions work on separately-created (and therefore more customized) lists of input motifs.
# Also, the code here is also updated from Github to handle any UTR region and hardshift as an input now.

# USAGE NOTES
# What it reports is essentially an enrichment score for every position type, ie every amino acid. It is very useful for determining if any particular amino acid or codon is stalling the ribosome in your sample.
# This function scans the entire transcriptome for all possible sequence motifs of a particular length (nt or amino acid).
# It returns an average plot binary for each and then computes the pause score for each.
# GFFfile - path and name of GFF file
# utr5gfffile - Not used at this time. Input ignored.
# utr3gfffile - Not used at this time. Input ignored.
# counts_filestring - This is the path and fileroot to the ribosome profiling data.
# motifsize - size of motif in nt or amino acids, depending on inframe setting
# inframe - Set to 0 if searching by nt in all frames, set to 1 if searching only main frame in nt, set to 2 if searching main frame by amino acid.
# thresh - This variable is not used at this time. 
# outfilestring - rootfile and path for output of script.
# mismatches - This variable is not used at this time. 
# shift - Amount to shift your ribosome profiling data to align it with the peak of interest.
# windowsize - Number of nt to add to either side of pause peak for computation of pause score. For example, 2 would give a peak of 2+2+1 = 5 nt total.
# avgwindow -  A list of 2, nt distance upstream and downstream of peak of interest to include in the average plot and to use for the denominator of the pause score. Note that motifs positioned at the ORF ends that do not "fit" these windows are excluded.
# hardshift - a value that can be used to filter which regions of gene are used. Set to shift for no adjustment. Set to 13 to avoid stop codon regions for 3' end alignment.
# typeregion - whether ORF, UTR5, or UTR3. Generally we use ORF, but others are available.
def motifavg_2_wf(GFFfile,utr5gfffile,utr3gfffile,counts_filestring,motifsize,inframe,thresh,outfilestring,mismatches,shift,windowsize,avgwindow,hardshift,typeregion):
	codons={} # Just a placeholder
	counts0=builddense.readcountsf(counts_filestring+"_plus_")		
	counts1=builddense.readcountsf(counts_filestring+"_minus_")
	counts=[counts0,counts1]
	GFFgen=GFF.parse(GFFfile)
	GFFlist=seqtools.makeGFFlist(GFFgen)
	utrgffgen=GFF.parse(utr5gfffile)
	utrtable5=seqtools.makeutrtable(utrgffgen)
	utrgffgen=GFF.parse(utr3gfffile)
	utrtable3=seqtools.makeutrtable(utrgffgen)
	GFFs=[GFFlist,utrtable5,utrtable3]
	
	avglist=[]
	outlist=[]
	outlist.append(["motif","na","na","na","na","hitsincluded","na","na","tothits","Pausescore"])
	
	motifdata=motifavg_2_simple(GFFs,motifsize,inframe,thresh,mismatches,codons,shift,counts,windowsize,avgwindow,hardshift,typeregion)
		
	for mm in list(motifdata.keys()):
		print(mm)
		motifdata[mm][1]=0 # Variable not used
		if motifdata[mm][5]>0:
			for i in range(sum(avgwindow)):
				motifdata[mm][0][i]/=float(motifdata[mm][5])	#Normalization taking place.
		motifdata[mm][3]=0	# Variable not used
		
		# Get pause scores:
		numerator=sum(motifdata[mm][0][avgwindow[0]-windowsize:avgwindow[0]+1+windowsize])
		denominator=sum(motifdata[mm][0])
		denominator/=len(motifdata[mm][0])
		numerator/=(2*windowsize+1)
		if denominator!=0:
			pause=numerator/denominator
		else:
			pause=0
		motifdata[mm][10]=pause

		outlist.append([mm]+motifdata[mm][1:9]+[motifdata[mm][10]])
		avglist+=(motifdata[mm][0])	#Concatenate average files.
	
	# WRite out list.
	gentools.writelisttoexcel(outlist,outfilestring)		#Includes new pause scores.
	# Write out avg file.
	favg=open(outfilestring+".bin","wb")
	gentools.writelistbinint(avglist,favg)
	favg.close()

def motifavg_2_simple(GFFs,motif,inframe,thresh,mismatches,codons,shift,countslist,windowsize,avgwindow,hardshift,typeregion):
	motifdata={}	
	shiftdif=shift-hardshift
	if inframe==2:
		motiflen=(int(motif))*3
	elif inframe<=1:
		motiflen=int(motif)
	else:	
		print("illegal input for inframe.") 
		return
	GFFlist=GFFs[0]
	print("Motif is: "+str(motif)+".")

	if (typeregion=="UTR5" or typeregion=="UTR3") and inframe!=0:
		print("Not looking at all frames in UTR region. Are you sure this is correct?")

# Call givegene for every gene in genome, one time.	
	for chrom in GFFlist:
		feat_num=0
		print(chrom)
		for feature in GFFlist[chrom].features:
			
			gg=seqtools.givegene(chrom,feat_num,GFFs,countslist,[0,0,shift],2)	
			gg_ref=seqtools.givegene(chrom,feat_num,GFFs,countslist,[0,0,hardshift],2)		# We do a second call to give gene to get the background reads for a gene. At some level, doing this is probably overly careful and might be eliminated in a future version.	
			if gg[1]==-1 or gg[1]==-2 or gg_ref[1]==-1 or gg_ref[1]==-2:		# For genes that are considered not allowed because they are dubious, nongene, wrong chrom, etc. 
				feat_num+=1
				continue
			bp5=gg[5]
			bp3=gg[6]
			
			if typeregion=="UTR5" and bp5!=0:
				gg[0]=gg[0][0:bp5]
				gg[1]=gg[1][0:bp5]
				gg_ref[0]=gg_ref[0][0:bp5]
				gg_ref[1]=gg_ref[1][0:bp5]
				
			elif typeregion=="UTR3" and bp3!=0:
				gg[0]=gg[0][-bp3:]
				gg[1]=gg[1][-bp3:]
				gg_ref[0]=gg_ref[0][-bp3:]
				gg_ref[1]=gg_ref[1][-bp3:]
			
			elif typeregion=="ORF":	#ORF
				gg[0]=gg[0][bp5:len(gg[0])-bp3]
				gg[1]=gg[1][bp5:len(gg[1])-bp3]
				gg_ref[0]=gg_ref[0][bp5:len(gg[0])-bp3]
				gg_ref[1]=gg_ref[1][bp5:len(gg[1])-bp3]
			
			else:	# no UTR info:
				feat_num+=1
				continue
	
			genelen=len(gg[0])
			if genelen==0:
				print("Zero length gene, ERROR!")
			seqlen=genelen
			
			position=0
			while position<(seqlen-motiflen):	
				currentsequence=gg[1][position:position+motiflen]
				if inframe==2:
					currentsequence=str(currentsequence.translate())
				else:
					currentsequence=str(currentsequence)
					
				# Check if motif made for this yet. And increment total count
				if currentsequence not in motifdata:
					motifdata[currentsequence]=[[0 for x in range(sum(avgwindow))],[],0,[],0,0,0,0,0,[],0,[],[]]
					motifdata[currentsequence][8]+=1
				else:
					motifdata[currentsequence][8]+=1
				
				if(sum(avgwindow)>0):
					# Add in counts for binary output of averaged hits.
						if (-avgwindow[0]+position)>=shiftdif and (-avgwindow[0]+position)>=0 and (avgwindow[1]+position)<=(genelen+shiftdif) and (avgwindow[1]+position)<=genelen:	
							normfactor=sum(gg[0][-avgwindow[0]+position:avgwindow[1]+position])	# So every component of average is weighted equally.
							if normfactor>0:
								for i in range(sum(avgwindow)):
									motifdata[currentsequence][0][i]+=((gg[0][i-avgwindow[0]+position])/normfactor)	#NOTE: output will have the first position of motif (shifted by shift) at the first point after midpoint, ie position 35 if avgwindow ranges 0-69.
								motifdata[currentsequence][5]+=1	# Increment count of averaged hits.
				if inframe>=1:
					position+=3
				if inframe==0:
					position+=1
			feat_num+=1
	return motifdata
