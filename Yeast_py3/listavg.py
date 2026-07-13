# This file contains code for creating average (metagene) plots and quantifying reads on genes (counting, pauses).
# There is a main piece of code, wrapper (workflow) functions, and helper functions.

from BCBio import GFF
import struct
import builddense
import seqtools
import gentools

# Workflow function to handle all the file operations. It calls totalquant to add up counts and compute pause scores in genes.
def totalquant_wf(filenames,output,GFFgen_filename,utrgfffilename,utr5gfffilename,bp5_0,bp3_0,ignoreutr5,ignoreutr3,shift,normalizepause,eliminate,start_offset,end_offset,startpauseoffset,termpauseoffset,termstackpauseoffset,pausehalfregion,covlength,future1,future2):
	outputs=[]		# This will be a list of merged dictionaries that will be written to file in the end.
	
	# This function is basically calling totalquant, merging the output dictionaries together, and making sure the general dictionary (which is provided every time) is added on once at the end
	firsttime=0
	for fname in filenames:		
		localoutput=totalquant(fname,GFFgen_filename,utrgfffilename,utr5gfffilename,bp5_0,bp3_0,ignoreutr5,ignoreutr3,shift,normalizepause,eliminate,start_offset,end_offset,startpauseoffset,termpauseoffset,termstackpauseoffset,pausehalfregion,covlength,future1,future2)
		if firsttime==0:		# First time through get the general dictionary and put the first sample's data values on outputs.
			numoutputs=len(localoutput)-1
			genedictionarygeneral=localoutput[0]	# Create the general dictionary with basic information about each gene.
			for i in range(numoutputs):
				outputs.append(localoutput[i+1])	# The first sample's data values are added in.
			firsttime=1
		else:					# Subsequent times through merge in new data values for each sample.
			for i in range(numoutputs):
				fillins=len(localoutput[i+1]["headers"])	# This counts out how many values there for each gene with a key so you can it can then put in -10s for any gene that has no key in the incoming dictionary.
				outputs[i]=gentools.mergedictionaries(outputs[i],localoutput[i+1],fillins,-10)
			
	for i in range(numoutputs):		# Write out each merged dictionary after first tacking on the generaldictionary information to each. 
		fillins=len(outputs[i]["headers"])
		outputs[i]=gentools.mergedictionaries(genedictionarygeneral,outputs[i],fillins,-10)	
		gentools.writedicttoexcel(outputs[i],output+"_output"+str(i+1))	


# Function that loops through the annotation GFF and counts reads on genes and pause scores at start/stop codons. 
# It can take other modules for added functions. See listavg_extra.py for examples of extra features.
# Variables: 
# counts_filestring is the input density file
# GFFgen_filename is the GFF file name
# utrgfffilename,utr5gfffilename are filenames for the UTR annotations.
# bp5_0,bp3_0 are how much to extend annotated UTRs, or the total length of the UTRs if in manual mode for UTRs.
# ignoreutr5,ignoreutr3 when set to 1 set UTR annotation to manual mode. 0 for annotated.
# shift for ribosome P sites, typically 12 or 13 for 5' aligned reads 
# normalizepause is set to 1 if you want to normalize the pause score to ORF (80S) or 0 to not do that (40S). 
# eliminate is 0 for all genes, 1 for genes not rejected by givegene (including overlaps), 2 for genes not rejected by givegene but includes overlaps.
# start_offset,end_offset - these are the number of nt to ignore at start/end of genes when computing main ORF counts. It eliminates start/stop codon artifact peaks.
# startpauseoffset,termpauseoffset,termstackpauseoffset - corrections to get the shifted P site over the pause. Typically, for 5' end alignments, 0 for start, but 5 for termination, and more for stacked (usually above 30).
# pausehalfregion is an arbitrary value that sets the window for the pause calculation.
# covlength is the length of the read if coverage was used. Otherwise 1. Required for getting raw counts.
# This function was completely revamped on August 8, 2020.

def totalquant(counts_filestring,GFFgen_filename,utrgfffilename,utr5gfffilename,bp5_0,bp3_0,ignoreutr5,ignoreutr3,shift,normalizepause,eliminate,start_offset,end_offset,startpauseoffset,termpauseoffset,termstackpauseoffset,pausehalfregion,covlength,future1,future2):
	GFFgen=GFF.parse(GFFgen_filename)
	GFFlist=seqtools.makeGFFlist(GFFgen)
	counts0=builddense.readcountsf(counts_filestring+"_plus_")		
	counts1=builddense.readcountsf(counts_filestring+"_minus_")
	counts=[counts0,counts1]
	utrgffgen=GFF.parse(utrgfffilename)
	utrtable=seqtools.makeutrtable(utrgffgen)
	utrgffgen=GFF.parse(utr5gfffilename)
	utrtable2=seqtools.makeutrtable(utrgffgen)
	
	# The dictionaries we will output.
	genedictionarygeneral={}	# This will get tacked on to the other dictionaries, and they in turn will be written as csv files.
	genedictionary1={}		
	genedictionary2={}
	#genedictionary3={}
	
	# Get header information in place.
	samplename=counts_filestring.split("/")[-1]
	genedictionarygeneral["headers"]=["alias","chrom","feat_num","note","UTR5seq","mainORFseq","UTR3seq"]
	genedictionary1["headers"]=["CDS_"+samplename,"CDS_RAW_"+samplename,"UTR5_"+samplename,"UTR5_RAW_"+samplename,"UTR3_"+samplename,"UTR3_RAW_"+samplename]
	genedictionary2["headers"]=["startpause_"+samplename,"termpause_"+samplename,"termstackpause_"+samplename]
	#genedictionary3["headers"]=["param1"]
	
	# Set all the counters to 0.
	illegalgenes=0
	overlapgenes=0
	genesinlist=0
	noutrgenes=0
	otherbadgenes=0
	
	# Safety check, make sure these aren't invalid.
	if bp5_0<0 or bp3_0<0:
		print("error in bp lengths")
		exit()

	# Get read depth by finding minimal counts at a position. And start preparing the general dictionary.
	mincountsglobal=1000
	for chrom in GFFlist:
		feat_num=0
		for feature in GFFlist[chrom].features:
			#print(feature)
			# Get Alias and note.
			if "Alias" in feature.qualifiers:
				alias = feature.qualifiers["Alias"][0]
			else:
				alias = "NA"	
			if "Note" in feature.qualifiers:
				note = feature.qualifiers["Note"][0]
			else:
				note = "NA"	
			genedictionarygeneral[feature.id]=[alias,chrom,feat_num,note," "," "," "]

			# Call givegene to get the sequence and counts for the spliced ORF.
			gg=seqtools.givegene(chrom,feat_num,[GFFlist,utrtable2,utrtable],counts,[0,0,shift],1)
			#print(gg)
			#print(gg[0])
			try:
				if gg[0]<0:
					print("success")
			except TypeError:
				print(gg)
				break
			if gg[0]<0:
				feat_num+=1
				continue
			zeroremoved=[i for i in gg[0] if i!=0]
			if len(zeroremoved)==0:
				feat_num+=1
				continue
			mincounts=min(zeroremoved)
			if mincounts<mincountsglobal:
				mincountsglobal=mincounts
			feat_num+=1

	readsmapped=int(round((1/(mincountsglobal*covlength))*1E6))
	rawconvert=1E6/float(readsmapped)		# Factor to be used for converting rpm to raw reads.
	print("Computed reads mapped = "+str(readsmapped)+" for "+counts_filestring)

	# For manual UTRs, the bp_0 values have to be negative so givegene knows these are manual UTRs, not extensions off annotated UTRs.
	# And for 0-length manual UTRs, then "none" has to be used for givegene.
	if ignoreutr3==1:
		bp3_0*=-1
		if bp3_0==0:
			bp3_0="none"
	if ignoreutr5==1:
		bp5_0*=-1
		if bp5_0==0:
			bp5_0="none"
			
	# Start looping through the genome.
	for chrom in GFFlist:
		feat_num=0		# Start the feature counter at 0.
		for feature in GFFlist[chrom].features:
					
			# Call givegene to get the sequence and counts for the spliced gene.
			gg=seqtools.givegene(chrom,feat_num,[GFFlist,utrtable2,utrtable],counts,[bp5_0,bp3_0,shift],eliminate)
			genesequence=gg[1]
			genecounts=gg[0]
			bp5=gg[5]
			bp3=gg[6]
			
			#####
			# These intron parameters are not used but can be analyzed with add-on functions. This information corresponds to the first intron in the main ORF.
			# These values will be -1 if the gene has no CDS intron.
			introncounts=gg[2]	# List of counts on intron.
			intronsequence=gg[3]	# Sequence of intron.
			intronframe=gg[4]	# This is the frame the intron begins with (0, 1, or 2).
			#####
			
			# Get rid of dubious genes, nongenes, genes with overlap of others and report appropriately. Note that if eliminate is set to 2, -2 will never be reported here.
			if (genesequence==-1 or genecounts ==-1 or genesequence==-2 or genecounts==-2):
				feat_num+=1
				if genesequence==-1:
					illegalgenes+=1
					continue
				if genesequence==-2:
					overlapgenes+=1
					continue
				else:
					print("Illegal output of givegene.")
					print(gg)
					exit()
			
			# Continue if the UTR annotation wasn't there (and we aren't in manual). This eliminates genes that don't have both UTR5 and UTR3 annotation.
			if bp5_0==bp5 and ignoreutr5!=1:
				feat_num+=1
				noutrgenes+=1
				continue
			if bp3_0==bp3 and ignoreutr3!=1:
				feat_num+=1
				noutrgenes+=1
				continue

			# Define ORF
			start=bp5
			end=len(genecounts)-bp3
			start_special=start+start_offset
			end_special=end-end_offset
			
			# To check for any cases with really short genes that are lost with the special offsets:
			if end_special<start_special:
				print("start_special or end_special too big for gene:")
				print(feature.id)
				feat_num+=1
				otherbadgenes+=1
				continue
			
			# Define UTRs.
			UTR5genesequence=genesequence[0:start]	
			UTR5genecounts=genecounts[0:start]
			UTR3genecounts=genecounts[end:len(genecounts)]
			UTR3genesequence=genesequence[end:len(genecounts)]	
			
			# DO THE COUNT
			genedictionary1[feature.id]=getcounts(genecounts,UTR5genecounts,UTR3genecounts,start_special,end_special,rawconvert)
			
			# DO PAUSE SCORES
#			stopcodon=genesequence[end-3:end]	# In case you want the stop codon, here it is.
			genedictionary2[feature.id]=getpausescores(genecounts,startpauseoffset,termpauseoffset,termstackpauseoffset,start_special,end_special,normalizepause,pausehalfregion,end,start)
			
			# OTHER FUTURE ADD-ONS in listavg_extra or other new add-on python file.
			#import newaddon
			#genedictionary3[feature.id]=newaddon.newfunction(genecounts,)	
			
			# Finish general dictionary with basic information about sequence.
			genedictionarygeneral[feature.id][4]=str(UTR5genesequence)
			genedictionarygeneral[feature.id][5]=str(genesequence[start:end])
			genedictionarygeneral[feature.id][6]=str(UTR3genesequence)

			feat_num+=1	
			genesinlist+=1
	
	# Report summary information.
	print("Genes with no UTR when requested = "+str(noutrgenes))
	print("Genes with overlap = "+str(overlapgenes))
	print("Genes dropped by givegene (non-mRNA, dubious, or known misannotation) = "+str(illegalgenes))
	print("Other bad genes (too short for special start/end) = "+str(otherbadgenes))
	print("Genes in list = "+str(genesinlist))
	
	return [genedictionarygeneral,genedictionary1,genedictionary2]		# You can add genedictionary3 for reporting other kinds of information if used above.
			

def getcounts(genecounts,UTR5genecounts,UTR3genecounts,start_special,end_special,rawconvert):
	totalgenereads=float(1000)*sum(genecounts[start_special:end_special])/(end_special-start_special)
	totalgenereadsRAW=int(round(sum(genecounts[start_special:end_special])/rawconvert))
	if len(UTR3genecounts)!=0:
		totalUTR3reads=float(1000)*sum(UTR3genecounts)/len(UTR3genecounts)
		totalUTR3readsRAW=int(round(sum(UTR3genecounts)/rawconvert))
	else:
		totalUTR3reads=0
		totalUTR3readsRAW=0
	if len(UTR5genecounts)!=0:
		totalUTR5reads=float(1000)*sum(UTR5genecounts)/len(UTR5genecounts)
		totalUTR5readsRAW=int(round(sum(UTR5genecounts)/rawconvert))
	else:
		totalUTR5reads=0
		totalUTR5readsRAW=0
	
	return [totalgenereads,totalgenereadsRAW,totalUTR5reads,totalUTR5readsRAW,totalUTR3reads,totalUTR3readsRAW]


def getpausescores(genecounts,startpauseoffset,termpauseoffset,termstackpauseoffset,start_special,end_special,normalizepause,pausehalfregion,end,start):
	lengc=len(genecounts)
	# Now define termination pause score. 
	if ((end-termpauseoffset)-pausehalfregion<0 or (end-termpauseoffset)+pausehalfregion+1>lengc):		# Check that pause half region isn't putting the window outside the gene - an issue if UTRs not used.
		termpausescore=-1
	else:
		termpausescore=sum(genecounts[(end-termpauseoffset)-pausehalfregion:(end-termpauseoffset)+pausehalfregion+1])/(2*pausehalfregion+1) 
		
	# Define the start codon pause score. 
	if ((start-startpauseoffset)-pausehalfregion<0 or (start-startpauseoffset)+pausehalfregion+1>lengc):
		startpausescore=-1
	else:
		startpausescore=sum(genecounts[(start-startpauseoffset)-pausehalfregion:(start-startpauseoffset)+pausehalfregion+1])/(2*pausehalfregion+1)
	
	# Define term stack pause score. 
	if ((end-termstackpauseoffset)-pausehalfregion<0 or (end-termstackpauseoffset)+pausehalfregion>lengc):
		termstackpausescore=-1
	else:
		termstackpausescore=sum(genecounts[(end-termstackpauseoffset)-pausehalfregion:(end-termstackpauseoffset)+pausehalfregion+1])/(2*pausehalfregion+1) 	
	
		# Normalize the pause scores.
	if normalizepause==0:
		ORFaverage=1		# New can set a denom of 1 if filtermodule is -1.
	else:
		ORFaverage=sum(genecounts[start_special:end_special])/len(genecounts[start_special:end_special])	# Now using average. So this is not actually a median anymore.
	
	if ORFaverage!=0:		# Normalize if ORF has reads and pause score is valid.
		if termpausescore>=0:
			termpausescore/=ORFaverage
		if termstackpausescore>=0:
			termstackpausescore/=ORFaverage
		if startpausescore>=0:
			startpausescore/=ORFaverage
	else:
		termpausescore=-1
		termstackpausescore=-1
		startpausescore=-1
	
	return [startpausescore,termpausescore,termstackpausescore]


# This is the workflow that creates the average (or "metagene") plots.
def totalavg_wf(filenames,output,GFFgen_filename,utr5gfffilename,utr3gfffilename,regionlength5,regionlength3,ignoreutr,shift,equalweight,eliminate,thresh,avgregion,goodzone):
	outputs={}
	for fname in filenames:
		samplename=fname.split("/")[-1]
		localoutput=makeavggene(fname,GFFgen_filename,utr5gfffilename,utr3gfffilename,regionlength5,regionlength3,ignoreutr,thresh,shift,eliminate,goodzone,equalweight,avgregion)
		outputs[samplename]=list(localoutput)
	gentools.writedicttoexcel(outputs,output+"_output"+str(avgregion))	
	gentools.transposecsv(output+"_output"+str(avgregion))

# Function that loops through the annotation GFF and averages reads into a metagene plot.
# Variables: 
# counts_filestring is the input density file
# GFFgen_filename is the GFF file name
# utrgfffilename,utr5gfffilename are filenames for the UTR annotations.
# avgregion determines the feature of interest: 0 for transcript start, 1 for start codon, 2 for stop codon, and 3 for start of polyA tail.
# regionlength_5,3 are the distance upstream and downstream of the feature of interest.
# ignoreutr when = 1 ignores UTR annotation (manual mode). 0 for annotated.
# shift for ribosome P sites, typically 12 or 13 if 5' aligned reads.
# eliminate is 0 for all genes, 1 for genes not rejected by givegene (including overlaps), 2 for genes not rejected by givegene but includes overlaps.
# thresh is the lowest acceptable value for rpkm in a mainORF to be included in the average.
# equalweight =1 means all genes in the average have equal weight (normalization done). For avgregion1,2 this is the ORF region only (i.e. the region set by goodzone).
# goodzone is a variable that excludes part of the ORF region for normalization to eliminate effects of start codon or stop codon peaks.
# In Dec 2015, added capability to filter out genes with random spikes. Currently set so nothing removed unless user changes code. However, it does report very high spikes for some cases.
# This function was completely revamped on August 8, 2020.
def makeavggene(counts_filestring,GFFgen_filename,utr5gfffilename,utr3gfffilename,regionlength5,regionlength3,ignoreutr,thresh,shift,eliminate,goodzone,equalweight,avgregion):
	GFFgen=GFF.parse(GFFgen_filename)
	GFFlist=seqtools.makeGFFlist(GFFgen)
	counts0=builddense.readcountsf(counts_filestring+"_plus_")		
	counts1=builddense.readcountsf(counts_filestring+"_minus_")
	counts=[counts0,counts1]
	utrgffgen=GFF.parse(utr5gfffilename)
	utrtable5=seqtools.makeutrtable(utrgffgen)
	utrgffgen=GFF.parse(utr3gfffilename)
	utrtable3=seqtools.makeutrtable(utrgffgen)
	missedthresh=0
	illegalgenes=0
	genesinlist=0
	tooshortlist=0
	averagegene=[0 for x in range(0,(regionlength5+regionlength3))]

	if ignoreutr==1 and (avgregion==0 or avgregion==3):
		print("UTR annotations must be used for averages at transcript ends.")
		exit()

	if ignoreutr==1:
		bp5_0=-regionlength5
		bp3_0=-regionlength3
	else:	
		bp5_0=0
		bp3_0=0
		
	# Need extensions if going over end of annotations for transcript start/stops. And no annotation for end not being examined.
	if avgregion==0:
		bp5_0+=regionlength5
		bp3_0="none"
	elif avgregion==3:
		bp3_0+=regionlength3
		bp5_0="none"
	elif avgregion==1:
		bp3_0="none"
	elif avgregion==2:
		bp5_0="none"	
	
	# Check goodzones to make sure they fit.
	if (regionlength3<=goodzone and avgregion==1) or (regionlength5<=goodzone and avgregion==2):
		print("Bad Goodzone. Doesn't fit.")
		exit()
	
	for chrom in GFFlist:
		feat_num=-1	
		for feature in GFFlist[chrom].features:
			feat_num+=1
			gg=seqtools.givegene(chrom,feat_num,[GFFlist,utrtable5,utrtable3],counts,[bp5_0,bp3_0,shift],eliminate)
			genecounts=gg[0]
			# Get rid of dubious genes, nongenes, genes with overlap of others.
			if (genecounts ==-1 or genecounts==-2):
				illegalgenes+=1
				continue

			# Continue if the UTR annotation isn't there (and we aren't in manual). This eliminates genes that don't have appropriate annotation (5' for avgregion0,1 or 3' for avgregion2,3.
			if ignoreutr!=1:
				if avgregion==0 or avgregion==1:
					if (gg[5]==bp5_0) and ignoreutr!=1:
						illegalgenes+=1
						continue
				else:
					if (gg[6]==bp3_0) and ignoreutr!=1:
						illegalgenes+=1
						continue
				
			# Define ORF
			bp5=gg[5]
			bp3=gg[6]
			start=bp5
			end=len(genecounts)-bp3

			############# Extra code that could be called.
			# Filter sequences for those of interest - specific sequences, or go term. Need to call supplementary code to be edited.
			# import listavg_extra
			# filterresult=listavg_extra.filterXXXX()
			#if filterresult==0:
			#	continue

			# Threshold on total gene counts: 
			totalgenereads=float(1000)*sum(genecounts[start:end])/len(genecounts[start:end])
			if totalgenereads<thresh:
				missedthresh+=1
				continue
			
			# Check if annotated UTRs or mainORFs are too short.
			if avgregion==0:		# Align on 5' ends.
				if regionlength3>(bp5-regionlength5):
					tooshortlist+=1
					continue
				else:
					genesinlist+=1
					countlist=genecounts[0:(regionlength5+regionlength3)]
					if equalweight==1:
						totcounts=sum(countlist)
					else:
						totcounts=1

			elif avgregion==1:	# Align at start codon.
				if regionlength5>bp5 or regionlength3>(end-start):
					tooshortlist+=1
					continue
				else:
					genesinlist+=1
					countlist=genecounts[bp5-regionlength5:(bp5+regionlength3)]			
					if equalweight==1:
						totcounts=sum(countlist[(regionlength5+goodzone):])
					else:
						totcounts=1
	
			elif avgregion==2:	#align at stop codon.
				if regionlength3>bp3 or regionlength5>(end-start):
					tooshortlist+=1
					continue
				else:
					genesinlist+=1
					countlist=genecounts[-(regionlength5+bp3):len(genecounts)-(bp3-regionlength3)]
					if equalweight==1:
						totcounts=sum(countlist[0:(regionlength5-goodzone)])
					else:
						totcounts=1

			elif avgregion==3:	# Align at 3' end.
				if regionlength5>(bp3-regionlength3):			
					tooshortlist+=1
					continue
				else:
					genesinlist+=1
					countlist=genecounts[-(regionlength5+regionlength3):]
					if equalweight==1:
						totcounts=sum(countlist)
					else:
						totcounts=1
			
			# Spike checking and averaging.
			totalcountspike=sum(countlist)
			if totalcountspike>0:
				if (sum(sorted(countlist)[-2:])/totalcountspike)>10:	# Currently the spike checker is disabled with this 10 value. 0.8 can be used.
					genesinlist-=1
					print("Spike removed. Gene "+feature.id)
					continue
			if totcounts!=0:	
				for i in range(len(countlist)): 
					if (countlist[i]/totcounts>100 and i<(regionlength5+goodzone) and avgregion==1) or (countlist[i]/totcounts>100 and (len(countlist)-i)<=(regionlength3+goodzone) and avgregion==2):
						print("Potential spike on avg"+str(avgregion)+" gene "+feature.id)
						
					averagegene[i]+=countlist[i]/totcounts	##### Do the averaging.

	for m in range(len(averagegene)):			# Divide by total number of genes.
		if genesinlist!=0:
			averagegene[m]/=genesinlist
		else:
			print("Error, no genes to average.")		

	print("Genes under threshold: "+str(missedthresh))
	print("Genes removed by givegene/UTR annotation concern: "+str(illegalgenes))
	print("Genes removed because feature(s) too short: "+str(tooshortlist))
	print("Genes in average: "+str(genesinlist))
	
	return averagegene








			
