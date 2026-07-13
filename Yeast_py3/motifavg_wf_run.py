import searchcode
import csv
from datetime import datetime
import os
thiscode=open("/Users/meydanmarksf2/Documents/dscode/GuydoshLabCode1.11/motifavg_wf_run.py")
codeprint=thiscode.read()
thiscode.close()

# motifavg is a very useful function that goes through the entire transcriptome and, depending on your input, computes an enrichment score
# at every amino acid, codon, nucleotide triplet, etc. 
# This is still very experimental - a lot of filters and cutoffs are hardcoded so you really need to be familiar with the underlying code in searchcode to use this effectively.
# Note that overlap checks not performed - can be changed in the hard code. Goodgene is 2.
# Note the old version of this is available in the "_extra" python file for searchcode.

# inframe = 0 is for nt positions (ie looking for sequence effects on ribosome entrance channel), 
# 1 is for main reading frame (i.e. triplets or codons), inframe 2 is for amino acids (translated.
inframe=str(2)  

motiflist=str(3)	# Just the size of the motif here, in units of bp or aa as appropriate to search.

mismatches=str(0)	#mismatches is not longer used - so doesn't matter what's here.

# This is the number of nt on either side of the position of interest to include when computing the pause scores (new or old) reported in the csv.
scorewin=str(0) 

avgwin=str([50,50])	#For the average window used in computing the binary output. List of 2 for up and downstream. Can enter 0 to skip the binary average output. This also skips new pause score since it is derived from it.
#avgwin=str([10,10])
#avgwin=str([0,0])

# This is the ORF threshold to be included in the old pause score. I believe new pause score is not thresholded.
thresh=str(10) 

# These are the shift values to use. Each counts file will have this applied. 
# Here things get a bit tricky because the motifs that are used are usually keyed to the first base. so that means for P site, you will probably want 12 most of the time and 15 for A site. 
# However, some datasets are underdigested, in which case you may want to use 13 and 16. For 3' end alignment, other factors should be considered.
shifts=[-18]

hardshift="13"	# This is the value where I think ribosome reads normally begin. Use care in choosing this.

gfffile="/Users/meydanmarksf2/Documents/dscode/yeast.gff"
gfffile5="/Users/meydanmarksf2/Documents/dscode/yeast_5UTRc.gff"
gfffile3="/Users/meydanmarksf2/Documents/dscode/yeast_3UTRc.gff"

# Say which type is being used.
#typeregion="UTR5"
typeregion="ORF"
#typeregion="UTR3"

# The stuff below is there to do multiple density files.

countsfiles= ["/Users/meydanmarksf2/Documents/dscode/density/end3SM061F/end3SM061F"]

#countsfiles=["/Volumes/pRAID/dsdata4/density/23F/23F"]
#countsfiles+=["/Volumes/pRAID/dsdatanew2/density/45F/45F","/Volumes/pRAID/dsdatanew2/density/46F/46F","/Volumes/pRAID/dsdatanew2/density/47F/47F","/Volumes/pRAID/dsdatanew2/density/48F/48F"]

# Here is where you setup your outfile strings.
outfilestring0s=[]
for filename in countsfiles:
	tpath="/Users/meydanmarksf2/Documents/Desktop/motifavg/"+filename.split("/")[-1]
	outfilestring0s.append(tpath+"/"+filename.split("/")[-1])
	if not os.path.exists(tpath):
		os.makedirs(tpath)

print("Operation begun at "+str(datetime.now()))

for i in shifts:
	shift=str(i)
	for (outs,countsf) in zip(outfilestring0s,countsfiles):
		outfilestring=outs+"_shift"+str(i)
		countsfile=countsf
		print("")
		commandstring="searchcode.motifavg_2_wf('"+gfffile+"','"+gfffile5+"','"+gfffile3+"','"+countsfile+"',"+motiflist+","+inframe+","+thresh+",'"+outfilestring+"',"+mismatches+","+shift+","+scorewin+","+avgwin+","+hardshift+",'"+typeregion+"')"
		print("Command is:")
		print(commandstring)
		eval(commandstring)	

print("Operation concluded at "+str(datetime.now()))



print("Command file is:")
print("*************************************")
print(codeprint)
print()
print("*************************************")
print("End command file.")
print()

