import genometools
from datetime import datetime


# Read in this file so it will be dumped out in case you don't have a good system for keeping track of your input parameters.
thiscode=open("/Users/marksjp/Documents/Scripts/GuydoshLabCode1.11/Yeast_py3/findseqinstance2_run_wf.py")
commandfile=thiscode.read()
thiscode.close()


#### DO not write above this line ####

# This script is for searching the yeast genome for a particular sequence motif and listing all occurrences of it.
# Note this program returns positions in 0-based coords.
# Note all mRNA positions referenced to start of region of interest (5'UTR, ORF, or 3'UTR), unlike the old findseqintance where 3'UTRs were referenced to main ORF start codon.
mlist=['_IP']
for x in mlist:
	outfilestring="/Users/marksjp/Documents/Scripts/test_output/test"+x
	gfffile="/Users/marksjp/Documents/Scripts/GuydoshLabCode1.11/yeast.gff"
	utrgfffile5="/Users/marksjp/Documents/Scripts/GuydoshLabCode1.11/cherry_UTRs_for_R64.1.1/yeast_5UTRc.gff"
	utrgfffile3="/Users/marksjp/Documents/Scripts/GuydoshLabCode1.11/cherry_UTRs_for_R64.1.1/yeast_3UTRc.gff"

# regionwindow is either -1 to skip it, or a list of 2 positions that demarcate a window within the region that is included in the search. For example, [0,5] is the first 5 nt only of the region.
	regionwindow="-1"
#regionwindow=[65,69]

	region=1		# 0 for 5'UTR, 1 for ORF, 2 for 3'UTR. 

	findstop=0 # Put 0 to skip this. Give the next frame stop downstream of motif that is no more than findstop nt away (short only). Put negative value to eliminate the short ones (long only).

# This filter prevents you from getting closely-spaced motifs. For example, if you are looking at single K motifs, you might not want to consider each K in the sequence KKKK as part of this.
# So put in the number of nucleotides that motifs should be separated by or just a really big number to take the first in each case.
# A negative value will give cases where there is only 1 per gene.
# Put 0 to keep everything.
	noconsec=0

	stopcodon="0"		# this filters out genes with other stop codons. Put in "0" to take all stops. Put in stop codon otherwise.
#stopcodon="TAA"

# inframe is a variable that tells you where to look. 0 is for nucleotide sequences in all frames, 1 is for nucleotide sequences only in the 0-frame, and 2 is for amino acid sequences.
# Use 2.1 for +1 frame and 2.2 for +2 frame aa searches.
	inframe=2		

# This is the motif you are looking for.
# Include a "_" character for having a space in the sequence. However this is not compatible with mismatches.
# Enter "PENULT" to get all the ORF penultimate codon matches in a 3'UTR.
# If motif is a filename, then the sites in that list will be filtered with the consec, stop codon, and region filters below. No new sites are added. The list should be sorted first so consecutive genes come one after another.
	motif=x
#motif='PENULT'
#motif='M'
#motif='G_S'
#motif="/Users/nick/genomics/loci of interesting sequences/UTR3/Stop_allframes_noconsec_sorted.csv"

# Number of mismatches to tolerate. 
	mismatches=0

######## Do not write below ############

	regionwindow=str(regionwindow)
	region=str(region)
	findstop=str(findstop)
	noconsec=str(noconsec)
	inframe=str(inframe)
	mismatches=str(mismatches)

	print("Operation begun at "+str(datetime.now()))

	print("")
	commandstring="genometools.findseqinstance2_wf('"+gfffile+"','"+utrgfffile5+"','"+utrgfffile3+"','"+motif+"',"+inframe+",'"+outfilestring+"',"+mismatches+","+noconsec+",'"+stopcodon+"',"+region+","+regionwindow+","+findstop+")"  


	print("Command is:")
	print(commandstring)
	eval(commandstring)	
 
 

	print("Operation concluded at "+str(datetime.now()))
	print("Command file is:")
	print("*************************************")
	print(commandfile)
	print()
	print("*************************************")
	print("End command file.")
	print()
