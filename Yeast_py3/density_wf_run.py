import builddense
from datetime import datetime
import os

thiscode=open("/Users/meydanmarksf2/Documents/dscode/GuydoshLabCode1.11/density_wf_run.py")
printcode=thiscode.read()
thiscode.close()

#### Do not edit above this line #######################################################
########################################################################################
# This script will tally the ribosome density in to 1) density files used by all the other scripts and 2) wig files for use in a genome browser. 
# Density files and wig files are interconvertible with another by using another script.

# Here is the rootpath where you are working.
rootpath="/Users/meydanmarksf2/Documents/dscode/GuydoshLabCode1.11/"

# Here are your detailed path settings.
### Be careful as to whether you include /coding/ or /codingandpolyAcombined/ here. Big difference depending on application.
# Typically for 3' end alignment, we use coding. But could be either for 5' end alignment, though usually coding as well.
#path1=rootpath+"/alignments/codingandpolyAcombined/chrom/"
#path2=rootpath+"/alignments/codingandpolyAcombined/splice/"
path1=rootpath+"/coding/chrom/"
path2=rootpath+"/coding/splice/"
#path2="-1" # For chromosome alignment (no splice) only. Use -1 in path1 for splice alignment only.

# Your file prefixes coming from your bowtie script.
#files=['DJY0910_1','DJY0910_2','DJY0708_1','DJY0708_2','DJY0708_3','DJY0708_4','68F69F','70F71F']
files=['SM128Fd','SM129Fd','SM130Fd','SM135Fd','SM136Fd','SM137Fd']

# Variables for doing manual or automatic normalization for each file (i.e. the per million in rpm).
# Putting in -1 for  makes sure the program computes the count itself for each file.
# Use actual values for each file you have an external normalization, ie for plasmid alignments.
#totcounts=['-1','-1','-1','-1','-1','-1','-1','-1']
totcounts=['-1', '-1', '-1', '-1', '-1', '-1', '1','-1', '-1', '-1', '-1', '-1', '-1', '1','-1', '-1', '-1', '-1', '-1', '-1', '1']

# The variable assignmode determines whether the read assignment is by 5' end, 3' end, coverage, or flipped coverage. (Flipped coverage is used for some kits). Options:
#assignmode="end5"
#assignmode="end3"
assignmode="end3"
#assignmode="covflip"

# These variables allow you to filter for certain read lengths. These are the smallest and largest (inclusive) read lengths that will be written out.
# Use caution here:
#	If you are doing RNA-Seq and the expected fragment size is > the readlength, these should be 0 and readlength, ie 0, 50 or 0, 100.
# 	If you ran a poly(A) removal previously in the pipeline, these should be the full range (or bigger), ie 0, 100.
# 	Usually, there is size selection earlier in cutadapt, so it may not be needed here.
smallestsize=57
largestsize=63

# Note that the prefix will be automatically set according to the assignment mode. 
# If you want to override it, do so below by changing auto to what you want. Generally, we don't do this but it might be helpful if you are doing something experimental. 
prefix="auto"
#prefix="put_manual_prefix_here"

# This variable set the species - S cerevisiae, S pombe, or E coli.
species="scer"
#species="pombe"
#species="coli"

##### DO not change below this line. Everything below is 100% automated. ###############
########################################################################################





# This sets the juncrange parameter automatically. 
# The value 99 in the command corresponds to the size of the splice junctions used in the splice alignment (juncrange). This is always 99 from July 2020 on. This means read lengths should never be >100 nt.
# However, for pombe, this should be 98, as it signifies something different.
if species == "scer":
	juncrange=str(99)
elif species == "pombe":
	juncrange=str(98)
elif species == "coli":
	juncrange=str(0)
else:
	print("Invalid species.")
	exit()
	
# This sets the auto prefix.
if prefix=="auto":
	prefix=assignmode
	suffix=""
	if assignmode=="covflip":
		prefix="cov"
		
smallestsize=str(smallestsize)
largestsize=str(largestsize)
		
wigpath=rootpath+"/wigfiles/"
if not os.path.exists(wigpath):
		os.makedirs(wigpath)
print("Operation begun at "+str(datetime.now()))
for i in range(len(files)):
	fname=files[i]
	cnts=totcounts[i]
	outpath=rootpath+"/density/"+prefix+fname+suffix+"/"
	if not os.path.exists(outpath):
		os.makedirs(outpath)
	#print ""
	commandstring="builddense.setdense_wf2('/Users/meydanmarksf2/Documents/dscode/yeast.gff','"+path1+fname+"_match.SAM','"+path2+fname+"_match.SAM',"+juncrange+",'"+outpath+prefix+fname+suffix+"','"+wigpath+prefix+fname+suffix+"',"+smallestsize+","+largestsize+",'"+species+"','"+assignmode+"',"+cnts+")"
	
	print("Command is:")
	print(commandstring)
	eval(commandstring)	

print("Operation concluded at "+str(datetime.now()))

# This appends this input file to the output so you have it for record keeping.
print("Command file is:")
print("*************************************")
print(printcode)
print()
print("*************************************")
print("End command file.")
print()
