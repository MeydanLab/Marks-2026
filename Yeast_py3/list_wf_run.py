import listavg
from datetime import datetime
from BCBio import GFF
import os

# Read in this file so it will be dumped out in case you don't have a good system for keeping track of your input parameters.
thiscode=open("/Users/marksjp/Documents/Scripts/GuydoshLabCode1.11/Yeast_py3/list_wf_run.py")
commandfile=thiscode.read()
thiscode.close()

####### DO not change above this line. ##################################################
#########################################################################################
# Input filenames and paths for up to 3 sets of files in different folders.
pathroot3="/Users/marksjp/Documents/Scripts/GuydoshLabCode1.11/Yeast_py3/density"
pathroot2="/Users/marksjp/Documents/Scripts/GuydoshLabCode1.11/Yeast_py3/density"
pathroot="/Users/marksjp/Documents/Scripts/GuydoshLabCode1.11/Yeast_py3/density"
files3=[]
files2=[]
files=['end5SM015Fdx5763','end5SM016Fdx5763','end5SM023Fdx5763','end5SM024Fdx5763']

# Input transcriptome sequence and annotations
gfffile="/Users/marksjp/Documents/Scripts/GuydoshLabCode1.11/yeast.gff"
utrgfffile="/Users/marksjp/Documents/Scripts/GuydoshLabCode1.11/cherry_UTRs_for_R64.1.1/yeast_3UTRc.gff"
utr5gfffile="/Users/marksjp/Documents/Scripts/GuydoshLabCode1.11/cherry_UTRs_for_R64.1.1/yeast_5UTRc.gff"

# Path to output files
outpath="/Users/marksjp/Documents/Scripts/test_output/"

# Output filenames prefix
outprefix="test"

# The shift value used to compute reads on genes. Typically 12 or 13 for 5' mapped reads.
# For mRNA-Seq, the shift should be 0 unless a strict control for profiling.
shift=-15
#shift=0

# This variable controls which genes are eliminated.
# All included = 0. No overlaps, dubious = 1. No dubious, overlaps OK = 2.
# Note that overlaps are determined by using the annotated UTR information on the neighboring genes, and whatever UTR information is below for the gene of interest.
eliminate=0

# This determines whether pause scores at start/stop are normalized to main ORF levels. 1 to normalize, 0 to not normalize. Generally 1 for 80S/disome, 0 for 40S.
normalizepause=1

# These variables allow you to ignore regions at start/end of main ORFs to avoid end artifacts. For example, start_special=15, would exclude 15 nt at start of ORFs in the count. Generally we keep these at 0 unless something special being done. 
start_offset=0
end_offset=0

# These nt offsets are here for pause scores because the peaks may not be where the P site shift is (usually 13 nt).
# Usually we put 0 for start, 5 for termination, and 35 for 40S stack or 31 for disome stack - if 5' aligned.
startpauseoffset=0
termpauseoffset=5
termstackpauseoffset=35

# This option allows you to ignore the UTR annotation and simply count reads a fixed distance up/down stream of the ORF.
# Typically we don't use this.
ignoreutr5=0
ignoreutr3=0

# If the annotated UTR information is not used, then this value will be the length of the UTR.
# Or, if annotation is used, this is the amount to go beyond upstream end of annotated 5'UTR or downstream of annotated 3'UTR. Typically we don't do this.
# One application is we sometimes add 25 nt onto the annotated end in order to make sure we got all the 3'UTR reads.
bp5=0
bp3=0

# This variable sets the number of nt downstream and upstream of pause peaks to include in pause score. Usually, 0, 1, or 2.
pausehalfregion=2

# This variable is 1 normally, but if you have coverage data, then it must be the readlength
covlength=1
#covlength=100

# Variable for future use. Currently does nothing
future1=0

# Variable for future use. Currently does nothing
future2=0

####### DO not change below this line. ##################################################
#########################################################################################

shift=str(shift)
eliminate=str(eliminate)
normalizepause=str(normalizepause)
start_offset=str(start_offset)
end_offset=str(end_offset)
startpauseoffset=str(startpauseoffset)
termpauseoffset=str(termpauseoffset)
termstackpauseoffset=str(termstackpauseoffset)
bp5=str(bp5)
bp3=str(bp3)
ignoreutr5=str(ignoreutr5)
ignoreutr3=str(ignoreutr3)
pausehalfregion=str(pausehalfregion)
covlength=str(covlength)
future1=str(future1)
future2=str(future2)

print("Operation begun at "+str(datetime.now()))

filelist=[]
for filename in files:
	filelist.append(pathroot+"/"+filename+"/"+filename)
for filename in files2:
	filelist.append(pathroot2+"/"+filename+"/"+filename)
for filename in files3:
	filelist.append(pathroot3+"/"+filename+"/"+filename)
filelist=str(filelist)

output=outpath+outprefix

print("")
commandstring="listavg.totalquant_wf("+filelist+",'"+output+"','"+gfffile+"','"+utrgfffile+"','"+utr5gfffile+"',"+bp5+","+bp3+","+ignoreutr5+","+ignoreutr3+","+shift+","+normalizepause+","+eliminate+","+start_offset+","+end_offset+","+startpauseoffset+","+termpauseoffset+","+termstackpauseoffset+","+pausehalfregion+","+covlength+","+future1+","+future2+")"
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
