import seqtools
import os
from datetime import datetime

thiscode=open("/Users/meydanmarksf2/Documents/dscode/GuydoshLabCode1.11/writegene2_wf_run.py")
thiscodetext=thiscode.read()
thiscode.close()

# This code an option for generating write genes of many different input files with varied extension for footprint length.
#densityfiles=[]
#for i in range(15,81):
#	densityfiles.append("/Users/guydoshnr/dsa/dsdatalab/sc/40S/sizerangeend5/density/end5_"+str(i)+"_DY182Ffx/end5_"+str(i)+"_DY182Ffx")
#densityfiles=str(densityfiles)

############## DO NOT CHANGE ABOVE THIS LINE.################

# Use this code normally to list your density files.
densityfiles=str(["/Users/meydanmarksf2/Documents/dscode/density/end3SM061F71F/end3SM061F71F",
"/Users/meydanmarksf2/Documents/dscode/density/end3SM062F72F/end3SM062F72F",
"/Users/meydanmarksf2/Documents/dscode/density/end3SM065F75F/end3SM065F75F",
"/Users/meydanmarksf2/Documents/dscode/density/end3SM066F76F/end3SM066F76F"])

# What genes you want:
feature=str(["GCN4"])		# Formal or informal names are okay.

gfffile="/Users/meydanmarksf2/Documents/dscode/yeast.gff"
utrgfffile3="/Users/meydanmarksf2/Documents/dscode/yeast_3UTRc.gff"
utrgfffile5="/Users/meydanmarksf2/Documents/dscode/yeast_5UTRc.gff"	

outfilepath="/Users/meydanmarksf2/Documents/gustavolab/writegene/output"

# Here is how much UTR to put on each gene. We used to call these variables "shift" values.
bp5=100
bp3=100

# If this value is set to 1, then the UTRs will be the exact lengths given in the bp values. 
# If set to 0, then the bp values will extend the annotated UTRs. 
# If you don't want UTRs at all, put -1.
manualutrs=1
#manualutrs=0
#manualutrs=-1

# Options for intron:
#intron=0; standard way with spliced transcript returned. (Sense direction, 5' to 3')
#intron=1; as before, unspliced transcript is returned. (Sense direction, 5' to 3')
#intron=2; unspliced transcript returned in 5' to 3' direction. Unspliced transcript returned on antisense strand in 5' to 3' direction.
#intron=3; same as intron = 2 except the orientation for both sense and antisense will be as they appear on the chromosome. Nothing is flipped. So one will be 5' to 3' and the other 
intron="0"		

############## DO NOT CHANGE BELOW THIS LINE.##################

if manualutrs==1:
	bp5=str(-bp5)
	bp3=str(-bp3)
elif manualutrs==0:
	bp5=str(bp5)
	bp3=str(bp3)
elif manualutrs==-1:
	bp5="'none'"
	bp3="'none'"
else:
	print("Error in manualutrs")
	exit()
	
print("Operation begun at "+str(datetime.now()))

commandstring="seqtools.writegene2_wf("+bp5+","+bp3+","+intron+","+densityfiles+","+feature+",'"+gfffile+"','"+utrgfffile5+"','"+utrgfffile3+"','"+outfilepath+"')"
print("Command is:")
print(commandstring)
eval(commandstring)	

file=open(outfilepath+"_output.txt","w")
file.write("Command is:\r\n")
file.write(commandstring)
file.write("\r\n\r\n")
file.write("Code in script:\r\n")
file.write(thiscodetext)


print("Operation concluded at "+str(datetime.now()))

