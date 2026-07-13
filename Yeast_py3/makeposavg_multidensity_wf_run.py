import searchcode
from datetime import datetime
import os


######### Do not change above line.
# This is the script that will average ribosome profiling data for a set region around a list of sites.

thiscode=open("/Users/meydanmarksf2/Documents/dscode/GuydoshLabCode1.11/makeposavg_multidensity_wf_run.py")
filecontents=thiscode.read()
thiscode.close()

GFFfile="/Users/meydanmarksf2/Documents/dscode/yeast.gff"
utrgfffile5="/Users/meydanmarksf2/Documents/dscode/yeast_5UTRc.gff"
utrgfffile3="/Users/meydanmarksf2/Documents/dscode/yeast_3UTRc.gff"	

densityfiles=[
"/Users/meydanmarksf2/Documents/dscode/density/end3SM128Fd/end3SM128Fd",
"/Users/meydanmarksf2/Documents/dscode/density/end3SM130Fd/end3SM130Fd",
"/Users/meydanmarksf2/Documents/dscode/density/end3SM135Fd/end3SM135Fd",
"/Users/meydanmarksf2/Documents/dscode/density/end3SM137Fd/end3SM137Fd"]

#densityfiles=['/Users/guydoshnr/pRAID_copy/dsdatalab//density/end3SM015F23F/end3SM015F23F','/Users/guydoshnr/pRAID_copy/dsdatalab//density/end3SM015Fd23FdX5763/end3SM015Fd23FdX5763']

riboshift=0	# Enter a shift value.

# This is the region to include in the average upstream, downstream, in nt. 
# If seqwin is negative, it allows "spillover" into other areas. For example, if you want to include parts of the gene outside the region (i.e. ORF reads just upstream of stop in a 3'UTR average), then make this value(s) negative.
seqwin=str([50,50])		

# 5'UTR = 0, ORF = 1, 3'UTR = 2. -2 is special for preAug2020 3'UTR lists that used ORF coordinates.
region = 1

# Set to 0 for equal weighting, 1 for ORF weighting, 2 for no additional weighting.
normmethod= 1

# ORFthresh is minimal threshold of the main ORF required, in rpkm units.
ORFthresh = 0

# thresh is the number of rpm counts required in the seqwin window of interest. Usually -1.
thresh = -1

# These are the sites you are averaging at.
genelist="/Users/meydanmarksf2/Documents/dscode/csvfiles/KKK_noconsec30.csv"
#genelist="/Users/guydoshnr/dsa/dsdatagen2/paper_oldnotused/posavg/csvfiles/polyK3.csv"

# Where your output csv will go and prefix for the outputfile.
outfilestring="/Users/meydanmarksf2/Documents/Disome_project/GL134/polyk3"

# IF you made your input list with nonstandard UTR settings, include them here. Typically these are 0 (so annotated UTRs).
bp5=0
bp3=0


####### Do not change below line.
region=str(region)
normmethod=str(normmethod)
ORFthresh=str(ORFthresh)
thresh=str(thresh)
densityfiles=str(densityfiles)
riboshift=str(riboshift)
bp5=str(bp5)
bp3=str(bp3)
print("")
print("Operation begun at "+str(datetime.now()))
commandstring="searchcode.makeposavg_wf0('"+genelist+"','"+GFFfile+"','"+utrgfffile5+"','"+utrgfffile3+"',"+seqwin+","+densityfiles+",'"+outfilestring+"',"+riboshift+","+thresh+","+ORFthresh+","+normmethod+","+region+","+bp5+","+bp3+")"
print("Command is:")
print(commandstring)
eval(commandstring)	

print("Operation concluded at "+str(datetime.now()))


print("Command file is:")
print("*************************************")
print(filecontents)


print()
print("*************************************")
print("End command file.")
print()


