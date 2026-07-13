import searchcode_extra
thiscode=open("/Users/meydanmarksf2/Documents/dscode/GuydoshLabCode1.11/findpauseratio_wf_run.py")
thiscodeprint=thiscode.read()
thiscode.close()

# This is a function that goes through the entire transcriptome of 2 datasets and computes for every nucleotide position:
# The pause score in both data sets (counts at read divided by avg in gene) and the ratio between reads datasets.
# This is a ribosome stall finding tool.

# Give the 2 samples that will be compared in the analysis.
#countsfile1="/Users/guydoshnr/Desktop/pombe_temp/dsdatapombe/density/end30mm64FsX1518/end30mm64FsX1518" # This is the numerator in the ratio
#countsfile="/Users/guydoshnr/Desktop/pombe_temp/dsdatapombe/density/end30mm63FsX1518/end30mm63FsX1518" # This is the denominator in the ratio
countsfile1="/Users/meydanmarksf2/Documents/dscode/density/end3SM062F72F/end3SM062F72F"
countsfile="/Users/meydanmarksf2/Documents/dscode/density/end3SM066F76F/end3SM066F76F"


# There is an option to offer "control" data to first weed out positions that do not meet some minimal ratio criterion in an alternative dataset.
# One example where this would be used is the case where you were looking at the ratio between disome and short-read data. Before checking these (typically noisy and biased) datasets, 
# you might first want to make sure that there at least some level reatio already there for a dom34KO to WT comparison. So dom34KO and WT might be included here.
# Or maybe you'd do like this: first screen for positions that had some enrichment in a dom34KO strain vs WT and then go check the situation in short-read or disome data vs 28-nt data.

# Set to -1 if not doing these.
#countsfile1_c="/Users/guydoshnr/Desktop/pombe_temp/dsdatapombe/density/0mm63FsX1518/0mm63FsX1518"
#countsfile_c="/Users/guydoshnr/Desktop/pombe_temp/dsdatapombe/density/0mm63F/0mm63F"
countsfile_c="-1"
countsfile1_c="-1"
# This is the threshold of the screen data - not needed if countsfile_c files above not provided.
ratio_cthresh="0.3"			

separation="0"	# Used when you want to look for reads shifted by different amounts in 2 samples, like 28mers to disome reads, for example. If both datasets are 28-nt footprints, then put 0 here.

#outfilestring="/Users/guydoshnr/Desktop/pombe_temp/pauseratio/64Fs_63Fs_end3_noratiothresh_lowcountthresh"
outfilestring="/Users/meydanmarksf2/Documents/gustavolab/findpauseratio/WTperoxidevsrad6peroxide_2"
GFFfilepath="/Users/meydanmarksf2/Documents/dscode/yeast.gff"
UTR5GFFfilepath="/Users/meydanmarksf2/Documents/dscode/yeast_5UTRc.gff"
UTR3GFFfilepath="/Users/meydanmarksf2/Documents/dscode/yeast_3UTRc.gff"
#GFFfilepath="/Users/guydoshnr/macpro/genomics/pombegenome/pombe.gff"
#UTR3GFFfilepath="/Users/guydoshnr/macpro/genomics/pombegenome/pombe_UTR3.txt"
#UTR5GFFfilepath="/Users/guydoshnr/macpro/genomics/pombegenome/pombe_UTR5.txt"

# Typically you don't need to use this because shifting can be done later in the csv output but if you want, you can shift the data for A site, P site, etc.
shift="0"

# This the number of nt on either side of the nt to include in computation of ratio and pause score (essentially its the width of a boxcar filter).
# Put 0 if you are looking at single nt positions. Put 1 if you want a 3-nt window to be used (1nt on either side).
halfwindow=str(0)

# Here is a threshold for each of the pause scores (for countsfile and countsfile1, respectively). Both thresholds need to be met to be included in the dataset.
# You can put in 0 so everything is included if you like.
pausethresh=str([0,0])

# This ratio threshold is how high the ratio has to be to be included. 
ratiothresh="0.8"
#ratiothresh="10"

# Counts thresh is how many counts you need in each dataset for a position to be considered. I will often put the value for some finite number of reads, say 3 reads or 5 reads, depending on application. Note that 1 read is usually around 0.01 - 0.1 rpm in decently deep dataset.
countsthresh=str([0.1,0.1])	
#medianthresh=str([.26,.33])


# This is a threshold for how many minimal counts should be in the entire gene (i.e. the denominator of the pause score). I usually just put a small, nonzero value here to essentially include all genes with any expression though you might want to go higher.
# Note that we no longer use the median and this is actually a mean.
medianthresh=str([.000001,.000001])	# For median of gene. (Often good to just put something nonzero.)



print("")
commandstring="searchcode_extra.findpauseratio_wf('"+outfilestring+"','"+countsfile+"','"+countsfile1+"','"+countsfile_c+"','"+countsfile1_c+"','"+GFFfilepath+"','"+UTR5GFFfilepath+"','"+UTR3GFFfilepath+"',"+shift+","+halfwindow+","+separation+","+ratiothresh+","+countsthresh+","+pausethresh+","+medianthresh+","+ratio_cthresh+")"
print("Command is:")
print(commandstring)
eval(commandstring)	





print("Command file is:")
print("*************************************")
print(thiscodeprint)
print()
print("*************************************")
print("End command file.")
print()

