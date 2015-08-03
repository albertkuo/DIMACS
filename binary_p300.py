# Albert Kuo, 2015
# The following code adds p300 feature as a new column to each chromosome file
# Note: It does not handle chromM or chromX
import glob
import os
import pandas as pd
import numpy as np
from random import randint, seed
from quicksect import IntervalNode

# read p300 data
p300 = pd.read_table('data/GM12878_distal_P300_align.hg19.txt.gz', compression = 'gzip', header=None)
p300.columns = ['chr','lower']
p300['upper'] = p300['lower']+200
p300['chr'] = p300['chr'].str[3:] # Get chr number only

def get_data(filename, chrom):
	dat = pd.read_table(filename,skiprows=1)
	nrows = dat.shape[0]
	dat['chr'] = chrom
	dat['lower'] = np.arange(1,200*(nrows),200)
	dat['upper'] = np.arange(200,200*(nrows+1),200).reshape(nrows,1)
	return dat

files = glob.glob('Spectacle/SAMPLEDATA_HG19_NEW/*.txt')  # create the list of files
for filename in files:
	## read file
	print filename
	chrom = os.path.split(filename)[1].split("_")[1][3:] # get chromosome number from file name
	dat = get_data(filename, chrom)
	chr_p300 = p300[p300['chr']==chrom]  # Select only rows with corresponding chromosome number

	## match intervals
	query = zip(dat['lower'],dat['upper'])
	data = zip(chr_p300['lower'],chr_p300['upper'])

	# Modified code from https://www.biostars.org/p/99/
	def find(start, end, tree):
		"Returns a list with the overlapping intervals"
		out = []
		tree.intersect( start, end, lambda x: out.append(x) )
		#return True if there is an intersection
		return not not out

	# start the root at the first element
	start, end = data[0]
	tree = IntervalNode( start, end )

	# build an interval tree from the rest of the data
	for start, end in data[1:]:
		tree = tree.insert( start, end )

	overlap = []
	for start, end in query:
		overlap.append(find(start, end , tree))

	dat['p300'] = map(int, overlap)
	dat = dat.ix[:, [0,1,2,3,4,5,6,7,8,12]]
	print dat['p300'].value_counts()

	## write file
	new_filename = "GM12878_chr"+chrom+"_binary_p300"+".txt"
	with open(new_filename, 'a') as the_file:
		the_file.write('GM12878\tchr'+chrom+'\n')
		dat.to_csv(the_file, sep='\t', index=False)
