# Adds binary TPM feature as a new column to each chromosome file consisting of 8 histone marks
import glob
import os
import pandas as pd
import numpy as np
from random import randint, seed
from quicksect import IntervalNode

df = pd.read_table('data/condensed2.txt')
TPM = df.iloc[:,0:4]
TPM = TPM[TPM['CNhs12331']>0] # Select only the rows with TPM>0

def get_data(filename, chrom):
	dat = pd.read_table(filename,skiprows=1)
	nrows = dat.shape[0]
	dat['chr'] = chrom
	dat['lower'] = np.arange(1,200*(nrows),200)
	dat['upper'] = np.arange(200,200*(nrows+1),200).reshape(nrows,1)
	#dat.head()
	return dat

files = glob.glob('Spectacle/SAMPLEDATA_HG19_NEW/*.txt')           # create the list of files
for filename in files:
	## read file
	print filename
	chrom = os.path.split(filename)[1].split("_")[1][3:] # get chromosome number from file name
	dat = get_data(filename, chrom)
	chr_TPM = TPM[TPM['chr']==chrom]  # Select only rows with corresponding chromosome number

	## match intervals
	query = zip(dat['lower'],dat['upper'])
	data = zip(chr_TPM['lower'],chr_TPM['upper'])

	#Code taken from: https://www.biostars.org/p/99/
	def find(start, end, tree):
		"Returns a list with the overlapping intervals"
		out = []
		tree.intersect( start, end, lambda x: out.append(x) )
		return not not out
		#return True if there is an intersection

	# start the root at the first element
	start, end = data[0]
	tree = IntervalNode( start, end )

	# build an interval tree from the rest of the data
	for start, end in data[1:]:
		tree = tree.insert( start, end )

	overlap = []
	for start, end in query:
		overlap.append(find(start, end , tree))

	dat['TPM'] = map(int, overlap)
	dat = dat.ix[:, [0,1,2,3,4,5,6,7,11]]
	print dat['TPM'].value_counts()

	## write file
	new_filename = "GM12878_chr"+chrom+"_binary_new"+".txt"
	with open(new_filename, 'a') as the_file:
		the_file.write('GM12878\tchr'+chrom+'\n')
		dat.to_csv(the_file, sep='\t', index=False)
