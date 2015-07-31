# Repository Description
This repository contains all the work that I did for the DIMACS REU. See my webpage [here](http://reu.dimacs.rutgers.edu/~albertk/) for more information on the project.

# Files Description

## Data Integration
### For Spectacle Input
- `binary.py` - given multiple files (one per chromosome), adds additional feature (eRNA) to each file
- `binary_p300.py` - same as binary.py but adds additional feature of p300 overlap

### For SVM/Logit Input
- `add_feature_v1.ipynb` - given multiple files (one per chromosome), input of additional feature (p300) with output in one file for all chromosomes
- `add_feature_v2.ipynb` - given one file (containing all chromosomes), adds additional feature (tfbs) to the file

## Data Integration + Analysis 
- `exp.R` - integrate enhancer prediction from histone marks and conduct exploratory analysis of eRNA data
- `tss.ipynb` - finds distances to nearest transcription start site for each state (labeled in Spectacle)
..- `tss.R` - plots output of tss.ipynb
- `p300.ipynb` - finds overlap with p300 for each state (labeled in Spectacle)
..- `p300.R` - plots output of p300.ipynb


## Machine Learning Models
- `svm_logit.ipynb` - runs SVM and logistic regression 