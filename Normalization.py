# This script describes the steps to pre-processed RNAseqmatrix counts derived from the IMmotion151 trial (log2TPM)
import pandas as pd
import numpy as np
from scipy.stats import iqr

############ Load data ############
c = pd.read_excel('mmc2(4).xlsx')
df = pd.read_csv('IMmotion151.expression.data.TPM.csv',sep=',') # Controlled access data. Processed RNAseq data (EGAD00001006618)
df.set_index('symbol',inplace=True)
df.drop(columns=['Unnamed: 0','gene_name'],inplace = True)

############ Data pre-processing ################
# Data preprocessing steps follow the steps described in Braun et al paper (https://doi.org/10.1038/s41591-020-0839-y)
interq = df.apply(lambda x: iqr(x), axis=0)
filtlowIQR = list(interq[interq >= 0.5].index) # Remove samples with IQR<0.5 (taking log2(TMP) as inputs)
df = df[filtlowIQR].copy()
df = df.apply(lambda x: 2**(x)-1) # Converto back to TPM
lowexpression = df.apply(lambda x: np.count_nonzero(x!=0),  axis = 0) # Filter samples with less than 15000 expressed genes (looking TPMs)
filtlowexpression = list(lowexpression[lowexpression >= 15000].index)
filteredsamples = list(lowexpression[lowexpression < 15000].index)
df = df[filtlowexpression].copy() # in principle, 4 samples are down, passing 819 samples this filter
df = df.loc[(df != 0).any(axis=1)] # Remove rows genes with zero expression in all genes (looking TPMs)
df.reset_index(inplace=True)
df.drop_duplicates(subset='symbol', keep='first', inplace = True) # remove duplicated genes
df.set_index('symbol',inplace = True)
df.to_csv('immotion151.filtered.TPMs.csv', sep='\t')
############# Apply UQ normalization ############
# Outpuf file ('immotion151.filtered.TPMs.csv') is loaded into R by UQ-normalization.R script to apply Upper-Quartile Normalization as described in Braun et al
# Ref: https://doi.org/10.1038/s41591-020-0839-y
# UQ-normalization.R will generate the output file 'immotion151.filtered.uqTPMs.csv', used for further analysis (Figure 1 and Supp Tables)
