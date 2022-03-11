# This script applies Upper-Quartile normalization from LPE library

# Define input and output files and read the data
myfile = "/mnt/64716603-5b56-4f9a-b195-c11560647a3a/data/IMMOTION151/immotion151.filtered.TPMs.csv"
outfile = "/mnt/64716603-5b56-4f9a-b195-c11560647a3a/data/IMMOTION151/immotion151.filtered.uqTPMs.csv"
mydata = read.csv(myfile, sep='\t')


# Install and load LPE library
#https://davetang.org/muse/2011/01/24/normalisation-methods-for-dge-data/
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("LPE")
library(LPE)

# Apply UQ normalization and save to output file
uqmydata = quartile.normalize(mydata[,2:dim(mydata)[2]], percent=75)
uqmydata <- cbind(mydata[,1], uqmydata)
colnames(uqmydata)[1] <- "symbol"
write.csv(uqmydata, file = outfile, row.names = FALSE)