#
# Rawcopy analysis of CNV data
#



# 1. Install necessary packages ----------------------------------------------------------

# bioconducter packages
source("http://bioconductor.org/biocLite.R")
biocLite("affxparser")
biocLite("DNAcopy")
biocLite("aroma.light")

# cran packages
install.packages(c('foreach','doMC','PSCBS','squash','digest','ape','SDMTools'))

# rawcopy package
install.packages('rawcopy',repos='http://array.medsci.uu.se/R',type='source')

# 2. rawcopy analysis --------------------------------------------------------------------

library(rawcopy)

# Run samples in a directory on multiple cores (~6GB of RAM is required per core):
 setwd("D:/SNP 6/Rohdaten/cel")
rawcopy()



