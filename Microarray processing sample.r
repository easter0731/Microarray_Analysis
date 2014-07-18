Skip to main content
ARK-Genomics Logo
Roslin Logo
Search form Search   
LoginRegister
Home
About Us
Services
Protocols
Outputs
Events
News
Contact Us
Working with Affymetrix CEL files

Microarrays are still used in research, though not as much as they use to be due to the rise of NGS.  Still, I get requests for bioinformatics support around analysis of microarray data, and the most common platform is Affymetrix. 

We need to get some data: download ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/series/GSE10940/GSE10940_RAW.tar which refers to a study in drosophila involving a mutant where the logjam gene has been knocked out.  Extract the contents to the default directory "GSE10940_RAW". 

In R, choose File -> Change dir and change the current working directory to the "GSE10940_RAW" directory.

Here is some R code that should help

# get the necessary libraries.  The first six lines only need to be executed the first time you run this code
source("http://www.bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("affyPLM")
biocLite("drosophila2cdf")
biocLite("org.Dm.eg.db")
biocLite("drosophila2.db")
library(affy)
library(affyPLM)
library(drosophila2cdf)
library(org.Dm.eg.db)
library(drosophila2.db)

# Data:
# We are using public data from GEO accession: GSE10940
#
# Details:
# We compared expression profiles from loj mutant females to those of controls.
# We compared mutant and control abdomen tissue as well as mutant and control
# tissue preparations from the remainder of the fly (the head/thorax). Three
# mutant and three control RNA samples for each set of tissue were used for array
# hybridization (12 total array samples).

# The samples look like this:
# ID          Name         
#  1 GSM277408     Control Abd_1
#  2 GSM277409     Control Abd_2
#  3 GSM277410     Control Abd_3
#  4 GSM277411     Control HT_1
#  5 GSM277412     Control HT_2
#  6 GSM277413     Control HT_3
#  7 GSM277414     Logjam Abd_1
#  8 GSM277415     Logjam Abd_2
#  9 GSM277416     Logjam Abd_3
# 10 GSM277417     Logjam HT_1
# 11 GSM277418     Logjam HT_2
# 12 GSM277419     Logjam HT_3

# 1 Read in probe level data

# The affy package will automatically download the appropriate array annotation
# when you require it. However, if you wish you may download and install the
# cdf environment you need from http://www.bioconductor.org/packages/release/data/annotation/
# manually. If there is no cdf environment currently built for your particular chip and you
# have access to the CDF file then you may use the makecdfenv package to create one
# yourself. To make the cdf packaes, Microsoft Windows users will need to use the tools
# described here: http://cran.r-project.org/bin/windows/rw-FAQ.html

# read the data
affydata <- ReadAffy()
affydata

# raw expression data
ed <- exprs(affydata)

samp <- sampleNames(affydata)
probes <- featureNames(affydata)

ed[1:10,]
probes[1:10]
samp

# 2 Normalizing Data   
#
# The Affy package has implementations of a number of normalization methods
# for single-channel arrays. This includes (among others):
#   - mas5() - Affymetrix's own normalization program
#   - rma() - 'Robust Multi-Chip' average
#   - gcrma() - A bias-corrected RMA
# GCRMA is good but takes significantly longer than RMA, so RMA is the
# most commonly used
nvals <- rma(affydata)

# normalised expression data
ned <- exprs(nvals)

nsamp <- sampleNames(nvals)
nprobes <- featureNames(nvals)

ned[1:10,]
nprobes[1:10]
nsamp

# the normalized data is on the log-scale

# 3 visualize the data
# We can plot a case vs a control
#    1 GSM277408     Control Abd_1
#    7 GSM277414     Logjam Abd_1

# we can figure out which columns these are from colnames
colnames(ned)
plot(ned[,"GSM277408.CEL.gz"], ned[,"GSM277414.CEL.gz"])

# tidy it up a bit
plot(ned[,"GSM277408.CEL.gz"], ned[,"GSM277414.CEL.gz"], pch=".",
     xlab="GSM277408", ylab="GSM277414",
     main="Logjam mutant vs Control (Abdomen)")

# add a line of y=x
abline(0,1,col="blue")

# plot lines at two fold up and down regulation
# hint: two fold upregulation is 2, so we draw the intercept at log(2)
# hint: two fold downregulation is .5, so we draw the intercept at log(.5)
abline(log(2),1,col="red")
abline(log(.5),1,col="red")

# use the identify function to label points
identify(ned[,"GSM277408.CEL.gz"], ned[,"GSM277414.CEL.gz"], nprobes)

# press escape to exit on windows, or right click to exit on Linux

# 4 Get some annotation
# get the package that contains the annotation for this array

# load it and do some R magic!
x <- drosophila2GENENAME
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

adf[1:10,]

# 5 Process and Export the data
# ned is our normalized expression matrix

# merge with our annotation
aned <- merge(ned, adf, by="row.names", all.x=TRUE, sort=FALSE)
write.csv(aned, "drosophila_normalized_expression_data.csv", row.names=FALSE)

# save data
save.image()

You should get a nice scatterplot that looks a little like this:
        
        
        
        And there we have it!  The .csv file should open in Excel, and you may now explore your data :)

Our Sponsors

Roslin Foundation RDSVS Roslin Foundation University of Edinburgh BBSRC Technology Strategy Board Europe & Scotland
Privacy and cookies policy