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


# read the data
affydata <- ReadAffy()

ed <- exprs(affydata)
samp <- sampleNames(affydata)
probes <- featureNames(affydata)


nvals <- rma(affydata)
ned <- exprs(nvals)

nsamp <- sampleNames(nvals)
nprobes <- featureNames(nvals)

plot(ned[,"GSM375719.CEL.gz"], ned[,"GSM375721.CEL.gz"],pch=".",
     +      xlab="GSM277408", ylab="GSM277414",
     +      main="Logjam mutant vs Control (Abdomen)")
abline(0,1,col="blue")
abline(log(2),1,col="red")
abline(log(.5),1,col="red")

identify(ned[,"GSM375719.CEL.gz"], ned[,"GSM375721.CEL.gz"], nprobes)


# 4 Get some annotation
# get the package that contains the annotation for this array

# load it and do some R magic!
x <- drosophila2GENENAME
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

adf[1:10,]

biocLite("hgu133plus2.db")
library("hgu133plus2.db")
x <- hgu133plus2ACCNUM
# Get the probe identifiers that are mapped to an ACCNUM
mapped_probes <- mappedkeys(x)
# Convert to a list

xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
        # Get the ACCNUM for the first five probes
        xx[1:5]
        # Get the first one
        xx[[1]]
}

vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

adf[1:10,]

# 5 Process and Export the data
# ned is our normalized expression matrix

# merge with our annotation
aned <- merge(ned, adf, by="row.names", all.x=TRUE, sort=FALSE)
write.csv(aned, "Her.csv", row.names=FALSE)


x <- hgu133plus2ENTREZID
# Get the probe identifiers that are mapped to an ENTREZ Gene ID
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
        # Get the ENTREZID for the first five probes
        xx[1:5]
        # Get the first one
        xx[[1]]
}
vals <- sapply(xx, as.vector)
adf2 <- data.frame(probe=names(vals), gene=vals)

adf2[1:10,]
aned2 <- merge(aned, adf2, by="row.names", all.x=TRUE, sort=FALSE)
write.csv(aned2, "Herceptin.csv", row.names=FALSE)


vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

adf[1:10,]
aned <- merge(ned, adf, by="row.names", all.x=TRUE, sort=FALSE)
write.csv(aned, "Her2.csv", row.names=FALSE)