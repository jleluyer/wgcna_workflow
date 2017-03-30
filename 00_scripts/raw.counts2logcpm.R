#!/usr/bin/Rscript

## test anova GxE interaction ###
#see post
ls()
rm(list=ls())
ls()
#https://support.bioconductor.org/p/56568/
## test Glm approahc with edger ##
library('edgeR')



counts <-read.table("02_data/count_matrix.txt",header=T,na.strings="NA")


# Create design
targets <-read.table("01_info_file/design.txt",header=T,na.strings="NA")
targets$genotype <- relevel(targets$genotype, ref=c("wild"))

targets


#reassignng column name
rownames(counts) <- counts[,1]
counts[,1]<- NULL
dim(counts)
str(counts)
colSums(counts) / 1e06

#nb gene with low count
table( rowSums( counts ) )[ 1:30 ]
#create object
Group <-  factor(paste(targets$environment,targets$genotype,sep="."))
Group
cds <- DGEList( counts , group = Group )
head(cds$counts)
cds$samples
levels(cds$samples$group)

#filtering low counts
cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 2, ]
dim( cds )

#calcultate normalization factors
cds <- calcNormFactors( cds )
cds$samples

#library size
cds$samples$lib.size*cds$samples$norm.factors

#Create design (can be created manually)
design <- model.matrix(~genotype *environment, data=targets)
colnames(design)
design
#library size
cds$samples$lib.size*cds$samples$norm.factors
cds$samples$norm.factors
#estimate dispersion
cds <-estimateDisp(cds,design) #estimate all dispersions in one run #add different prior values (prior.df = INT)
names( cds )

#create a logcpm matrix
logcpm <- cpm(cds, prior.count=2, log=TRUE)
write.table(logcpm,file="02_data/logcpm_edger.csv",sep=";",row.names=T,quote=F)
