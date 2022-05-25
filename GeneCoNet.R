#!/usr/bin/env Rscript

##
## DESCRIPTION and information to run on linux/R
##

#Construct gene co-expression network

## Example to run in linux: Rscript GeneCoNet.R TPM_matrix diffExprGenes condition
## If running directly in R/Rstudio, run code from lines 27-93 only + change arg1, arg2, and arg3 (lines 35-37) to file names, as defined on lines 12-14

# Arg1 = TPM_matrix (tab delimited; TPM values for each gene across samples of interest) - RNAseq_AM_labelled_PBS_TPM_matrix.txt
# Arg2 = diffExprGenes (list of differentially expressed genes) - RNAseq_AM_labelled_PBS_diffExpGenes
# Arg3 = condition (Condition prefix from TPM_matrix to extract (e.g. PBS or Infected)

# There will be 2 output files : condition.corData.txt = network edges & condition.modules.txt = network modules

#args = commandArgs(trailingOnly=TRUE)

# Test if there is three arguments: if not, return an error
#if (length(args)!=3) {
#  cat(DESCRIPTION)
#  stop("Three arguments must be supplied!", call.=FALSE)
#}


## Start network construction code

#Load libraries
suppressMessages(suppressWarnings(library(DGCA, quietly = TRUE)))
suppressMessages(suppressWarnings(library(WGCNA, quietly = TRUE)))
suppressMessages(suppressWarnings(library(matrixStats, quietly = TRUE)))

# Set file variables
expression.file <- "RNAseq_PMN_labelled_PBS.csv"
diffExprGenes <- "RNAseq_PMN_labelled_PBS_diffExpGenes.csv"
condition <- "Uninfected"
output <- "PMN_labelled_PBS"

# Parameters
maxPval=0.05

# Read data
tpm.data=as.matrix(read.table(expression.file,h=T,row.names=1,sep=","))
deg.data=read.table(diffExprGenes,h=T,row.names=1,sep=",")
deg.ids=rownames(deg.data)

# Remove non numerical data and reformat the matrix
indx <- grepl(condition, colnames(tpm.data))
TP=which(indx==T)
tpm.data=tpm.data[deg.ids,TP]
C=colnames(tpm.data)
R=rownames(tpm.data)
dims=dim(tpm.data)
tpm.data=as.numeric(tpm.data)
dim(tpm.data)=dims
colnames(tpm.data)=C
rownames(tpm.data)=R

### Start computing co-expression networks

t.data=t(tpm.data)
cor.data=matCorr(t.data,corrType="pearson")
pairs.cor.data=data.frame(n1=rownames(cor.data)[row(cor.data)],n2=colnames(cor.data)[col(cor.data)],cor=c(cor.data))
cor.data.row=rownames(cor.data)
cor.data.col=colnames(cor.data)
nsample.data=matNSamp(t.data)
cor.data.pval=matCorSig(cor.data,nsample.data)
colnames(cor.data.pval)=cor.data.col
rownames(cor.data.pval)=cor.data.row
pairsPval=data.frame(n1=rownames(cor.data.pval)[row(cor.data.pval)],n2=colnames(cor.data.pval)[col(cor.data.pval)],pval=c(cor.data.pval))
cor.data.pval.vec=as.vector(cor.data.pval)
cor.data.adjPval=adjustPVals(cor.data.pval.vec,adjust="BH")
cor.data.adjPval=as.numeric(format.pval(cor.data.adjPval,digits=2,nsmall=3))
dim(cor.data.adjPval)=dim(cor.data.pval)
colnames(cor.data.adjPval)=cor.data.col
rownames(cor.data.adjPval)=cor.data.row
pairsAdjPval=data.frame(n1=rownames(cor.data.adjPval)[row(cor.data.adjPval)],n2=colnames(cor.data.adjPval)[col(cor.data.adjPval)],adjPval=c(cor.data.adjPval))
#
cor.data.val=cbind(pairs.cor.data,pval=pairsPval$pval,adjPval=pairsAdjPval$adjPval)
cor.data.val.final=cor.data.val[complete.cases(cor.data.val),]
cor.data.val.final.filtered=cor.data.val.final[cor.data.val.final$adjPval <= maxPval,]
#
corData.file=paste(output,".corData.txt",sep = "")
write.table(cor.data.val.final.filtered,file=corData.file,quote=FALSE,sep="\t",row.names=FALSE)
#
gene_tree<-hclust(as.dist(1-cor.data),method="average")
module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=10,deepSplit=TRUE)
m <- as.data.frame(module_labels)
rownames(m)=rownames(cor.data)
modules.file=paste(output,".modules.txt",sep = "")
write.table(m,file=modules.file,quote=FALSE,sep="\t",row.names=T,col.names=F)

#END
