# Loading library
library(ggfortify)
library(DESeq2)
library(ggplot2)
library(scales)
library(devtools)
library(ggpubr)
library(magrittr)
library(stats)

output.file <-  "AEC_uninfected_network_info.csv"
  
# Load files
data <- read.table("RNAseq_all_cell_types_labelled_matrix.csv", header=T, check.names=FALSE, row.names = 1, sep=',')
sample.info <- read.table("Sample_info.csv", header=T, sep=',', row.names = 1)

# Check data for correct data mode and no NA values
sapply(data, mode)
sum(is.na(data))

# Transpose the data to have variables (genes) as columns
t.data <- t(data)
dim(data_for_PCA)

# Annotate the data with condition group as labels
data_for_PCA <- merge(sample.info, t.data, by = 0)
rownames(data_for_PCA)=data_for_PCA$Row.names
data_for_PCA$Row.names <- NULL

# Create PCA object
data.pca <- prcomp(data_for_PCA[,c(6:35855)],
                   center = TRUE,
                   scale. = TRUE)

# Summary of the prcomp object
summary(data.pca)

# Structure of the PCA object
str(data.pca)

# PCA plot
ggbiplot(data.pca,
         data = data_for_PCA,
         obs.scale = 1, 
         var.scale=1,
         ellipse=T,
         ellipse.prob = 0.68, 
         labels = NULL,
         var.axes = FALSE,
         circle=F,
         varname.size=3,
         var.axes=T,
         groups = data_for_PCA$Treatment) 


write.csv(filtered.deg.genenames, output.file, row.names = FALSE)

