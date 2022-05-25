# Install programs and load libraries
library(devtools)

#BiocManager::install("goseq")
#BiocManager::install("GO.db")
#install.packages("ape", dependencies = T)
#install.packages("graphlayouts", dependencies = T)
#BiocManager::install("clusterProfiler",dependencies = T)
#install.packages("pathview", dependencies = T)
#install.packages("qusage", dependencies = T)
#install.packages("ggnewscale")

library(dplyr)
library(ape)
library(ggplot2)
library(tidygraph)
library(data.table)
library(goseq)
library(GO.db)
library(graphlayouts)
library(clusterProfiler)
library(pathview)
library(qusage)
library(magrittr)
library(AnnotationHub)
library(AnnotationDbi)
library(enrichplot)
library(ggnewscale)
library(GenomicRanges)

# Load files
deseq.data <- read.table("RNAseq_PMN_labelled_INF.csv", header=T, check.names = FALSE, sep=',')
deg.list <- read.table("RNAseq_PMN_labelled_INF_diffExpGenes.csv", header=T, check.names = FALSE, sep=',')

genename_output <- "RNAseq_PMN_labelled_INF_genenames.csv"
clusterprofiler_output <- "RNAseq_PMN_labelled_INF_clusterprofiler_results.csv"
enrichGO_output <- "RNAseq_PMN_labelled_INF_enrichGO_results.csv"

### START CLUSTERPROFILER
search_kegg_organism('Mus musculus', by = "scientific_name")

# Make OrgDB
hub = AnnotationHub()
head(unique(hub$species))
orgs <- subset(hub, hub$rdataclass == "OrgDb")
species_hub <- subset(hub, species == "Mus musculus")
orgdb <- query(orgs, "Mus musculus")
org.DB <- hub[["AH95960"]]
keytypes <- keytypes(org.DB)
columns(org.DB)

# Convert biological IDs using OrgDb object via the bitr function to get gene names
orgDb.genename.IDs <- deseq.data %>% pull(Gene)
genename.conversion <- bitr(orgDb.genename.IDs, fromType='ENSEMBL', toType='GENENAME', org.DB)
filtered.deg.genenames <- deseq.data %>% left_join(genename.conversion, 
                                                    by = c('Gene' = 'ENSEMBL'))
write.csv(filtered.deg.genenames, file = genename_output, row.names = FALSE)

# Convert biological IDs using OrgDb object via the bitr function to get entrez IDs
orgDb.ensemble.IDs <- deseq.data %>% pull(Gene)
entrez.conversion <- bitr(orgDb.ensemble.IDs, fromType='ENSEMBL', toType='ENTREZID', org.DB)
deg.entrezIDs <- deseq.data %>% inner_join(entrez.conversion, 
                                                       by = c('Gene' = 'ENSEMBL'))

# Filter data for DEGs and match organism gene IDs for clusterProfiler
deg.IDs <- deg.list %>% left_join(deg.entrezIDs)
filtered.degs <- deg.IDs[!duplicated(deg.IDs$Gene), ]
na.strings=c(""," ","NA")
filtered.degs <- na.omit(filtered.degs, c("pvalue"))


# Use BioMart to convert biological IDs (if bitr doesn't return 1:1 mapping between keys)
#require(biomaRt)
#mart <- useMart("ENSEMBL_MART_ENSEMBL")
#ensembl <- useEnsembl(biomart = "ensembl")
#searchDatasets(mart = ensembl, pattern = "(M|m)ouse")
#mart <- useDataset("mmusculus_gene_ensembl", mart)
#bm.geneID.conversion <- getBM(mart=mart, attributes=c("ensembl_gene_id","entrezgene_id"),
#      filter="ensembl_gene_id", values=orgDb.genename.IDs, uniqueRows=TRUE
#filtered.deg.IDs <- filtered.deg.IDs %>% inner_join(bm.geneID.conversion, by = c('Gene' = 'ensembl_gene_id'))


# Set universe background data
universe <- deg.entrezIDs %>% pull(ENTREZID)

# Set significant genes and extract IDs
siGenes <- deg.entrezIDs %>% 
  filter(pvalue < 0.05, !is.na(ENTREZID)) %>% pull(ENTREZID)

# Run enrichKEGG
enrichKEGG <- enrichKEGG(gene = siGenes,
                         organism = "mmu",
                         keyType = "ncbi-geneid",
                         pvalueCutoff = 0.05)

enrichKEGG_tidy <- enrichKEGG %>% 
  slot("result") %>% 
  tibble::as_tibble() 

head(enrichKEGG_tidy)
write.csv(enrichKEGG_tidy, file = clusterprofiler_output, row.names = FALSE)

# Run enrichGO
enrich_go <- enrichGO(
  gene= siGenes,
  OrgDb = org.DB,
  keyType = "ENTREZID",
  ont = "BP",
  qvalueCutoff = 0.05,
  readable=TRUE
)

enrich_go_tidy <- enrich_go %>% 
  slot("result") %>% 
  tibble::as_tibble() 

head(enrich_go_tidy)
write.csv(enrich_go_tidy, file = enrichGO_output, row.names = FALSE)

dotplot(enrich_go,
        color = "p.adjust",
        showCategory = 25,
        size = NULL,
        split = NULL,
        font.size = 8,
        title = "RNAseq_PMN_labelled_Ifected enrichment results",)


# View KEGG maps - replace second argument with KEGG pathway ID (e.g. mmu00100) and the KEGG map will open in a browser
browseKEGG(enrichKEGG, 'mmu04974')

# End