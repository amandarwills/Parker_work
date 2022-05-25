# Install programs and load libraries
library(devtools)
library(magrittr)
library(ape)
library(ggplot2)
library(tidygraph)
library(fgsea)
library(data.table)
library(graphlayouts)
library(clusterProfiler)
library(pathview)
library(qusage)
library(AnnotationHub)
library(enrichplot)
library(ggnewscale)
library(GenomicRanges)

#install_github("ctlab/fgsea")
#BiocManager::install("goseq")
#BiocManager::install("GO.db")
#install.packages("ape", dependencies = T)
#install.packages("graphlayouts", dependencies = T)
#BiocManager::install("clusterProfiler",dependencies = T)
#install.packages("pathview", dependencies = T)
#install.packages("qusage", dependencies = T)
#install.packages("ggnewscale")


# Load files
deseq.data <- read.table("RNAseq_AM_labelled_INF_AM.csv", header=T, check.names = FALSE, sep=',')
deg.list <- read.table("RNAseq_AM_labelled_INF_AM_diffExpGenes.txt", header=T, check.names = FALSE, sep='\t')
pathways <- read.gmt("c2.cp.kegg.v7.4.entrez.gmt")
pathways <- Mm.c2

# Make OrgDB
hub = AnnotationHub()
head(unique(hub$species))
orgs <- subset(hub, hub$rdataclass == "OrgDb")
species_hub <- subset(hub, species == "Mus musculus")
orgdb <- query(orgs, "Mus musculus")
org.DB <- hub[["AH95960"]]
keytypes <- keytypes(org.DB)
columns(org.DB)

#Format dataframe to GSEA requirements (set columns to keep, remove all other columns, reorder columns, remove rows without matching KEGG ID)
deg.IDs <- deg.list %>% left_join(deseq.data)
filtered.deg.IDs <- deg.IDs[!duplicated(deg.IDs$Gene.ID), ]
keeps <- c("Gene.ID","Gene","log2FoldChange","stat","pvalue","padj")
filtered.deg.IDs <- filtered.deg.IDs[ , keeps, drop = FALSE]
na.strings=c(""," ","NA")
filtered.deg.IDs <- na.omit(filtered.deg.IDs, c("pvalue"))

# Convert biological IDs using OrgDb object via the bitr function to get gene names
orgDb.genename.IDs <- filtered.deg.IDs %>% pull(Gene)
genename.conversion <- bitr(orgDb.genename.IDs, fromType='ENSEMBL', toType='ENTREZID', org.DB)
gseaData <- filtered.deg.IDs %>% inner_join(genename.conversion, 
                                                    by = c('Gene' = 'ENSEMBL'))
gseaData <- gseaData[!duplicated(gseaData$Gene), ]

write.csv(gseaData, file = "RNAseq_AM_RefSeqIDs.csv")

### Start FGSEA ###

# Create ranks
ranks <- gseaData$log2FoldChange
names(ranks) <- gseaData$ENTREZID
head(ranks)

# Run FGSEA
fgseaRes <- fgsea(pathways = pathways, 
                  stats = ranks,
                  minSize  = 1,
                  maxSize  = 500)

head(fgseaRes[order(pval, -abs(NES)), ], n=40)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(pval)

# Make table plots for top pathways and remove pathways not applicable
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=40), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=40), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
topPathways
remove <- c("KEGG_PARKINSONS_DISEASE","KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS","KEGG_PATHWAYS_IN_CANCER",
            "KEGG_PROSTATE_CANCER","KEGG_DILATED_CARDIOMYOPATHY","KEGG_HYPERTROPHIC_CARDIOMYOPATHY_HCM",
            "KEGG_ARRHYTHMOGENIC_RIGHT_VENTRICULAR_CARDIOMYOPATHY_ARVC","KEGG_RENAL_CELL_CARCINOMA",
            "KEGG_EPITHELIAL_CELL_SIGNALING_IN_HELICOBACTER_PYLORI_INFECTION","KEGG_ASTHMA",
            "KEGG_VIRAL_MYOCARDITIS","KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION","KEGG_ALLOGRAFT_REJECTION",
            "KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION","KEGG_VASOPRESSIN_REGULATED_WATER_REABSORPTION",
            "KEGG_HUNTINGTONS_DISEASE","KEGG_ARRHYTHMOGENIC_RIGHT_VENTRICULAR_CARDIOMYOPATHY_ARVC",
            "KEGG_MELANOGENESIS","KEGG_GLIOMA","KEGG_MELANOMA","KEGG_INSULIN_SIGNALING_PATHWAY",
            "KEGG_CARDIAC_MUSCLE_CONTRACTION","KEGG_CHRONIC_MYELOID_LEUKEMIA","KEGG_AUTOIMMUNE_THYROID_DISEASE",
            "KEGG_TYPE_I_DIABETES_MELLITUS","KEGG_LONG_TERM_DEPRESSION","KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
            "KEGG_ALZHEIMERS_DISEASE","KEGG_ARACHIDONIC_ACID_METABOLISM","KEGG_GRAFT_VERSUS_HOST_DISEASE","KEGG_AXON_GUIDANCE",
            "KEGG_LONG_TERM_POTENTIATION","KEGG_SMALL_CELL_LUNG_CANCER","KEGG_GAP_JUNCTION","KEGG_HEMATOPOIETIC_CELL_LINEAGE",
            "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY","KEGG_TGF_BETA_SIGNALING_PATHWAY","KEGG_LEISHMANIA_INFECTION")
topPathways <- topPathways[! topPathways %in% remove]
topPathways

# Make data frame for top pathways 
topUp <- fgseaRes %>% 
  filter(ES > 0) %>% 
  top_n(20, wt=-padj)
topDown <- fgseaRes %>% 
  filter(ES < 0) %>% 
  top_n(20, wt=-padj)
topPathways.df <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
topPathways.df.Tidy <- topPathways.df %>%
  as_tibble() %>%
  arrange(desc(NES))

#write.csv(complete_topPathways, "MC_TP3_PSEAres.csv")

# Create traditional enrichment plot
plotGseaTable(pathways[topPathways.df$pathway], 
              ranks, 
              fgseaRes, 
              gseaParam = 0.5)

# Make enrichment plot for pathways of interest
plotEnrichment(pathways[["00480__Glutathione_metabolism"]],
               ranks)

# Create enrichment plot with ranked NES
ggplot(topPathways.df.Tidy, aes(reorder(topPathways, NES), NES)) +
  geom_col(aes(fill=NES<0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG pathways NES from PSEA for M. capitata TP 3")





