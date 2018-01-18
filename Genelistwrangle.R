################################################################################
# Data Wrangling of excel tables obatined from Deseq2 and EdgeR tables
#
# Can do ...
#   - read tables and csv files
#   - select rquired columns and rows
#   - plot ggplot2 volcano plot
#   - subset data
#
################################################################################

###############################################################################################################################

#---------------
# LOAD PACKAGES
#---------------
library(ggplot2)
library(gridExtra)
library(tidyverse)
library("AnnotationDbi")
library("org.Mm.eg.db")
library(tidyr)
library(data.table)
library(ggrepel)
###############################################################################################################################

mtc <- read_tsv(file = "/Users/bachum/Desktop/OneDrive/Pan-Richard/tables/KOP2vsWtP2.complete.txt", col_names = T)
colnames(mtc)
glimpse(mtc)
dim(mtc)
#remove all NA rows
mtc <- mtc %>% drop_na()
glimpse(mtc)
dim(mtc)
colnames(mtc)
#mtc <- read.csv(file="Scatter_Plot_NT.csv", header = T)
tail(mtc)
#Select only the required columns
mtc <- dplyr::select(mtc,Id, WtP2, KOP2, KOP3, FC, log2FoldChange,pvalue,padj)
head(mtc)
dim(mtc)

###############################################################################################################################

#Remove dot from ensemble annotation using command and add an additional variable Gene_ID
mtc <- mutate(mtc, Gene_Id = gsub("\\..*","",mtc$Id))
mtc <- dplyr::select(mtc,Gene_Id, WtP2, KOP2, KOP3, FC, log2FoldChange,pvalue,padj)
dim(mtc)
head(mtc)
###############################################################################################################################
#Add Annotation and add gene symbols and entrez IDs

columns(org.Mm.eg.db)
mtc$symbol <- mapIds(org.Mm.eg.db,
                     keys=mtc$Gene_Id,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
mtc$entrez <- mapIds(org.Mm.eg.db,
                     keys=mtc$Gene_Id,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

mtc$GeneName <- mapIds(org.Mm.eg.db,
                     keys=mtc$Gene_Id,
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")
colnames(mtc)
head(mtc)
mtc <- dplyr::select(mtc,symbol, WtP2, KOP2, KOP3, FC, log2FoldChange,pvalue,padj,entrez,Gene_Id,GeneName)
write.csv(mtc, file = "/Users/bachum/Desktop/OneDrive/Pan-Richard/tables/KOP2vsWtP2.complete.csv", row.names = T)

damup <- read_csv("/Users/bachum/Desktop/OneDrive/Pan-Richard/tables/Dam-Genes.csv")
glimpse(damup)
head(damup)
class(damup)

###############################################################################################################################
#How to subset a column in data frame based on another data frame/list
#A better option would be data.table
library(data.table)
setDT(mtc)[symbol %chin% damup$ID]
#or
subset(mtc, symbol %in% damup$ID)
#or
mtc %>%
  filter(symbol %in% damup$ID)
damup_list <- mtc %>%
  filter(symbol %in% damup$ID)
head(damup_list)
head(mtc)
dim(damup_list)
write.csv(damup_list, file = "/Users/bachum/Desktop/OneDrive/Pan-Richard/tables/KOP2vsWtP2.damGene_list.csv", row.names = T)


###############################################################################################################################
#check row for row if a combination exists in another dataframe and add annotation saying yes if found and no if not
mtc$Dam.Gene <- ifelse(is.na(match(paste0(mtc$symbol), 
                                paste0(damup$ID))),"No", "Yes")
filter(mtc, Dam.Gene == "Yes")
colnames(mtc)
p <- ggplot(mtc, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=Dam.Gene), size =1, alpha = 1) + scale_color_manual(values=c("grey", "red")) + 
  theme_bw(base_size = 10) 
p
#plot over points
q <- ggplot(mtc, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=Dam.Gene), size =1, alpha = 1, color = 'grey') + 
geom_point(data = subset(mtc, Dam.Gene == 'Yes'), color = 'red') + theme_bw(base_size = 10)
q

library(ggrepel)
#Select your gene of interest using the row number of excel sheet
r <- q+geom_text_repel(data=dplyr::filter(mtc, symbol %in% c("Ctsz",
                                                             "Ctss",
                                                             "Apoe",
                                                             "Apoc1",
                                                             "Apoc4",
                                                             "Npc2",
                                                             "Ch25h",
                                                             "Lpl",
                                                             "Trem2",
                                                             "Axl",
                                                             "Tyrobp",
                                                             "Igf1",
                                                             "Spp1",
                                                             "Gpnmb",
                                                             "Itgax",
                                                             "Csf1",
                                                             "Clec7a",
                                                             "Lyz2",
                                                             "Lgals3bp",
                                                             "Hifa",
                                                             "Plin2",
                                                             "Dpp7",
                                                             "Hexa",
                                                             "Ank",
                                                             "Acaca",
                                                             "Soat1")), aes(label= symbol))
r + geom_vline(xintercept = c(-1,1)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

###############################################################################################################################
# FACETING 
# To do ...
#   - multiple small multiple charts

library(ggplot2)
# Basic barplot
colnames(damup_list)
forbar <- dplyr::select(damup_list, symbol, WtP2, KOP2, KOP3)
glimpse(forbar)

#Creating a long table from a wide table
library(tidyr)
forbar.l <- gather(forbar, key = "Genotype", value = "CPM" , WtP2,KOP2,KOP3)
head(forbar.l)
tail(forbar.l)
colnames(forbar.l)
library(ggplot2)

plot <- ggplot(forbar.l, aes(x = symbol, y = CPM, fill = Genotype)) + geom_bar(stat = "identity", position=position_dodge())
plot + facet_wrap(~ symbol, nrow = 20, ncol = 20)
