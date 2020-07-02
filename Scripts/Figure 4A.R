##Translation heatmap
##Resistant

library(anota2seq)
library(reshape2)
library(pheatmap)
library(readr)
library(ggplot2)
library(GO.db)
library(ggrepel)

source("polysome_profiling_functions.R")

counts <- list()
  
counts$Input <- read.delim("Resistant_input.txt")
counts$Poly <- read.delim("Resistant_poly.txt")

use_only_coding_genes <- TRUE
path_gene_classes <- "gene_classes.txt"
  
if(use_only_coding_genes == TRUE){
  gene_classes <- read.csv(path_gene_classes, header=FALSE, sep=";")
  gene_classes <- gene_classes[,1:2]
  colnames(gene_classes) <-c("Gene", "Class")
  gene_classes$Class <- trimws(gene_classes$Class) #trims white spaces
  
  coding_genes <- gene_classes[which(gene_classes$Class == "protein_coding"),]
  
  counts$Input <- counts$Input[which(rownames(counts$Input) %in% coding_genes$Gene),]
  counts$Poly <- counts$Poly[which(rownames(counts$Poly) %in% coding_genes$Gene),]
  print("Number of coding genes:")
  print(nrow(counts$Input))
}

#Make sure that the samples are paired correctly
sample_types_input <-  colnames(counts$Input)
sample_types_poly <-  colnames(counts$Poly)

sample_types_input2 <- gsub("INPUT", "", sample_types_input, ignore.case=TRUE)
sample_types_poly2 <- gsub("POLY", "", sample_types_poly, ignore.case=TRUE)

### Create an annotation file matrix
if(length(sample_types_input2) == length(sample_types_input2)){
  print("Ok. The number of input and poly samples is the same!")
}else{
  print("Error. Number of input and poly samples are not the same! You need to do something about it.")
}

#Getting input indices based on poly ordering
input_indices <- match(sample_types_poly2, sample_types_input2)
counts$Input <- counts$Input[,input_indices]

sample_types_input <-  colnames(counts$Input)
sample_types_input2 <- gsub("INPUT", "", sample_types_input, ignore.case=TRUE)

for(i in 1:length(sample_types_input2)){
  if(sample_types_input2[i] != sample_types_poly2[[i]]){
    print("Error. Samples do not match up! Need to check why input and poly samples are not paired.")
  }
}


counts$Input <- counts$Input[,!(colnames(counts$Input) %in% c("EV26.Input",
                                                              "EV4.Input",
                                                              "EV12.Input"))]
counts$Poly <- counts$Poly[,!(colnames(counts$Poly) %in% c("EV26.Poly",
                                                              "EV4.Poly",
                                                              "EV12.Poly"))]

##Normalising using voom

preProcess <- anota2seqRNAseqPreProcessing(
  dataP=counts$Poly,     # raw polysome-associated mRNA/RPF counts
  dataT=counts$Input,     # raw total mRNA counts
  transformation = "voom",  # data transformation algorithm
  filterZeroGenes = FALSE # remove genes with 0 counts
)

counts_all_voom <- data.frame(preProcess$dataT, preProcess$dataP)

upper_gene_names <- toupper(rownames(counts_all_voom))
counts_all_voom$Genes <- upper_gene_names
counts <- aggregate(. ~ Genes, counts_all_voom, mean)
rownames(counts) <- counts$Genes
counts$Genes <- NULL

genes_to_appear <- as.character(read.delim("Translation_genes_for_heatmap.txt",
                                           header=FALSE)[,1])
local_counts <- counts[rownames(counts) %in% genes_to_appear,]
local_counts$Gene <- rownames(local_counts)

local_counts_melt <- melt(local_counts)
local_counts_melt$Sample <- substr(local_counts_melt$variable, 1, 2)
local_counts_melt$Source <- NA
local_counts_melt$Source[grep("Input", local_counts_melt$variable, ignore.case=TRUE)] <- "Input"
local_counts_melt$Source[grep("Poly", local_counts_melt$variable, ignore.case=TRUE)] <- "Poly"
local_counts_melt$Sample_source <- paste(local_counts_melt$Sample, local_counts_melt$Source, sep="_")

local_counts_mean <- aggregate(value ~ Gene+Sample_source, local_counts_melt, median)

local_counts_mean_cast <- dcast(local_counts_mean, Gene~Sample_source)
rownames(local_counts_mean_cast) <- local_counts_mean_cast$Gene
local_counts_mean_cast$Gene <- NULL

local_counts_mean_cast <- local_counts_mean_cast[,c("CX_Poly","Ct_Poly","CM_Poly")]

sig <- "translation"
experiment <- "resistant"

#italics
newnames <- lapply(
  rownames(local_counts_mean_cast),
  function(x) bquote(italic(.(x))))

pheatmap(local_counts_mean_cast, main=paste(toupper(experiment), "-", sig),
         fontsize=4, scale="row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         border_color = "black",
         labels_row = as.expression(newnames))
