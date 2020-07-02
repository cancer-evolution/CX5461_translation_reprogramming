##Translation heatmap
##Resistant

library(anota2seq)
library(reshape2)
library(ggplot2)

source("~/Documents/CX5461/Paper_2020_scripts/polysome_profiling_functions.R")

# files <- c("~/Documents/CX5461/0318_design_matrix_HTSeq_Expression_Matrix.txt",
#            "~/Documents/CX5461/0319_design_matrix_HTSeq_Expression_Matrix.txt")
# 
# sample_types <- c("CMB", "CX", "EV", "Ctrl")
# contrasts <- c("CMB_Ctrl", "CMB_CX", "CX_Ctrl")
# 
# experiment <- "resistant"
# 
# counts <- read_counts2(files)
# 
# 
# new_files <- c("~/Documents/CX5461/Counts/polysome_design_matrix_nov_HTSeq_Expression_Matrix.txt")
# txtList = lapply(new_files, function(x){ read.delim(x)} )
# 
# txtList = lapply(txtList, function(x){
#   colnames(x) <- as.character(unname(unlist(x[1,])))
#   x <- x[-1,]
#   rownames(x) <- x[,1]
#   x <- x[,-1]
#   return(x)
# })
# txtList <- txtList[[1]]
# 
# all_genes <- rownames(txtList)
# all_genes <- all_genes[!(all_genes %in% c("__ambiguous", "__no_feature"))]
# 
# txtList <- txtList[rownames(txtList) %in% all_genes,]
# 
# matrix_input <- txtList[,grep("Input", colnames(txtList), ignore.case=TRUE)]
# 
# matrix_poly <- txtList[,grep("Poly", colnames(txtList), ignore.case=TRUE)]
# 
# matrix_input2 <- apply(matrix_input, 2, function(x){ as.numeric(as.character(x))})
# matrix_poly2 <- apply(matrix_poly, 2, function(x){ as.numeric(as.character(x))})
# 
# rownames(matrix_input2) <- all_genes
# rownames(matrix_poly2) <- all_genes
# 
# colnames(matrix_input2) <- paste(colnames(matrix_input2), "NEW", sep="_")
# colnames(matrix_poly2) <- paste(colnames(matrix_poly2), "NEW", sep="_")
# 
# print("Input matrix column names")
# print(colnames(matrix_input2))
# 
# print("Poly matrix column names")
# print(colnames(matrix_poly2))
# 
# 
# counts2 <- list(Input = matrix_input2,
#                 Poly = matrix_poly2)
# 
# 
# all_genes_both <- sort(unique(c(all_genes, rownames(counts$Input))))
# 
# counts$Input <- counts$Input[match(all_genes_both, rownames(counts$Input)),]
# counts2$Input <- counts2$Input[match(all_genes_both, rownames(counts2$Input)),]
# 
# counts$Poly <- counts$Poly[match(all_genes_both, rownames(counts$Poly)),]
# counts2$Poly <- counts2$Poly[match(all_genes_both, rownames(counts2$Poly)),]
# 
# counts$Input <- cbind(counts$Input, counts2$Input)
# counts$Poly <- cbind(counts$Poly, counts2$Poly)
# 
# counts$Input[is.na(counts$Input)] <- 0
# counts$Poly[is.na(counts$Poly)] <- 0
# 
# #Select samples to use:
# counts$Input <-  counts$Input[,c("CMB8-Input","CMB9-Input","CMB31-Input",
#                                  "CX2-Input","CX14-Input_NEW","CX39-Input_NEW",
#                                  "Ctrl33-Input","Ctrl38-Input","Ctrl6-Input")]
# 
# counts$Poly <- counts$Poly[,c("CMB8-Poly", "CMB9-Poly", "CMB31-Poly",
#                               "CX2-Poly", "CX14-Poly_NEW", "CX39-Poly_NEW",
#                               "Ctrl33-Poly", "Ctrl38-Poly", "Ctrl6-Poly")]

counts <- list()
counts$Input <- read.delim("~/Documents/CX5461/Paper_2020_counts/Resistant_input.txt")
counts$Poly <- read.delim("~/Documents/CX5461/Paper_2020_counts/Resistant_poly.txt")
  
  
use_only_coding_genes <- TRUE
path_gene_classes <- "~/Documents/CX5461/Paper_2020_text_files/gene_classes.txt"
  
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


preProcess <- anota2seqRNAseqPreProcessing(
  dataP=counts$Poly,     # raw polysome-associated mRNA/RPF counts
  dataT=counts$Input,     # raw total mRNA counts
  transformation = "voom",  # data transformation algorithm
  filterZeroGenes = FALSE # remove genes with 0 counts
)

counts_all_voom <- data.frame(preProcess$dataT, preProcess$dataP)
#colnames(counts_all_voom) <- gsub(".Input", "", colnames(counts_all_voom))
#colnames(counts_all_voom) <- gsub(".Poly", "", colnames(counts_all_voom))

upper_gene_names <- toupper(rownames(counts_all_voom))
counts_all_voom$Genes <- upper_gene_names
counts <- aggregate(. ~ Genes, counts_all_voom, mean)
rownames(counts) <- counts$Genes
  
local_counts_melt <- melt(counts)
local_counts_melt$Sample <- substr(local_counts_melt$variable, 1, 2)
local_counts_melt$Source <- NA
local_counts_melt$Source[grep("Input", local_counts_melt$variable, ignore.case=TRUE)] <- "Input"
local_counts_melt$Source[grep("Poly", local_counts_melt$variable, ignore.case=TRUE)] <- "Poly"
local_counts_melt$Sample_source <- paste(local_counts_melt$Sample, local_counts_melt$Source, sep="_")

local_counts_melt <- local_counts_melt[local_counts_melt$Sample != "EV",]
local_counts_melt$Sample <- factor(local_counts_melt$Sample, levels=c("Ct", "CX", "CM"))

pdf("~/Documents/CX5461/Paper_2020_3_resistant_controls/Expression of Rap1a - resistant.pdf",
    height=4, width=4)
g <- ggplot(subset(local_counts_melt, Source=="Poly" & Genes == "RAP1A"), aes(x=Sample, y =value))+
  geom_boxplot()+
  ylab("Expression")+
  theme_bw()
print(g)
dev.off()

