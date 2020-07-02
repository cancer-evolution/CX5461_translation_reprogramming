library(DESeq2)
library(limma)
library(edgeR)
library(reshape2)
source("polysome_profiling_functions.R")

htseq_matrix_input <- read.delim("Resistant_input.txt")
htseq_matrix_poly <- read.delim("Resistant_poly.txt")

contrasts_input <- c("designEV.Input-designCt.Input",
                     "designCX.Input-designCt.Input",
                     "designCM.Input-designCt.Input",
                     "designCM.Input-designEV.Input",
                     "designCM.Input-designCX.Input",
                     "designCX.Input-designEV.Input")

contrasts_poly <- c("designEV.Poly-designCt.Poly",
                    "designCX.Poly-designCt.Poly",
                    "designCM.Poly-designCt.Poly",
                    "designCM.Poly-designEV.Poly",
                    "designCM.Poly-designCX.Poly",
                    "designCX.Poly-designEV.Poly")

transcription_DE <- calculate_DE_include_batches(htseq_matrix_input, contrasts_input)
translation_DE <- calculate_DE_include_batches(htseq_matrix_poly, contrasts_poly)

transcription_new_results <- list()
translation_new_results <- list()
transcription_sig_new_results <- list()
translation_sig_new_results <- list()

for(i in 1:6){
  local_contrast_input <- paste(substr(contrasts_input[i], 7,8), substr(contrasts_input[i], 22,23), sep="_")
  local_contrast_poly <- paste(substr(contrasts_poly[i], 7,8), substr(contrasts_poly[i], 21,22), sep="_")
  transcription_new_results[[local_contrast_input]] <- topTable(transcription_DE,coef=i, n=nrow(htseq_matrix_input),p.value=1)
  translation_new_results[[local_contrast_poly]] <- topTable(translation_DE,coef=i, n=nrow(htseq_matrix_poly),p.value=1)
  
  transcription_sig_new_results[[local_contrast_input]] <- transcription_new_results[[local_contrast_input]][intersect(which(abs(transcription_new_results[[local_contrast_input]]$logFC) > log2(1.1)), which(transcription_new_results[[local_contrast_input]]$adj.P.Val < 0.05)),]
  
  translation_sig_new_results[[local_contrast_poly]] <- translation_new_results[[local_contrast_poly]][intersect(which(abs(translation_new_results[[local_contrast_poly]]$logFC) > log2(1.1)), which(translation_new_results[[local_contrast_poly]]$adj.P.Val < 0.05)),]
  
  n_DE <- nrow(transcription_sig_new_results[[local_contrast_input]])
  n_DT <- nrow(translation_sig_new_results[[local_contrast_poly]])
  print(paste(local_contrast_input, " DE:", n_DE, " DT:", n_DT, sep=""))
}



#Print out DE/DT genes

transcription_new_results <- list()
translation_new_results <- list()
transcription_sig_new_results <- list()
translation_sig_new_results <- list()

gene_classes <- read.csv("gene_classes.txt", header=FALSE, sep=";")
gene_classes <- gene_classes[,1:2]
colnames(gene_classes) <-c("Gene", "Class")
gene_classes$Class <- trimws(gene_classes$Class)

gene_classes_a <- aggregate( .~ Gene, gene_classes, function(x) toString(unique(x)))


for(i in 1:6){
  local_contrast_input <- paste(substr(contrasts_input[i], 7,8), substr(contrasts_input[i], 22,23), sep="_")
  local_contrast_poly <- paste(substr(contrasts_poly[i], 7,8), substr(contrasts_poly[i], 21,22), sep="_")
  transcription_new_results[[local_contrast_input]] <- topTable(transcription_DE,coef=i, n=nrow(htseq_matrix_input),p.value=1, sort.by="p")
  translation_new_results[[local_contrast_poly]] <- topTable(translation_DE,coef=i, n=nrow(htseq_matrix_poly),p.value=1, sort.by="p")
  
  transcription_new_results[[local_contrast_input]]$Gene_class <- gene_classes_a[match(transcription_new_results[[local_contrast_input]]$genes, gene_classes_a$Gene), "Class"]
  translation_new_results[[local_contrast_poly]]$Gene_class <- gene_classes_a[match(translation_new_results[[local_contrast_poly]]$genes, gene_classes_a$Gene), "Class"]
  write.table(transcription_new_results[[local_contrast_input]], file=paste("Sig_genes", local_contrast_input, "input_all.txt", sep="_"), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(translation_new_results[[local_contrast_poly]], file=paste("Sig_genes", local_contrast_poly, "poly_all.txt", sep="_"), row.names=FALSE, quote=FALSE, sep="\t")

}
