##Script to find DEG/DET genes using anota/limma
#Do PCA plot to being with

library(ggplot2)
library(limma)
library(edgeR)
library(anota2seq)
library(DESeq2)
library(ggrepel)

source("polysome_profiling_functions.R")

files <- c("Acute_input.txt",
           "Acute_poly.txt")

sample_types <- c("CMB", "CX", "EV", "CTRL")
contrasts <- c("CMB_CX", "CMB_CTRL", "CMB_EV")

omni = TRUE

#First column should be gene names, other columns counts
#The names of counts columns are assumed to be of the format:
#CTRL-sample_number-INPUT

use_only_coding_genes <- TRUE

counts <- read_counts2(files)
  
if(use_only_coding_genes == TRUE){
  gene_classes <- read.csv("gene_classes.txt", header=FALSE, sep=";")
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

if(length(sample_types_input2) == length(sample_types_input2)){
  print("Ok. Input and poly number of samples are of the same!")
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

colnames(counts$Poly)
colnames(counts$Input)

##Normalising using voom
preProcess <- anota2seqRNAseqPreProcessing(
  dataP=counts$Poly,     # raw polysome-associated mRNA/RPF counts
  dataT=counts$Input,     # raw total mRNA counts
  transformation = "voom",  # data transformation algorithm
  filterZeroGenes = TRUE # remove genes with 0 counts
)

counts_all_voom <- data.frame(preProcess$dataT, preProcess$dataP)
colnames(counts_all_voom) <- gsub("_cyto", "", colnames(counts_all_voom))
colnames(counts_all_voom) <- gsub("_poly", "", colnames(counts_all_voom))



#####Anota analysis with anota2seq following the step-by-step method

sample_classes <- rep(NA, length(ncol(preProcess$dataT)))
for(local_type in sample_types){
  sample_classes[grep(local_type, colnames(preProcess$dataT))] <- local_type
}

qc <- anota2seqPerformQC(dataT= preProcess$dataT,
                         dataP= preProcess$dataP,
                         phenoVec = sample_classes,
                         useProgBar = FALSE)

#########OMNIBUS
##Select set of omnibus genes, and then re-run anota

# Get the genes that have a FDR < 0.15 from the omnibus test
omniSig <- qc$omniGroupStats[qc$omniGroupStats[,"groupRvmPAdj"] < .15,] 

omni_genes <- rownames(omniSig)

write.table(omni_genes, file="Omni_genes.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)

counts_omni <- counts_all_voom[rownames(counts_all_voom) %in% omni_genes, ]

sd_omni <- apply(counts_omni,1,sd)

plot(density(sd_omni))

counts_omni_filt <- counts_omni[sd_omni > .3,]


pca_omni_filt<- generatePCA(count_data = counts_omni_filt,nameApp = "OMNI")

counts_omni_list <- counts
counts_omni_list$Input <- counts$Input[rownames(counts$Input) %in% omni_genes,]
counts_omni_list$Poly <- counts$Poly[rownames(counts$Poly) %in% omni_genes,]

##Normalising using voom
preProcess_omni <- anota2seqRNAseqPreProcessing(
  dataP=counts_omni_list$Poly,     # raw polysome-associated mRNA/RPF counts
  dataT=counts_omni_list$Input,     # raw total mRNA counts
  transformation = "voom",  # data transformation algorithm
  filterZeroGenes = TRUE # remove genes with 0 counts
)


all.equal(dim(preProcess_omni$dataP),dim(preProcess_omni$dataT))

sample_classes_omni <- NA
for(local_type in sample_types){
  sample_classes_omni[grep(local_type, colnames(preProcess_omni$dataP), ignore.case=TRUE)] <- local_type
}

qc_omni <- anota2seqPerformQC(dataT= preProcess_omni$dataT,
                              dataP= preProcess_omni$dataP,
                              phenoVec = sample_classes_omni,
                              useProgBar = TRUE)

residOutlier_omni <- anota2seqResidOutlierTest(anota2seqQcObj = qc_omni,
                                               useProgBar = T)

n_contrasts <- length(contrasts)

contrast_matrix <- matrix(0, ncol=n_contrasts,
                          nrow=length(qc$phenoClasses))
rownames(contrast_matrix) <- qc$phenoClasses
colnames(contrast_matrix) <- contrasts

for(loc_con in contrasts){
  loc_con_treatment <- unlist(strsplit(loc_con, "_"))[1]
  loc_con_control <- unlist(strsplit(loc_con, "_"))[2]
  contrast_matrix[rownames(contrast_matrix) == loc_con_treatment,
                  colnames(contrast_matrix) == loc_con] <- 1
  
  contrast_matrix[rownames(contrast_matrix) == loc_con_control,
                  colnames(contrast_matrix) == loc_con] <- -1
}

print("Contrast matrix")
print(contrast_matrix)

## differential expression analysis polysome-associated mRNA
anota2seqPolyOut_omni <- anota2seqAnalyse(
  dataT = preProcess_omni$dataT,
  dataP = preProcess_omni$dataP,
  phenoVec = sample_classes_omni,
  analyse = "polysomeassociated mRNA",
  contrasts = contrast_matrix,
  useProgBar = FALSE)

## differential expression analysis total mRNA
anota2seqTotalOut_omni <- anota2seqAnalyse(
  dataT = preProcess_omni$dataT,
  dataP = preProcess_omni$dataP,
  phenoVec = sample_classes_omni,
  analyse = "total mRNA",
  contrasts = contrast_matrix,
  useProgBar = FALSE)

## analysis of differential translation
anota2seqTranslationOut_omni <- anota2seqAnalyse(
  dataT = preProcess_omni$dataT,
  dataP = preProcess_omni$dataP,
  phenoVec = sample_classes_omni,
  analyse = "translation",
  contrasts = contrast_matrix,
  useProgBar = FALSE)

## analysis of translational buffering
anota2seqBufferingOut_omni <- anota2seqAnalyse(
  dataT = preProcess_omni$dataT,
  dataP = preProcess_omni$dataP,
  phenoVec = sample_classes_omni,
  analyse = "buffering",
  contrasts = contrast_matrix,
  useProgBar = FALSE)

combineOutput_omni <- list("qualityControl" = NULL,
                           "residOutlierTest" = NULL,
                           "polysomeassociatedmRNA" =anota2seqPolyOut_omni,
                           "totalmRNA" = anota2seqTotalOut_omni,
                           "translation" = anota2seqTranslationOut_omni,
                           "buffering" = anota2seqBufferingOut_omni)


for(i in 1:length(contrasts)){
  contrast <- contrasts[i]
  
  temp_poly <- as.data.frame(anota2seqPolyOut_omni$apvStatsRvm[[i]])
  temp_poly <- temp_poly[order(temp_poly$apvRvmPAdj, decreasing=FALSE),]
  write.table(temp_poly, file=paste("Anota_poly_results_",contrast, "_omni.txt", sep=""),
              quote=FALSE, sep="\t")
  
  temp_input <- as.data.frame(anota2seqTotalOut_omni$apvStatsRvm[[i]])
  temp_input <- temp_input[order(temp_input$apvRvmPAdj, decreasing=FALSE),]
  write.table(temp_input, file=paste("Anota_input_results_",contrast, "_omni.txt", sep=""),
              quote=FALSE, sep="\t")
  
  temp_trans <- as.data.frame(anota2seqTranslationOut_omni$apvStatsRvm[[i]])
  temp_trans <- temp_trans[order(temp_trans$apvRvmPAdj, decreasing=FALSE),]
  write.table(temp_trans, file=paste("Anota_differential_translation_results_",contrast, "_omni.txt", sep=""),
              quote=FALSE, sep="\t")
  
  temp_buff <- as.data.frame(anota2seqBufferingOut_omni$apvStatsRvm[[i]])
  temp_buff <- temp_buff[order(temp_buff$apvRvmPAdj, decreasing=FALSE),]
  write.table(temp_buff, file=paste("Anota_buffering_results_",contrast, "_omni.txt", sep=""),
              quote=FALSE, sep="\t")
  print(contrast)
  
  df <- data.frame(Gene = rownames(temp_input), Input = temp_input$apvEff)
  df <- data.frame(df, Poly=temp_poly$apvEff[match(rownames(temp_input), rownames(temp_poly))])
  
  diff_translated <- rownames(temp_trans)[intersect(which(abs(temp_trans$apvEff) > log2(1.25)),
                                                    which(temp_trans$apvRvmPAdj < 0.15))]
  buffered <- rownames(temp_buff)[intersect(which(abs(temp_buff$apvEff) > log2(1.25)),
                                            which(temp_buff$apvRvmPAdj < 0.15))]
  sig_input <- rownames(temp_input)[intersect(which(abs(temp_input$apvEff) > log2(1.25)),
                                              which(temp_input$apvRvmPAdj < 0.15))]
  sig_poly <- rownames(temp_poly)[intersect(which(abs(temp_poly$apvEff) > log2(1.25)),
                                            which(temp_poly$apvRvmPAdj < 0.15))]
  
  df$Class <- "None"
  df$Class[df$Gene %in% buffered] <- "Buffered"
  df$Class[df$Gene %in% diff_translated] <- "Differentially_translated"
  
  maximum_value <- max(c(df$Input, df$Poly))
  
  
  df[df$Class == "Differentially_translated","Class"] <- ifelse(df[df$Class == "Differentially_translated","Poly"] >0, "Diff_trans_up",
                                                                "Diff_trans_down")
  
  df[df$Class == "Buffered","Class"] <- ifelse(df[df$Class == "Buffered","Input"] >0, "Buff_down",
                                               "Buff_up")
  
  mRNA <- sig_input[sig_input %in% sig_poly]
  
  df$Class[df$Gene %in% mRNA] <- "mRNA"
  
  df[df$Class == "mRNA","Class"] <- ifelse(df[df$Class == "mRNA","Input"] >0, "mRNA_up", "mRNA_down")
  
  
  pdf(paste("Total mRNA - Polysome plot ", contrast, ".pdf", sep=""), height=6, width=8)
  
  g <- ggplot(df, aes(x=Input, y =Poly))+
    geom_point(aes(fill=Class), pch=21)+
    geom_vline(xintercept=0, linetype="dashed")+
    geom_hline(yintercept=0, linetype="dashed")+
    geom_vline(xintercept=log2(1.25))+
    geom_hline(yintercept=log2(1.25))+
    geom_vline(xintercept=-log2(1.25))+
    geom_hline(yintercept=-log2(1.25))+
    geom_abline(intercept=0, slope=1, linetype="dashed")+
    geom_abline(intercept=log2(1.25), slope=1)+
    geom_abline(intercept=-log2(1.25), slope=1)+
    scale_fill_manual(values=c(Diff_trans_up = "orange",
                               Diff_trans_down = "darkred",
                               Buff_down="royalblue", Buff_up="lightblue",
                               None="black", mRNA_up="lightgreen",
                               mRNA_down="darkgreen"))+
    ggtitle(contrast)+
    xlab("Total mRNA")+
    ylab("Polysome-associated mRNA")+
    xlim(-maximum_value,maximum_value)+
    ylim(-maximum_value,maximum_value)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  print(g)
  dev.off()
}
  