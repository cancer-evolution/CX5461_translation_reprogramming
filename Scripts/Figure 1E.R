###Expected signatures in samples

library(GSVA)
library(anota2seq)
library(reshape2)
library(pheatmap)
library(readr)
library(ggplot2)
library(GO.db)

source("polysome_profiling_functions.R")

files <- c("Acute_input.txt","Acute_poly.txt")

sample_types <- c("CTRL", "CX", "EV", "CMB")

counts <- read_counts2(files)


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

##Normalising using voom
preProcess <- anota2seqRNAseqPreProcessing(
  dataP=counts$Poly,     # raw polysome-associated mRNA/RPF counts
  dataT=counts$Input,     # raw total mRNA counts
  transformation = "voom",  # data transformation algorithm
  filterZeroGenes = FALSE # remove genes with 0 counts
)

counts_all_voom <- data.frame(preProcess$dataT, preProcess$dataP)
colnames(counts_all_voom) <- gsub("_cyto", "", colnames(counts_all_voom))
colnames(counts_all_voom) <- gsub("_poly", "", colnames(counts_all_voom))


##Load signatures of interest
load("gene_sets.Rdata")

interest <- as.character(read.delim("GO interest.txt",header=FALSE)$V1)

signature1 <- gene_sets[interest]
names(signature1) <- interest

#####REMOVE NA ONES
interest[sapply(signature1, is.null)] <- NA
interest <- interest[!is.na(interest)]
signature1[sapply(signature1, is.null)] <- NULL
names(signature1) <- interest

terms_not_original_set <- c("anaerobic glutamate catabolic process","glutamate catabolic process via 2-hydroxyglutarate",
                            "mitochondrial electron transport, succinate to ubiquinone", "protein folding",
                            "protein acylation", "protein acetylation", "ncRNA metabolic process", "RNA splicing",
                            "histone modification", "covalent chromatin modification", "DNA damage induced protein phosphorylation")

gene_GO1 <- read_delim("gene_association.mgi", 
                       "\t", escape_double = FALSE, trim_ws = TRUE, col_names=FALSE)
gene_GO <- gene_GO1[gene_GO1[,7] != "IEA",]
gene_GO <- gene_GO[,c(3,5)]
gene_GO <- unique(gene_GO)  #230736
gene_GO <- as.data.frame(gene_GO)
colnames(gene_GO) <- c("Gene", "GOID")

gonames <- Term(GOTERM)
gene_GO$GOname <- unname(gonames[match(gene_GO$GOID, names(gonames))])
load("GOBOFFSPRING_is_a.Rdata")

signature2 <- lapply(terms_not_original_set, function(x){
  id <- names(gonames[gonames == trimws(x)])
  all_offspring <- c(GOBOFFSPRING_is_a[[id]], id)
  
  return(toupper(gene_GO[gene_GO$GOID %in% all_offspring, "Gene"]))
})
names(signature2) <- terms_not_original_set

names_sig2 <- names(signature2)
names_sig2[sapply(signature2, is.null)] <- NA
names_sig2 <- names_sig2[!is.na(names_sig2)]

signature2[sapply(signature2, is.null)] <- NULL
names(signature2) <- names_sig2


signatures <- c(signature1, signature2)

signatures <- lapply(signatures, toupper)

signatures <- signatures[!(names(signatures) %in% c("GO_SODIUM_AMINO_ACID_SYMPORTER_ACTIVITY","GO_AROMATIC_AMINO_ACID_FAMILY_METABOLIC_PROCESS",
                                                    "GO_PYRIMIDINE_CONTAINING_COMPOUND_TRANSMEMBRANE_TRANSPORT","HALLMARK_MTORC1_SIGNALING", "GO_L_AMINO_ACID_IMPORT",
                                                    "GO_AMINO_ACID_TRANSMEMBRANE_TRANSPORTER_ACTIVITY","GO_PEPTIDE_TRANSPORT",
                                                    "GO_POSITIVE_REGULATION_OF_AMINO_ACID_TRANSPORT"))]

upper_gene_names <- toupper(rownames(counts_all_voom))
counts_all_voom$Genes <- upper_gene_names
counts <- aggregate(. ~ Genes, counts_all_voom, mean)
rownames(counts) <- counts$Genes
counts$Genes <- NULL

counts_poly <- counts[,grep("POLY", colnames(counts))]

ssGSEA_signatures <- gsva(as.matrix(counts_poly), signatures, 
                          method="ssgsea", verbose=TRUE, rnaseq=FALSE, ssgsea.norm=TRUE)

ssGSEA_signatures_melt <- melt(ssGSEA_signatures)
colnames(ssGSEA_signatures_melt) <- c("Signature", "Sample", "ssGSEA")
ssGSEA_signatures_melt$Source <- ifelse(grepl("Input", ssGSEA_signatures_melt$Sample, ignore.case=TRUE), "Input", "Poly")
ssGSEA_signatures_melt$Source <- toupper(ssGSEA_signatures_melt$Source)
ssGSEA_signatures_melt$Condition <- substr(ssGSEA_signatures_melt$Sample, 1,2)

ssGSEA_signatures_melt$Condition <- replace(ssGSEA_signatures_melt$Condition, ssGSEA_signatures_melt$Condition=="CM", "CMB")
ssGSEA_signatures_melt$Condition <- replace(ssGSEA_signatures_melt$Condition, ssGSEA_signatures_melt$Condition=="Ct", "CTRL")
ssGSEA_signatures_melt$Condition <- replace(ssGSEA_signatures_melt$Condition, ssGSEA_signatures_melt$Condition=="CT", "CTRL")


#Calculate averages of ssGSEA scores
ssGSEA_signatures_melt_ag <- aggregate(ssGSEA~Signature+Condition+Source, ssGSEA_signatures_melt, mean)

colnames(ssGSEA_signatures_melt_ag) <- c("Pathway", "Treatment", "Source", "ssGSEA_score")

ssGSEA_signatures_melt_ag$ssGSEA_score <- rangecustom(ssGSEA_signatures_melt_ag$ssGSEA_score)
ssGSEA_signatures_melt_ag$Source <- toupper(ssGSEA_signatures_melt_ag$Source)


###Calculate the percentage of difference

ssGSEA_signatures_differences_vsCTRL_per <- calculate_percentage_difference_ssGSEA(ssGSEA_signatures_melt_ag, "CTRL")
ssGSEA_signatures_differences_CXvsCTRL_cast_per <- select_comparison_and_format(ssGSEA_signatures_differences_vsCTRL_per, "CX")
ssGSEA_signatures_differences_CXvsCTRL_cast_per <- data.frame(CX_CTRL = ssGSEA_signatures_differences_CXvsCTRL_cast_per)
rownames(ssGSEA_signatures_differences_CXvsCTRL_cast_per) <- unique(ssGSEA_signatures_melt$Signature)

ssGSEA_signatures_differences_EVvsCTRL_cast_per <- select_comparison_and_format(ssGSEA_signatures_differences_vsCTRL_per, "EV")
ssGSEA_signatures_differences_EVvsCTRL_cast_per <- data.frame(EV_CTRL = ssGSEA_signatures_differences_EVvsCTRL_cast_per)
rownames(ssGSEA_signatures_differences_EVvsCTRL_cast_per) <- unique(ssGSEA_signatures_melt$Signature)

ssGSEA_signatures_differences_CMBvsCTRL_cast_per <- select_comparison_and_format(ssGSEA_signatures_differences_vsCTRL_per, "CMB")
ssGSEA_signatures_differences_CMBvsCTRL_cast_per <- data.frame(CMB_CTRL=ssGSEA_signatures_differences_CMBvsCTRL_cast_per)
rownames(ssGSEA_signatures_differences_CMBvsCTRL_cast_per) <- unique(ssGSEA_signatures_melt$Signature)

temp <- data.frame(ssGSEA_signatures_differences_CXvsCTRL_cast_per,
                   ssGSEA_signatures_differences_EVvsCTRL_cast_per,
                   ssGSEA_signatures_differences_CMBvsCTRL_cast_per)

rownames(temp) <- rownames(ssGSEA_signatures_differences_CMBvsCTRL_cast_per)

real_results <- temp
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "tomato"))(paletteLength)
myBreaks <- c(seq(-14, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(16/paletteLength, 16, length.out=floor(paletteLength/2)))

order <- c("HALLMARK_PI3K_AKT_MTOR_SIGNALING",                                            
           "BIOCARTA_P53_PATHWAY",                                               
           "HALLMARK_P53_PATHWAY",                                                            
           "GO_REGULATION_OF_CELLULAR_SENESCENCE",                                            
           "GO_CELLULAR_SENESCENCE",                                                          
           "HALLMARK_APOPTOSIS",                                                              
           "KEGG_APOPTOSIS",                                                                  
           "REACTOME_AKT_PHOSPHORYLATES_TARGETS_IN_THE_CYTOSOL",                              
           "GO_APOPTOTIC_DNA_FRAGMENTATION",                                                  
           "REACTOME_APOPTOSIS_INDUCED_DNA_FRAGMENTATION",                                    
           "GO_ASPARTATE_FAMILY_AMINO_ACID_METABOLIC_PROCESS",                                
           "GO_BRANCHED_CHAIN_AMINO_ACID_METABOLIC_PROCESS",                                  
           "GO_CELLULAR_AMINO_ACID_BIOSYNTHETIC_PROCESS",                                     
           "GO_CELLULAR_AMINO_ACID_CATABOLIC_PROCESS",                                        
           "GO_GLUTAMINE_FAMILY_AMINO_ACID_METABOLIC_PROCESS",                                
           "GO_SULFUR_AMINO_ACID_METABOLIC_PROCESS",                                          
           "lipid catabolic process",                                                         
           "lipid biosynthetic process",                                                      
           "triglyceride metabolic process",                                                  
           "glycolytic process",                                                              
           "pyruvate metabolic process",                                                      
           "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",                                            
           "ATP biosynthetic process",                                                        
           "oxidative phosphorylation",                                                       
           "respiratory electron transport chain",                                            
           "mitochondrial electron transport, NADH to ubiquinone", 
           "mitochondrial electron transport, succinate to ubiquinone",
           "mitochondrial electron transport, ubiquinol to cytochrome c",                     
           "mitochondrial electron transport, cytochrome c to oxygen",                        
           "rRNA processing",                                                                 
           "oxidative RNA demethylation",                                                     
           "GO_RIBOSOME_BIOGENESIS",                                                          
           "translation",                                                                     
           "REACTOME_PEPTIDE_CHAIN_ELONGATION",                                               
           "mitochondrial translation",                                                       
           "lipid oxidation",                                                                 
           "protein oxidation",                                                               
           "ubiquitin-dependent protein catabolic process",                                   
           "protein methylation",                                                             
           "GO_REGULATION_OF_TRANSCRIPTION_INVOLVED_IN_G1_S_TRANSITION_OF_MITOTIC_CELL_CYCLE",
           "GO_REGULATION_OF_CELL_CYCLE_G1_S_PHASE_TRANSITION",                               
           "GO_POSITIVE_REGULATION_OF_CELL_CYCLE_ARREST",                                     
           "GO_MITOTIC_CELL_CYCLE_CHECKPOINT",                                                
           "GO_NEGATIVE_REGULATION_OF_CELL_CYCLE",                                            
           "GO_POSITIVE_REGULATION_OF_CELL_CYCLE",                                            
           "GO_CELL_CYCLE_PHASE_TRANSITION",                                                  
           "PID_ATM_PATHWAY",                                                                 
           "KEGG_CELL_CYCLE",                                                                 
           "GO_MITOTIC_G2_M_TRANSITION_CHECKPOINT",                                           
           "GO_G2_DNA_DAMAGE_CHECKPOINT",                                                     
           "GO_SPINDLE_CHECKPOINT", 
           "GO_PYRIMIDINE_CONTAINING_COMPOUND_BIOSYNTHETIC_PROCESS",                                           
           "GO_PYRIMIDINE_CONTAINING_COMPOUND_CATABOLIC_PROCESS",    
           "purine nucleobase biosynthetic process",                                          
           "purine nucleobase metabolic process",
           "DNA repair",                                                                      
           "protein folding",                                                                 
           "protein acylation",                                                               
           "protein acetylation",                                                             
           "ncRNA metabolic process",                                                        
           "RNA splicing",                                                                    
           "histone modification",                                                            
           "covalent chromatin modification",                                                 
           "DNA damage induced protein phosphorylation")   

real_results <- real_results[match(order, rownames(real_results)),]

pheatmap(real_results,
         color=myColor, breaks=myBreaks, fontsize=6, cluster_cols=FALSE, cluster_rows=FALSE)


