calculate_DE_include_batches <- function(counts, contrasts){
  y <- DGEList(counts = counts, genes = rownames(counts))
  
  #Keep genes with a count of more than 25
  A <- rowSums(y$counts)
  isexpr <- A > 25

  y <- y[isexpr, , keep.lib.size = FALSE]
  y <- calcNormFactors(y)
  
  y$samples$lib.size <- colSums(y$counts)
  
  sample_types <- substr(colnames(counts), 1,2)
  method <- unname(sapply(colnames(counts), function(x){
    unlist(strsplit(unlist(strsplit(x, "\\."))[2], "_"))[1]
  }))
  
  design_names <- paste(sample_types, method, sep=".")
  
  design <- factor(design_names, unique(design_names))
  
  batches <- rep(NA, length(colnames(counts)))
  batches <- replace(batches, grep(paste(c("CX14.Poly", "CX39.Poly"), collapse="|"), colnames(counts)), "N")
  batches <- replace(batches, is.na(batches), "O")
  batches <- as.factor(batches)
  
  if(length(levels(batches)) == 1){
    design_matrix <- model.matrix(~ 0+design)
  }else{
    design_matrix <- model.matrix(~ 0+design+batches)
    
  }
  
  v <- voom(y,design_matrix,plot=TRUE)
  
  fit <- lmFit(v,design_matrix)
  
  astr=paste(contrasts, collapse=",")
  commandstr=paste("makeContrasts(",astr,", levels=design_matrix)",sep="")
  contrast.matrix <- eval(parse(text=commandstr))
  
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  
  fit2 <- eBayes(fit2)
  
  return(fit2)
}


trim <- function (x) gsub("^\\s+|\\s+$", "", x)


read_counts2 <- function(files){

  txtList = lapply(files, function(x){ read.delim(x)} ) 
  
  txtList = lapply(txtList, function(x){
    colnames(x) <- as.character(unname(unlist(x[1,])))
    x <- x[-1,]
    rownames(x) <- x[,1]
    x <- x[,-1]
    return(x)
  })
  
  #Since not all files have the same genes
  #I will complete the matrices, adding counts of zero
  #in cases when a gene is not expressed at all in only one of the input types
  
  all_genes <- sort(unique(do.call("c", sapply(txtList, rownames))))
  all_genes <- all_genes[!(all_genes %in% c("__ambiguous", "__no_feature"))]
  
  txtList <- lapply(txtList, function(x){
    x <- x[match(all_genes, rownames(x)),]
    x[is.na(x)] <- 0
    return(x)
  })
  
  matrix_input <- do.call("cbind", lapply(txtList, function(x){
    return(x[,grep("Input", colnames(x), ignore.case=TRUE)])
  }))
  
  matrix_poly <- do.call("cbind", lapply(txtList, function(x){
    return(x[,grep("Poly", colnames(x), ignore.case=TRUE)])
  }))

  matrix_input2 <- apply(matrix_input, 2, function(x){ as.numeric(as.character(x))})
  matrix_poly2 <- apply(matrix_poly, 2, function(x){ as.numeric(as.character(x))})
  
  rownames(matrix_input2) <- all_genes
  rownames(matrix_poly2) <- all_genes
  
  print("Input matrix column names")
  print(colnames(matrix_input2))
  
  print("Poly matrix column names")
  print(colnames(matrix_poly2))
  
  return(list(Input = matrix_input2,
              Poly = matrix_poly2))
}


calculate_percentage_difference_ssGSEA <- function(melt_ag, reference){
  differences <- vector()
  for(path in unique(melt_ag$Pathway)){
    temp <- melt_ag[melt_ag$Pathway == path, ]
    temp_input <- temp[temp$Source == "INPUT",]
    temp_poly <- temp[temp$Source == "POLY",]
    
    temp_input_remove <- temp_input[temp_input$Treatment == reference, "ssGSEA_score"]
    temp_input$ssGSEA_score <- ((temp_input$ssGSEA_score - temp_input_remove)/temp_input_remove)*100
    temp_input <- temp_input[temp_input$Treatment != reference,]
    
    temp_poly_remove <- temp_poly[temp_poly$Treatment == reference, "ssGSEA_score"]
    temp_poly$ssGSEA_score <- ((temp_poly$ssGSEA_score - temp_poly_remove)/temp_poly_remove)*100
    temp_poly <- temp_poly[temp_poly$Treatment != reference,]
    
    differences <- rbind(differences, temp_input, temp_poly)
  }
  return(differences)
}


generatePCA <- function(count_data,nameApp=""){
  pca <- prcomp(t(count_data))
  pca_ggplot <- as.data.frame(pca$x[,1:4])
  pca_ggplot$Samples <- rownames(pca_ggplot)
  
  pca_ggplot$Source <- sapply(strsplit(pca_ggplot$Samples,"\\."), `[`, 3)
  
  pca_ggplot$Type <- NA
  for(local_type in sample_types){
    pca_ggplot$Type[grep(local_type, pca_ggplot$Samples, ignore.case=TRUE)] <- local_type
  }
  
  pdfName <- paste("PCA_",nameApp,".pdf",sep="")
  
  pdf(pdfName, height=5)
  par(mfrow=c(2,2))
  g1 <- ggplot(pca_ggplot, aes(x=PC1, y=PC2))+
    geom_point(aes(colour=Type, shape=Source))+
    geom_text_repel(aes(label=Samples), size=2)+
    theme_bw()
  g2 <- ggplot(pca_ggplot, aes(x=PC1, y=PC3))+
    geom_point(aes(colour=Type, shape=Source))+
    geom_text_repel(aes(label=Samples), size=2)+
    theme_bw()
  g3 <- ggplot(pca_ggplot, aes(x=PC2, y=PC3))+
    geom_point(aes(colour=Type, shape=Source))+
    geom_text_repel(aes(label=Samples), size=2)+
    theme_bw()
  g4 <- ggplot(pca_ggplot, aes(x=PC1, y=PC4))+
    geom_point(aes(colour=Type, shape=Source))+
    geom_text_repel(aes(label=Samples), size=2)+
    theme_bw()
  print(g1)
  print(g2)
  print(g3)
  print(g4)
  dev.off()
  return(pca)
}


rangecustom <- function(x){ (x - min(x))/(max(x)-min(x)) * (10 - 1) + 1 }


plot_p_logFC <- function(genes_df, genes){
  
  genes_df$Genes <- c(genes, genes)
  genes_df$Genes <- factor(genes_df$Genes, levels=rev(genes))
  
  genes_df$adj.P.Val_plot <- 1/genes_df$adj.P.Val
  genes_df$Significance_p <- genes_df$adj.P.Val_plot#+(genes_df$adj.P.Val_plot*3)
  genes_df$Significance_p[genes_df$adj.P.Val_plot < 1/0.05] <- NA
  
  g <- ggplot(genes_df, aes(x=1, y=Genes))+
    geom_point(aes(size=adj.P.Val_plot, colour=logFC))+
    geom_point(aes(size=Significance_p), pch=21, stroke = 1.25)+
    #geom_point(aes(size=adj.P.Val_plot, colour=logFC))+
    scale_size("adj. p-value", labels=c(0.0005, 0.005, 0.05, 0.5, 1), breaks=c(1/0.0005, 1/0.005, 1/0.05, 1/0.5, 1))+
    scale_colour_gradient2(low="blue", mid="white", high="red")+
    theme_bw()+
    ylab("")+
    facet_grid(.~Comparison, scale="free_y")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))
  return(g)
}