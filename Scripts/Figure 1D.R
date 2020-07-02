##Visualising Metacore results

library(ggplot2)

input_elongation_list <- list()
input_initiation_list <- list()

for(contrast in c("CMB_CTRL", "EV_CTRL")){
  
  if(contrast == "CMB_CTRL"){
    input_elongation <- read.table("Translation_elongation_both 29 out of 133 CMB vs_CTRL.txt",header=TRUE)
    input_initiation <- read.table("Translation_initiation_both 31 out of 146 CMB vs_CTRL 2nd.txt", header=TRUE)
  }else if(contrast=="EV_CTRL"){
    input_elongation <- read.table("Translation_elongation_both out of 145 EV vs CTRL.txt",header=TRUE)
    input_initiation <- read.table("Translation_initiation_both out of 145 EV vs CTRL.txt", header=TRUE)
  }
  
  rownames(input_elongation) <- input_elongation$Gene.Symbol
  input_elongation <- input_elongation[,c("Gene.Symbol", "Signal", "p.value")]
  colnames(input_elongation) <- c("Genes", "logFC", "adj.P.Val")
  
  genes <- rownames(input_elongation)
  
  input_elongation$adj.P.Val_plot <- 1/input_elongation$adj.P.Val
  input_elongation$Significance_p <- input_elongation$adj.P.Val_plot
  input_elongation$Significance_p[input_elongation$adj.P.Val_plot < 1/0.05] <- NA
  
  input_elongation <- input_elongation[!is.na(input_elongation$Significance_p),]
  
  input_elongation$Genes <- factor(input_elongation$Genes, levels=rev(input_elongation$Genes))
  input_elongation$Translation <- "Elongation"
  
  
  rownames(input_initiation) <- input_initiation$Gene.Symbol
  input_initiation <- input_initiation[,c("Gene.Symbol", "Signal", "p.value")]
  colnames(input_initiation) <- c("Genes", "logFC", "adj.P.Val")
  
  genes <- rownames(input_initiation)
  
  input_initiation$adj.P.Val_plot <- 1/input_initiation$adj.P.Val
  input_initiation$Significance_p <- input_initiation$adj.P.Val_plot
  input_initiation$Significance_p[input_initiation$adj.P.Val_plot < 1/0.05] <- NA
  
  input_initiation <- input_initiation[!is.na(input_initiation$Significance_p),]
  
  input_initiation$Translation <- "Initiation"
  
  input_initiation$Contrast <- contrast
  
  input_elongation$Contrast <- contrast
  
  input_initiation_list[[contrast]] <- input_initiation
  input_elongation_list[[contrast]] <- input_elongation
  
}



input_initiation_both <- rbind(input_initiation_list$CMB_CTRL, input_initiation_list$EV_CTRL)
input_elongation_both <- rbind(input_elongation_list$CMB_CTRL, input_elongation_list$EV_CTRL)

input_initiation_both$Genes <- factor(input_initiation_both$Genes, levels=sort(unique(input_initiation_both$Genes)))
input_elongation_both$Genes <- factor(input_elongation_both$Genes, levels=sort(unique(input_elongation_both$Genes)))


input_initiation_both$Translation <- as.character(input_initiation_both$Translation)
input_initiation_both$Contrast <- as.character(input_initiation_both$Contrast)
input_initiation_both$Genes <- as.character(input_initiation_both$Genes)
input_initiation_both$Genes <- factor(input_initiation_both$Genes)

input_initiation_both$Genes <- factor(input_initiation_both$Genes,
                                      level=c("Eif1","Eif1a","Eif2s2","Eif3f","Eif4e","Gm4705","Rpl11","Rpl14",
                                              "Rpl17","Rpl22","Rpl23","Rpl24","Rpl26","Rpl27a","Rpl30","Rpl34",   
                                              "Rpl34-ps1","Rpl35","Rpl35a","Rpl37","Rpl6","Rps13","Rps15a","Rps16",    
                                              "Rps18","Rps19","Rps21","Rps25","Rps3a1","Rps6","Rps8",
                                              "Rps9","Stau1","Rps24"))

input_initiation_both$adj.P.Val_plot <- as.numeric(as.character(input_initiation_both$adj.P.Val_plot))
input_initiation_both$Significance_p <- as.numeric(as.character(input_initiation_both$Significance_p))
input_initiation_both$logFC <- as.numeric(as.character(input_initiation_both$logFC))

#Translation initiation
ggplot(input_initiation_both, aes(x=1, y=Genes))+
  geom_point(aes(size=adj.P.Val_plot, colour=logFC))+
  geom_point(aes(size=Significance_p), pch=21, stroke = 1.25)+
  scale_size("adj. p-value", labels=c(0.02, 0.04, 0.1, 1), 
             breaks=c(1/0.02, 1/0.04, 1/0.1, 1))+
  scale_colour_gradient2(low="blue", mid="white", high="red")+
  theme_bw()+
  ylab("")+
  ggtitle("Translation initiation")+
  facet_grid(.~Contrast, scale="free_y")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))


##Translation elongation
ggplot(input_elongation_both, aes(x=1, y=Genes))+
  geom_point(aes(size=adj.P.Val_plot, colour=logFC))+
  geom_point(aes(size=Significance_p), pch=21, stroke = 1.25)+
  #geom_point(aes(size=adj.P.Val_plot, colour=logFC))+
  scale_size("adj. p-value", labels=c(0.0005, 0.005, 0.05, 0.5, 1), breaks=c(1/0.0005, 1/0.005, 1/0.05, 1/0.5, 1))+
  scale_colour_gradient2(low="blue", mid="white", high="red")+
  theme_bw()+
  ylab("")+
  ggtitle("Translation elongation")+
  facet_grid(.~Contrast, scale="free_y")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

