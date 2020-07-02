###Cluster profiler visualization of resistant GeneGO results

library(ggplot2)
source("helper_functions.R")


genego_up <- read.table("CX_vs_CTRL_diff_translated (w bg list) - ALL SIG up only.txt",
                        sep="\t", header=TRUE)
genego_down <- read.table("CX_vs_CTRL_diff_translated (w bg list) - ALL SIG down only.txt",
                          sep="\t", header=TRUE)
genego_both <- read.table("CX_vs_CTRL_diff (with bg list) - ALL SIG both up and down.txt",
                          sep="\t", header=TRUE)


genego_up <- genego_up[,c("Networks", "FDR", "In.Data", "Total")]
genego_down <- genego_down[,c("Networks", "FDR", "In.Data", "Total")]
genego_both <- genego_both[,c("Networks", "FDR", "In.Data", "Total")]

all_pathways <- unique(c(as.character(genego_up$Networks), as.character(genego_down$Networks), 
                         as.character(genego_both$Networks)))
pathway_number <- data.frame(Pathways = all_pathways,
                             Number_both = genego_both[match(all_pathways, genego_both$Networks), "In.Data"],
                             Number_up = genego_up[match(all_pathways, genego_up$Networks), "In.Data"],
                             Number_down = genego_down[match(all_pathways, genego_down$Networks), "In.Data"])

pathway_number_complete <- vector()
for(i in 1:nrow(pathway_number)){
  local_row <- unlist(pathway_number[i,])
  
  if(is.na(local_row[3])){
    local_row[3] <- local_row[2]-local_row[4]
  }else if(is.na(local_row[4])){
    local_row[4] <- local_row[2]-local_row[3]
  }
  pathway_number_complete <- rbind(pathway_number_complete, local_row)
}

pathway_number_complete <- as.data.frame(pathway_number_complete)
pathway_number_complete$Pathways <- pathway_number$Pathways

pathway_number_complete <- pathway_number_complete[!is.na(pathway_number_complete$Number_both),]

pathway_number_complete$Total <- genego_both[match(pathway_number_complete$Pathways, genego_both$Networks), "Total"]
pathway_number_complete$FDR <- genego_both[match(pathway_number_complete$Pathways, genego_both$Networks), "FDR"]

pathway_number_complete <- pathway_number_complete[pathway_number_complete$Number_up >= 0,]

genego <- pathway_number_complete[order(pathway_number_complete$FDR, decreasing = FALSE),]

genego$Ratio <- genego$Number_both/genego$Total

genego$Per_down <- (genego$Number_down/genego$Number_both)*100

genego$Pathways <- factor(genego$Pathways, levels=genego$Pathways[order(genego$FDR, decreasing=TRUE)])

genego$Direction <- genego$Per_down

##If only one direction

genego_up <- read.table("CX-res vs Ctrl naive (SIG Up only FC 1.5).txt",
                        sep="\t", header=TRUE)

genego_up <- genego_up[,c("Networks", "FDR", "In.Data", "Total")]
genego_up$Ratio <- genego_up$In.Data/genego_up$Total
colnames(genego_up)[1] <- "Pathways"

genego_up$Pathways <- factor(genego_up$Pathways, levels=rev(genego_up$Pathways))

top_10 <- genego_up$Pathways[order(genego_up$FDR, decreasing=FALSE)][1:10]

genego_up_top_10 <- genego_up[genego_up$Pathways %in% top_10,]

g <- ggplot(genego_up_top_10, aes(x=FDR, y =Pathways))+
  geom_point(aes(size=Ratio, fill=FDR),pch=21)+
  scale_fill_gradient(low="red", high="blue")+
  ylab("")+
  theme_bw()
print(g)

