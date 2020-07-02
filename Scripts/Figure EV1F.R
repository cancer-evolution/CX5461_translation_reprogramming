###Cluster profiler visualization of resistant GeneGO results

library(ggplot2)
source("helper_functions.R")

##EV Ctrl comparison
genego_up <- read.table("EV_vs_CTRL_diff_translated (w bg list) - ALL SIG up only.txt",
                        sep="\t", header=TRUE)
genego_down <- read.table("EV_vs_CTRL_diff_translated (w bg list) - ALL SIG down only.txt",
                          sep="\t", header=TRUE)
genego_both <- read.table("EV_vs_CTRL_diff (with bg list) - ALL SIG both up and down.txt",
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

genego <- pathway_number_complete[order(pathway_number_complete$FDR, decreasing = FALSE),]
genego <- genego[1:15,]

genego$Ratio <- genego$Number_both/genego$Total

genego$Per_down <- (genego$Number_down/genego$Number_both)*100

genego$Pathways <- factor(genego$Pathways, levels=genego$Pathways[order(genego$FDR, decreasing=TRUE)])

genego$Direction <- genego$Per_down

genego$Direction <- 100-genego$Per_down

ggplot(genego, aes(x=FDR, y =Pathways))+
  geom_point(aes(size=Ratio, fill=Direction),pch=21)+
  scale_fill_gradient2(low="blue", high="red", midpoint=50, mid="white")+
  geom_point(aes(size=Ratio, fill=Direction), pch=21)+
  ylab("")+
  theme_bw()


