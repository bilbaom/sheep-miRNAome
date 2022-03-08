setwd("sheep-miRNAome")
library("edgeR")
library("pheatmap")
library("ggplot2"
library("Rtsne")
library("UpSetR")
#The vector of colors for the plots 
c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2", 
         "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", 
         "green1", "yellow4", "yellow3", "darkorange4", "brown")

############################################################################################################################################################################################
# Introduce the results of miRDeep2 quantifier tool and the information table of samples in two variables
miRNACount <- read.csv("miRNAs_expressed_all_samples_now.csv", sep = "\t") # The quantified expression matrix
info <- read.csv("sample_info.csv", row.names = 1, stringsAsFactors = FALSE, colClasses = rep("character", 6)) # The samples information
colnames(miRNACount) <- gsub("X", "", colnames(miRNACount))

# The novel and known miRNAs are separated
knowmiRNAC <- miRNACount[grep("oar", miRNACount$.miRNA), ]
newmiRNAc <- miRNACount[!grepl("oar", miRNACount$.miRNA), ]
  
# Multicopies in novel miRNAs are removed 
allmirnas <- as.data.frame(read.table("all_mirnas.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
allmirnas <- allmirnas[allmirnas["V10"] == "q",]
allmirnasF <- data.frame(do.call('rbind', strsplit(as.character(allmirnas$V4),'.',fixed=TRUE)))
allmirnasF <- as.vector(allmirnasF[['X1']])
allmirnas["V11"] <- allmirnasF
newmiRNAcF <-  data.frame(do.call('rbind', strsplit(as.character(newmiRNAc$.miRNA),'.',fixed=TRUE)))
newmiRNAc["pref"] <- newmiRNAcF$X1
newmiRNAc2 <- subset(newmiRNAc, pref %in% allmirnasF)

# Change the names of novel miRNAs and set the correct suffix 
newmiRNAc2["miRNA"] <- NA

for (i in allmirnas$V11) {    
    t <- newmiRNAc2$pref == i 
    p <- allmirnas$V11 == i
    n <- allmirnas[p, "V7"]
    newmiRNAc2[t, "miRNA"] <- n
   }
    
newmiRNAc2["subf"] <- NA

for (i in 1:nrow(newmiRNAc2)) {
    if (endsWith(as.character(newmiRNAc2[i, ".miRNA"]), "-5p")) {
    newmiRNAc2[i, "subf"] <- "-5p"
    } else if (endsWith(as.character(newmiRNAc2[i, ".miRNA"]), "-3p")) {
    newmiRNAc2[i, "subf"] <- "-3p"
    }
   }

newmiRNAc2["miRNA2"] <- paste(newmiRNAc2$miRNA, newmiRNAc2$subf, sep = "")
rownames(newmiRNAc2) <- newmiRNAc2$miRNA2
miRNACNW <- newmiRNAc2[info$Code]
    
# In know miRNAs there are 2 multicopies: oar-miR-181a and oar-miR-29b. We take only one of them.
row.names(knowmiRNAC) <- NULL
repeats <- knowmiRNAC[knowmiRNAC$.miRNA == "oar-miR-29b", ]
repeats <- rbind(repeats, knowmiRNAC[knowmiRNAC$.miRNA == "oar-miR-181a", ])
knowmiRNAC <- knowmiRNAC[-c(37, 59), ] # Those two numbers are the row coordinates of the miRNAs that we have to remove
rownames(knowmiRNAC) <- knowmiRNAC$.miRNA
miRNACKN <- knowmiRNAC[info$Code]
   
# The expression of known and novel miRNAs are joined and the definitive names for the  samples are set
miRNACount2 <- rbind(miRNACNW, miRNACKN)
colnames(miRNACount2) <- rownames(info)

# miRNAs that do not have 10 or more reads on average in at least one tissue are removed
Tissue_average <- data.frame(row.names = rownames(miRNACount2))
for(i in unique(info$Tissue)) { 
  samples <- rownames(info[info["Tissue"] == i,])
  row_means <- rowMeans(miRNACount2[ ,samples])
  Tissue_average[i] <- row_means
  }

x <- miRNACount2[rowSums(Tissue_average > 10) > 0,]
expressionMatrix <- unique(x)

# Normalized expression    
normalizedEM <- cpm(expressionMatrix) 
        
############################################################################################################################################################################################
#Correlation heatmap of the tissues:
c21 <- c25[c(1:21)]
names(c21) <- unique(info["Tissue"])$Tissue
plot_colors <- list(Tissue = c21)
sample_correlation_log <- cor(log(normalizedEM + 1))

while (!is.null(dev.list())) dev.off()
pdf(file = "plots/correlation_heatmap.pdf", width = 7, height = 5.5)
pheatmap(sample_correlation_log, 
         annotation = info["Tissue"],
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         scale = "none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_colors = plot_colors)
dev.off()
  
############################################################################################################################################################################################
#tSNE
Tissues <- info$Tissue
transp_normalizedEM <- t(normalizedEM)
log_transp_normalizedEM <- (log2(transp_normalizedEM + 0.1))
tsne <- Rtsne(log_transp_normalizedEM, dim = 2, perplexity = 20, verbose = TRUE, check_duplicates = FALSE)
tsne_color <- data.frame(x = tsne$Y[ ,1], y = tsne$Y[ ,2], Tissues = info$Tissue) 

pdf(file = "plots/tSNE.pdf", width = 6.5, height = 4)
ggplot(tsne_color, aes(x, y, colour = Tissues)) + geom_point() + scale_color_manual(values = c21) + labs(x = "t-SNE-1", y = "t-SNE-2") + 
      theme_bw() + 
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) 
    dev.off()
        
############################################################################################################################################################################################
#UPSET:

filtr <- function(counts, tissue) {
  samples <- rownames(info[info["Tissue"] == tissue,])
  return(counts[ ,samples])
}

count_expressed <- function(counts, cpm)
   tissue_vector <- c(rep(0, length(rownames(counts))))
   for(i in 1:nrow(counts)) {
   row <- counts[i,]
   if (sum(as.numeric(row) > cpm) >= (length(as.data.frame(counts))/2)) {
   tissue_vector[i] <- 1
    }   
  }
   return(tissue_vector)
}

upset_df <- data.frame(normalizedEM)

for(i in unique(info$Tissue)) {
  upset_df[i] <- count_expressed(filtr(normalizedEM, i), 5)
}

upset_df <- upset_df[ ,as.character(unique(info$Tissue))]

while (!is.null(dev.list())) dev.off()
pdf(file = "plots/Upset_5.pdf", width = 7, height = 3.5)
upset(upset_df, 
      order.by = "freq", 
      nsets = 21, 
      nintersects = 30, # The number of intersects have to be changed to remove the sets with only one or two miRNAs
      point.size = 1, 
      text.scale = 1, 
      line.size = 0.6, 
      mb.ratio = c(0.4, 0.6), 
      mainbar.y.label = "Intersection Size", 
      sets.x.label = "Set Size", 
      cutoff = 2,
      set_size.numbers_size = 10,
      sets.bar.color = c21,
      keep.order = T,
      sets = names(c21))
dev.off()
        
#********************************************************************************************#
# This can be performed to determine how much miRNA we are considering in this test: 
upset2_df <- upset_df
upset2_df <- upset2_df[rowSums(upset_df) > 0,]
# The number of rows of upset2 is the number of miRNAs. 
#********************************************************************************************#
    
############################################################################################################################################################################################
# Tissue specificity index 
miRNA_tau_values <- data.frame(row.names = row.names(normalizedEM))
miRNA_tau_values["tau"] <- NA

# The input of this function (x) is a vector of the normalized expression values of each tissue of one miRNA. The output is the Tau value of this miRNA.  
  fTau <- function(x)
  {
    if(all(!is.na(x)))
    {
      if(min(x, na.rm=TRUE) >= 0)
      {
        if(max(x)!=0)
        {
          x <- (1-(x/max(x)))
          res <- sum(x, na.rm=TRUE)
          res <- res/(length(x)-1)
        } else {
          res <- 0
        }
      } else {
        res <- NA
        print("Expression values have to be positive!")
      } 
    } else {
      res <- NA
      print("No data for this gene avalable.")
    } 
    return(res)
  }
  
Tissue_average2 <- upset_df
for(i in unique(info$Tissue)) { 
  samples <- rownames(info[info["Tissue"] == i,])
  rmeans <- rowMeans(normalizedEM[ ,samples])
  Tissue_average2[i] <- rmeans
 }
  
for (i in rownames(Tissue_average2)) {
  p <- as.vector(log(Tissue_average2[i, ]+1))
  miRNA_tau_values[i, "tau"] <- fTau(p) 
 }

ggplot(miRNA_tau_values, aes(x=tau)) + geom_histogram()
miRNA_tau_values["miRNAs"] <- rownames(miRNA_tau_values)
t_specific <- as.data.frame(miRNA_tau_values[miRNA_tau_values$tau >= 0.90,]) # miRNAs with A tau value equal or greater than 0.90 are considered as tissue specific miRNAs
t_specific <- na.omit(t_specific)
t_specific_ave <- subset(Tissue_average2, rownames(Tissue_average2) %in% rownames(t_specific))
t_specific["top_tissue"] <- NA
t_specific["top-value"] <- max.col(t_specific_ave)
t_specific["top_tissue"] <- colnames(Tissue_average2)[t_specific$`top-value`]
t_specific_norm <-subset(normalizedEM, rownames(normalizedEM) %in% rownames(t_specific))

while (!is.null(dev.list())) dev.off()
pdf(file = "plots/tissue_specific_heatmap.pdf", width = 7, height = 7)
pheatmap(log(t_specific_norm + 1), 
         annotation = info["Tissue"],
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         scale = "row",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_colors = plot_colors)
dev.off()
  
