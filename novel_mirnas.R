setwd("sheep-miRNAome")

library("seqRFLP")

# Create fasta files of unfiltered miRNAs with >= 5 score 

needed_colums <- c("miRDeep2.score","provisional.id","example.miRBase.miRNA.with.the.same.seed","consensus.mature.sequence","consensus.star.sequence","consensus.precursor.sequence","precursor.coordinate")

mirdeep_results <- read.csv("result_26_07_2021_t_13_36_10.csv", sep="\t", skip=26) # The output of mirdeep2 core command

novel_mirnas <- mirdeep_results[mirdeep_results["miRBase.miRNA"] == "-", ]



filtered_novel_mirnas <- novel_mirnas[ ,needed_colums]

filtered_novel_mirnas$miRDeep2.score <- as.numeric(as.character(filtered_novel_mirnas$miRDeep2.score))

filtered_novel_mirnas <- filtered_novel_mirnas[filtered_novel_mirnas["miRDeep2.score"] >= 5, ]

filtered_novel_mirnas$names <- paste(filtered_novel_mirnas$provisional.id, filtered_novel_mirnas$example.miRBase.miRNA.with.the.same.seed, sep = ".")


filtered_novel_mirnas$names <- sub("^",">", filtered_novel_mirnas$names)


df_mature <- filtered_novel_mirnas[ ,c("names","consensus.mature.sequence")]

df_star <- filtered_novel_mirnas[ ,c("names","consensus.star.sequence")]

df_premiRNA <- filtered_novel_mirnas[ ,c("names","consensus.precursor.sequence")]

D_mature <- do.call(rbind, lapply(seq(nrow(df_mature)), function(i) t(df_mature[i, ])))

write.table(D_mature, file = "mature.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)

D_star <- do.call(rbind, lapply(seq(nrow(df_star)), function(i) t(df_star[i, ])))

write.table(D_star, file = "star.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)

D_premirna <- do.call(rbind, lapply(seq(nrow(df_premiRNA)), function(i) t(df_premiRNA[i, ])))

write.table(D_premirna, file = "premirna.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)



########################################################################################################


filtered_novel_mirnas <- novel_mirnas[ ,needed_colums]

filtered_novel_mirnas$miRDeep2.score <- as.numeric(as.character(filtered_novel_mirnas$miRDeep2.score))

filtered_novel_mirnas <- filtered_novel_mirnas[filtered_novel_mirnas["miRDeep2.score"] >= 5, ]

filtered_novel_mirnas$names <- paste(filtered_novel_mirnas$provisional.id, filtered_novel_mirnas$example.miRBase.miRNA.with.the.same.seed, sep = ".")

rownames(filtered_novel_mirnas) <- filtered_novel_mirnas$names


# Open the filtered miRNAs table 

mirnas <- read.csv("all_mirnas.bed", header = FALSE) 

mirnas <- as.vector(mirnas$V4)

clustered_filtered_novel_mirnas <- filtered_novel_mirnas[mirnas, ]

clustered_filtered_novel_mirnas <- na.omit(clustered_filtered_novel_mirnas)



clustered_filtered_novel_mirnas$example.miRBase.miRNA.with.the.same.seed <- as.character(clustered_filtered_novel_mirnas$example.miRBase.miRNA.with.the.same.seed)

clustered_filtered_novel_mirnas["star_sequences_names"] <- clustered_filtered_novel_mirnas$example.miRBase.miRNA.with.the.same.seed

clustered_filtered_novel_mirnas["premirna_names"] <- clustered_filtered_novel_mirnas$example.miRBase.miRNA.with.the.same.seed

clustered_filtered_novel_mirnas$consensus.precursor.sequence <- as.character(clustered_filtered_novel_mirnas$consensus.precursor.sequence)

clustered_filtered_novel_mirnas$consensus.star.sequence <- as.character(clustered_filtered_novel_mirnas$consensus.star.sequence)

clustered_filtered_novel_mirnas$consensus.mature.sequence <- as.character(clustered_filtered_novel_mirnas$consensus.mature.sequence)

sapply(clustered_filtered_novel_mirnas, class)



###################################################################################################


# Set correct name and suffix for each mature miRNA based on the premiRNA sequence

# -5p or -3p suffixes are excised from the name to name the premiRNA. The opposite suffix is set to the name of star sequence (-3p --> -5p).

for (i in 1:nrow(clustered_filtered_novel_mirnas)) {
  
   if (endsWith(as.character(clustered_filtered_novel_mirnas[i, "example.miRBase.miRNA.with.the.same.seed"]), "-5p")) {
     
     clustered_filtered_novel_mirnas[i, "star_sequences_names"] <- sub("-5p", "-3p", clustered_filtered_novel_mirnas[i, "example.miRBase.miRNA.with.the.same.seed"])
     clustered_filtered_novel_mirnas[i, "premirna_names"] <- sub("-5p", "", clustered_filtered_novel_mirnas[i, "example.miRBase.miRNA.with.the.same.seed"])
     
     } else if (endsWith(as.character(clustered_filtered_novel_mirnas[i, "example.miRBase.miRNA.with.the.same.seed"]), "-3p")) {
       
       clustered_filtered_novel_mirnas[i, "star_sequences_names"] <- sub("-3p", "-5p", clustered_filtered_novel_mirnas[i, "example.miRBase.miRNA.with.the.same.seed"])
       clustered_filtered_novel_mirnas[i, "premirna_names"] <- sub("-3p", "", clustered_filtered_novel_mirnas[i, "example.miRBase.miRNA.with.the.same.seed"])
   }

  
}


# Some miRNAs didn't have suffix in the mature name, so they were given manually.
# split premiRNA sequences and search the mature sequence in each arm.


substr(clustered_filtered_novel_mirnas[3, "consensus.precursor.sequence"], 0, round(nchar(clustered_filtered_novel_mirnas[3, "consensus.precursor.sequence"])/2 ))
  
round(nchar(clustered_filtered_novel_mirnas[3, "consensus.precursor.sequence"])/2)
  

# premiRNAs are split
for (i in 1:nrow(clustered_filtered_novel_mirnas)) {

clustered_filtered_novel_mirnas[i,"5p"] <- substr(clustered_filtered_novel_mirnas[i, "consensus.precursor.sequence"], 0, round(nchar(clustered_filtered_novel_mirnas[i, "consensus.precursor.sequence"])/2 ))
clustered_filtered_novel_mirnas[i,"3p"] <- substr(clustered_filtered_novel_mirnas[i, "consensus.precursor.sequence"], round(nchar(clustered_filtered_novel_mirnas[i, "consensus.precursor.sequence"])/2 ), nchar(clustered_filtered_novel_mirnas[i, "consensus.precursor.sequence"]))

}


clustered_filtered_novel_mirnas["mature_sense"] <- NA
clustered_filtered_novel_mirnas["star_sense"] <- NA


for (i in 1:nrow(clustered_filtered_novel_mirnas)) {
 
   star_length <- nchar(clustered_filtered_novel_mirnas[i, "consensus.star.sequence"])
   
   mature_length <- nchar(clustered_filtered_novel_mirnas[i, "consensus.mature.sequence"])
   
   mature <- substr(clustered_filtered_novel_mirnas[i, "consensus.mature.sequence"], 4, mature_length -4)
   
   star <- substr(clustered_filtered_novel_mirnas[i, "consensus.star.sequence"], 4, star_length -4)
  
  if (grepl(mature, clustered_filtered_novel_mirnas[i, "5p"])) {
    
    clustered_filtered_novel_mirnas[i, "mature_sense"] <- "5p"
 
  } else if (grepl(mature, clustered_filtered_novel_mirnas[i, "3p"])) {
    
    clustered_filtered_novel_mirnas[i, "mature_sense"] <- "3p"
   
  }
   
   
   if (grepl(star, clustered_filtered_novel_mirnas[i, "3p"])) {
  
    clustered_filtered_novel_mirnas[i, "star_sense"] <- "3p"
    
  } else if (grepl(star, clustered_filtered_novel_mirnas[i, "5p"])) {
    
    clustered_filtered_novel_mirnas[i, "star_sense"] <- "5p"
    
  }

}

 
 
##########################################################################################################3


# novel mature and premiRNA are written into two fasta files and are joined withe the mature and premiRNAs of mirbase
   
clustered_filtered_novel_mirnas["premirna_names"][clustered_filtered_novel_mirnas["premirna_names"] == "-"] <- "_"

clustered_filtered_novel_mirnas["prem_I_D"] <- paste(clustered_filtered_novel_mirnas$provisional.id, clustered_filtered_novel_mirnas$premirna_names, sep = ".")

clustered_filtered_novel_mirnas["mat_I_D"] <- paste(clustered_filtered_novel_mirnas$prem_I_D, clustered_filtered_novel_mirnas$mature_sense, sep = "-")

clustered_filtered_novel_mirnas["str_I_D"] <- paste(clustered_filtered_novel_mirnas$prem_I_D, clustered_filtered_novel_mirnas$star_sense, sep = "-")


clustered_filtered_novel_mirnas$prem_I_D <- sub("^",">", clustered_filtered_novel_mirnas$prem_I_D)
clustered_filtered_novel_mirnas$mat_I_D <- sub("^",">", clustered_filtered_novel_mirnas$mat_I_D)
clustered_filtered_novel_mirnas$str_I_D <- sub("^",">", clustered_filtered_novel_mirnas$str_I_D)


df_mature2 <- clustered_filtered_novel_mirnas[ ,c("mat_I_D","consensus.mature.sequence")]

df_star2 <- clustered_filtered_novel_mirnas[ ,c("str_I_D","consensus.star.sequence")]

df_premiRNA2 <- clustered_filtered_novel_mirnas[ ,c("prem_I_D","consensus.precursor.sequence")]

D_mature <- do.call(rbind, lapply(seq(nrow(df_mature2)), function(i) t(df_mature2[i, ])))

write.table(D_mature, file = "Qmature.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)

D_star <- do.call(rbind, lapply(seq(nrow(df_star2)), function(i) t(df_star2[i, ])))

write.table(D_star, file = "Qstar.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)

D_premirna <- do.call(rbind, lapply(seq(nrow(df_premiRNA2)), function(i) t(df_premiRNA2[i, ])))

write.table(D_premirna, file = "Qpremirna.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)


# merge mirbase and novel mature and premiRNAs

system("cat Qmature.fa Qstar.fa > Qsmmirna.fa")

system("cat Qsmmirna.fa sheep_mature_mirna.fa > Qallmirna.fa")

system("cat Qpremirna.fa sheep_premirna.fa > Qallpremirna.fa")
                                    
