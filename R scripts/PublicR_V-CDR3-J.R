#!/usr/bin/env Rscript
##############################################################################
##                           SetUp                                   ##
##############################################################################
#Title: shared CDR3
#Description: 
#Installing and loading required packages
if(!require(stringr))
{install.packages("stringr")
}
library(stringr)
if(!require(ggplot2))
{install.packages("ggplot2")
}
library(ggplot2)

if(!require(dplyr))
{install.packages("dplyr")
}
library(dplyr)
#################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#READ files in the directory:
file_names_txt <- as.data.frame(list.files(getwd(), pattern = ".txt"))
colnames(file_names_txt) <- "files"
#WT
file_names_WT <- as.data.frame(file_names_txt[grep(pattern = "WT", file_names_txt$files),])
colnames(file_names_WT) <- "files"
#IL2
file_names_IL2 <- as.data.frame(file_names_txt[grep(pattern = "IL2", file_names_txt$files),])
colnames(file_names_IL2) <- "files"
#IL15
file_names_IL15 <- as.data.frame(file_names_txt[grep(pattern = "IL15", file_names_txt$files),])
colnames(file_names_IL15) <- "files"
rm(file_names_txt)
##############################################################################
##                           Initial Setup                                 ##
##############################################################################
### WT
m=0
for(i in file_names_WT$files){
  m=m+1
  file_tmp <- read.table(i, header = TRUE)
  file_tmp <-as.data.frame(file_tmp[,c(1,3:11)])
  for (j in 2:10){file_tmp[,j] <- as.character(file_tmp[,j])}
  file_tmp$V_CDR3_J <- as.character(paste( file_tmp[,4], file_tmp[,3],file_tmp[,6],sep = "_"))
  assign(paste("WT_", m, sep = ""), file_tmp)
}
rm(i, file_tmp, m)

### IL2
m=0
for(i in file_names_IL2$files){
  m=m+1
  file_tmp <- read.table(i, header = TRUE)
  file_tmp <-as.data.frame(file_tmp[,c(1,3:11)])
  for (j in 2:10){file_tmp[,j] <- as.character(file_tmp[,j])}
  file_tmp$V_CDR3_J <- as.character(paste( file_tmp[,4], file_tmp[,3],file_tmp[,6],sep = "_"))
  assign(paste("IL2_", m, "_KO", sep = ""), file_tmp)
}
rm(i, file_tmp, m)

### IL15
m=0
for(i in file_names_IL15$files){
  m=m+1
  file_tmp <- read.table(i, header = TRUE)
  file_tmp <-as.data.frame(file_tmp[,c(1,3:11)])
  for (j in 2:10){file_tmp[,j] <- as.character(file_tmp[,j])}
  file_tmp$V_CDR3_J <- as.character(paste( file_tmp[,4], file_tmp[,3],file_tmp[,6],sep = "_"))
  assign(paste("IL15_", m, "_KO", sep = ""), file_tmp)
}
rm(i, file_tmp, m)

###############################################################################
##                          Process Data                                     ##
###############################################################################

###########
#  WT     #
###########
###Compare
#WT_shared
tmp_1 <- intersect(WT_1$V_CDR3_J, WT_2$V_CDR3_J)
tmp_1 <- intersect(tmp_1, WT_2$V_CDR3_J)
tmp_1 <- intersect(tmp_1, WT_3$V_CDR3_J)
tmp_1 <- intersect(tmp_1, WT_4$V_CDR3_J)
WT_shared <- as.data.frame(tmp_1)
colnames(WT_shared) <- "CDR3"
#Do table shared
indx <- match(WT_shared$CDR3, WT_1$V_CDR3_J,  nomatch = 0) 
WT_shared <- WT_1[indx,c(4,3,6,11)]
#WT_shared <- WT_shared %>% distinct() ## eliminate duplicated rows

### Add data 
#WT_1
#aggregate:
WT_1 <- WT_1[,c(1,11)]
WT_1 <- aggregate(count~V_CDR3_J,data=WT_1,FUN=sum)
#add counts
indx <- match(WT_shared$V_CDR3_J, WT_1$V_CDR3_J, nomatch = 0) 
WT_shared$WT_1 <- WT_1[indx,2]  
#WT_2
#aggregate:
WT_2 <- WT_2[,c(1,11)]
WT_2 <- aggregate(count~V_CDR3_J,data=WT_2,FUN=sum)
#add counts
indx <- match(WT_shared$V_CDR3_J, WT_2$V_CDR3_J, nomatch = 0) 
WT_shared$WT_2 <- WT_2[indx,2]  
#WT_3
#aggregate:
WT_3 <- WT_3[,c(1,11)]
WT_3 <- aggregate(count~V_CDR3_J,data=WT_3,FUN=sum)
#add counts
indx <- match(WT_shared$V_CDR3_J, WT_3$V_CDR3_J, nomatch = 0) 
WT_shared$WT_3 <- WT_3[indx,2]  
#WT_4
#aggregate:
WT_4 <- WT_4[,c(1,11)]
WT_4 <- aggregate(count~V_CDR3_J,data=WT_4,FUN=sum)
#add counts
indx <- match(WT_shared$V_CDR3_J, WT_4$V_CDR3_J, nomatch = 0) 
WT_shared$WT_4 <- WT_4[indx,2]  
#Percentages
WT_shared$mean_count_WT <- (WT_shared[,5] + WT_shared[,6] + WT_shared[,7] + WT_shared[,8])/4
WT_shared$Percentages_WT <- (WT_shared$mean_count/sum(WT_shared$mean_count))*100

###########
# IL_2KO  #
###########
###Compare
#IL2KO_shared
tmp_1 <- intersect(IL2_1_KO$V_CDR3_J, IL2_2_KO$V_CDR3_J)
tmp_1 <- intersect(tmp_1, IL2_2_KO$V_CDR3_J)
tmp_1 <- intersect(tmp_1, IL2_3_KO$V_CDR3_J)
tmp_1 <- intersect(tmp_1, IL2_4_KO$V_CDR3_J)
IL2KO_shared <- as.data.frame(tmp_1)
colnames(IL2KO_shared) <- "CDR3"
#Do table shared
indx <- match(IL2KO_shared$CDR3, IL2_1_KO$V_CDR3_J,  nomatch = 0) 
IL2KO_shared <- IL2_1_KO[indx,c(4,3,6,11)]
#IL2KO_shared <- IL2KO_shared %>% distinct() ## eliminate duplicated rows

### Add data 
#IL2_1_KO
#aggregate:
IL2_1_KO <- IL2_1_KO[,c(1,11)]
IL2_1_KO <- aggregate(count~V_CDR3_J,data=IL2_1_KO,FUN=sum)
#add counts
indx <- match(IL2KO_shared$V_CDR3_J, IL2_1_KO$V_CDR3_J, nomatch = 0) 
IL2KO_shared$IL2_1_KO <- IL2_1_KO[indx,2]  
#IL2_2_KO
#aggregate:
IL2_2_KO <- IL2_2_KO[,c(1,11)]
IL2_2_KO <- aggregate(count~V_CDR3_J,data=IL2_2_KO,FUN=sum)
#add counts
indx <- match(IL2KO_shared$V_CDR3_J, IL2_2_KO$V_CDR3_J, nomatch = 0) 
IL2KO_shared$IL2_2_KO <- IL2_2_KO[indx,2]  
#IL2_3_KO
#aggregate:
IL2_3_KO <- IL2_3_KO[,c(1,11)]
IL2_3_KO <- aggregate(count~V_CDR3_J,data=IL2_3_KO,FUN=sum)
#add counts
indx <- match(IL2KO_shared$V_CDR3_J, IL2_3_KO$V_CDR3_J, nomatch = 0) 
IL2KO_shared$IL2_3_KO <- IL2_3_KO[indx,2]  
#IL2_4_KO
#aggregate:
IL2_4_KO <- IL2_4_KO[,c(1,11)]
IL2_4_KO <- aggregate(count~V_CDR3_J,data=IL2_4_KO,FUN=sum)
#add counts
indx <- match(IL2KO_shared$V_CDR3_J, IL2_4_KO$V_CDR3_J, nomatch = 0) 
IL2KO_shared$IL2_4_KO <- IL2_4_KO[indx,2]  
#Percentages
IL2KO_shared$mean_count_IL2KO <- (IL2KO_shared[,5] + IL2KO_shared[,6] + IL2KO_shared[,7] + IL2KO_shared[,8])/4
IL2KO_shared$Percentages_IL2KO <- (IL2KO_shared$mean_count/sum(IL2KO_shared$mean_count))*100

###########
# IL_15KO #
###########
###Compare
#IL15KO_shared
tmp_1 <- intersect(IL15_1_KO$V_CDR3_J, IL15_2_KO$V_CDR3_J)
tmp_1 <- intersect(tmp_1, IL15_2_KO$V_CDR3_J)
tmp_1 <- intersect(tmp_1, IL15_3_KO$V_CDR3_J)
tmp_1 <- intersect(tmp_1, IL15_4_KO$V_CDR3_J)
IL15KO_shared <- as.data.frame(tmp_1)
colnames(IL15KO_shared) <- "CDR3"
#Do table shared
indx <- match(IL15KO_shared$CDR3, IL15_1_KO$V_CDR3_J,  nomatch = 0) 
IL15KO_shared <- IL15_1_KO[indx,c(4,3,6,11)]
#IL15KO_shared <- IL15KO_shared %>% distinct() ## eliminate duplicated rows

### Add data 
#IL15_1_KO
#aggregate:
IL15_1_KO <- IL15_1_KO[,c(1,11)]
IL15_1_KO <- aggregate(count~V_CDR3_J,data=IL15_1_KO,FUN=sum)
#add counts
indx <- match(IL15KO_shared$V_CDR3_J, IL15_1_KO$V_CDR3_J, nomatch = 0) 
IL15KO_shared$IL15_1_KO <- IL15_1_KO[indx,2]  
#IL15_2_KO
#aggregate:
IL15_2_KO <- IL15_2_KO[,c(1,11)]
IL15_2_KO <- aggregate(count~V_CDR3_J,data=IL15_2_KO,FUN=sum)
#add counts
indx <- match(IL15KO_shared$V_CDR3_J, IL15_2_KO$V_CDR3_J, nomatch = 0) 
IL15KO_shared$IL15_2_KO <- IL15_2_KO[indx,2]  
#IL15_3_KO
#aggregate:
IL15_3_KO <- IL15_3_KO[,c(1,11)]
IL15_3_KO <- aggregate(count~V_CDR3_J,data=IL15_3_KO,FUN=sum)
#add counts
indx <- match(IL15KO_shared$V_CDR3_J, IL15_3_KO$V_CDR3_J, nomatch = 0) 
IL15KO_shared$IL15_3_KO <- IL15_3_KO[indx,2]  
#IL15_4_KO
#aggregate:
IL15_4_KO <- IL15_4_KO[,c(1,11)]
IL15_4_KO <- aggregate(count~V_CDR3_J,data=IL15_4_KO,FUN=sum)
#add counts
indx <- match(IL15KO_shared$V_CDR3_J, IL15_4_KO$V_CDR3_J, nomatch = 0) 
IL15KO_shared$IL15_4_KO <- IL15_4_KO[indx,2]  
#Percentages
IL15KO_shared$mean_count_IL15KO <- (IL15KO_shared[,5] + IL15KO_shared[,6] + IL15KO_shared[,7] + IL15KO_shared[,8])/4
IL15KO_shared$Percentages_IL15KO <- (IL15KO_shared$mean_count/sum(IL15KO_shared$mean_count))*100

#################################
###       WT vs IL2           ###
#################################
### Shared
tmp_1 <- as.data.frame(intersect(WT_shared$V_CDR3_J, IL2KO_shared$V_CDR3_J))
colnames(tmp_1) <- "V_CDR3_J"
#WT
indx <- match(tmp_1$V_CDR3_J, WT_shared$V_CDR3_J, nomatch = 0) 
WT_IL2_shared <- WT_shared[indx,]
#IL2
indx <- match(WT_IL2_shared$V_CDR3_J, IL2KO_shared$V_CDR3_J, nomatch = 0) 
WT_IL2_shared <- cbind(WT_IL2_shared, IL2KO_shared[indx,5:10])
### Not shared
#WT
tmp_1 <- as.data.frame(setdiff(WT_shared$V_CDR3_J, IL2KO_shared$V_CDR3_J))
colnames(tmp_1) <- "V_CDR3_J"
indx <- match(tmp_1$V_CDR3_J, WT_shared$V_CDR3_J, nomatch = 0) 
WT_not_shared_comp_IL2 <- WT_shared[indx,]
#IL2
tmp_1 <- as.data.frame(setdiff(IL2KO_shared$V_CDR3_J, WT_shared$V_CDR3_J))
colnames(tmp_1) <- "V_CDR3_J"
indx <- match(tmp_1$V_CDR3_J, IL2KO_shared$V_CDR3_J, nomatch = 0) 
IL2_not_shared_comp_WT <- IL2KO_shared[indx,]

#################################
###       WT vs IL15           ###
#################################
### Shared
tmp_1 <- as.data.frame(intersect(WT_shared$V_CDR3_J, IL15KO_shared$V_CDR3_J))
colnames(tmp_1) <- "V_CDR3_J"
#WT
indx <- match(tmp_1$V_CDR3_J, WT_shared$V_CDR3_J, nomatch = 0) 
WT_IL15_shared <- WT_shared[indx,]
#IL15
indx <- match(WT_IL15_shared$V_CDR3_J, IL15KO_shared$V_CDR3_J, nomatch = 0) 
WT_IL15_shared <- cbind(WT_IL15_shared, IL15KO_shared[indx,5:10])
### Not shared
#WT
tmp_1 <- as.data.frame(setdiff(WT_shared$V_CDR3_J, IL15KO_shared$V_CDR3_J))
colnames(tmp_1) <- "V_CDR3_J"
indx <- match(tmp_1$V_CDR3_J, WT_shared$V_CDR3_J, nomatch = 0) 
WT_not_shared_comp_IL15 <- WT_shared[indx,]
#IL15
tmp_1 <- as.data.frame(setdiff(IL15KO_shared$V_CDR3_J, WT_shared$V_CDR3_J))
colnames(tmp_1) <- "V_CDR3_J"
indx <- match(tmp_1$V_CDR3_J, IL15KO_shared$V_CDR3_J, nomatch = 0) 
IL15_not_shared_comp_WT <- IL15KO_shared[indx,]


#################################
###       IL2KO vs IL15           ###
#################################
### Shared
tmp_1 <- as.data.frame(intersect(IL2KO_shared$V_CDR3_J, IL15KO_shared$V_CDR3_J))
colnames(tmp_1) <- "V_CDR3_J"
#IL2KO
indx <- match(tmp_1$V_CDR3_J, IL2KO_shared$V_CDR3_J, nomatch = 0) 
IL2KO_IL15_shared <- IL2KO_shared[indx,]
#IL15
indx <- match(IL2KO_IL15_shared$V_CDR3_J, IL15KO_shared$V_CDR3_J, nomatch = 0) 
IL2KO_IL15_shared <- cbind(IL2KO_IL15_shared, IL15KO_shared[indx,5:10])
### Not shared
#IL2KO
tmp_1 <- as.data.frame(setdiff(IL2KO_shared$V_CDR3_J, IL15KO_shared$V_CDR3_J))
colnames(tmp_1) <- "V_CDR3_J"
indx <- match(tmp_1$V_CDR3_J, IL2KO_shared$V_CDR3_J, nomatch = 0) 
IL2KO_not_shared_comp_IL15 <- IL2KO_shared[indx,]
#IL15
tmp_1 <- as.data.frame(setdiff(IL15KO_shared$V_CDR3_J, IL2KO_shared$V_CDR3_J))
colnames(tmp_1) <- "V_CDR3_J"
indx <- match(tmp_1$V_CDR3_J, IL15KO_shared$V_CDR3_J, nomatch = 0) 
IL15_not_shared_comp_IL2KO <- IL15KO_shared[indx,]

#####################
# Tables to graph
#####################
#### WT vs IL2  #########
WT_vs_IL2_table <- WT_IL2_shared
WT_vs_IL2_table$comparation <- rep("WT_IL2_shared", nrow(WT_IL2_shared))
tmp_1 <- WT_not_shared_comp_IL2
tmp_2 <- as.data.frame(matrix(data = runif(nrow(tmp_1), min=0.00025, max=0.001), nrow = nrow(tmp_1), ncol = 4))
tmp_2 <- cbind(tmp_1, tmp_2, tmp_2[,4], tmp_2[,4], rep("IL2_dep"))
colnames(tmp_2) <- colnames(WT_vs_IL2_table)  
tmp_3 <- IL2_not_shared_comp_WT
tmp_4 <- as.data.frame(matrix(data = runif(nrow(tmp_3), min=0.00025, max=0.001), nrow = nrow(tmp_3), ncol = 4))
tmp_4 <- cbind(tmp_3[,1:4], tmp_4, tmp_4[,4],tmp_4[,4] , tmp_3[,5:10], rep("IL2_indep"))
colnames(tmp_4) <- colnames(WT_vs_IL2_table)  
WT_vs_IL2_table <- rbind(WT_vs_IL2_table, tmp_2, tmp_4)

#### WT vs IL15  #########
WT_vs_IL15_table <- WT_IL15_shared
WT_vs_IL15_table$comparation <- rep("WT_IL15_shared", nrow(WT_IL15_shared))
tmp_1 <- WT_not_shared_comp_IL15
tmp_2 <- as.data.frame(matrix(data = runif(nrow(tmp_1), min=0.00025, max=0.001), nrow = nrow(tmp_1), ncol = 4))
tmp_2 <- cbind(tmp_1, tmp_2, tmp_2[,4], tmp_2[,4], rep("IL15_dep"))
colnames(tmp_2) <- colnames(WT_vs_IL15_table)  
tmp_3 <- IL15_not_shared_comp_WT
tmp_4 <- as.data.frame(matrix(data = runif(nrow(tmp_3), min=0.00025, max=0.001), nrow = nrow(tmp_3), ncol = 4))
tmp_4 <- cbind(tmp_3[,1:4], tmp_4, tmp_4[,4],tmp_4[,4] , tmp_3[,5:10], rep("IL15_indep"))
colnames(tmp_4) <- colnames(WT_vs_IL15_table)  
WT_vs_IL15_table <- rbind(WT_vs_IL15_table, tmp_2, tmp_4)

#### IL2 vs IL15  #########
IL2_vs_IL15_table <- IL2KO_IL15_shared
IL2_vs_IL15_table$comparation <- rep("IL2_IL15_shared", nrow(IL2KO_IL15_shared))
tmp_1 <- IL2KO_not_shared_comp_IL15
tmp_2 <- as.data.frame(matrix(data = runif(nrow(tmp_1), min=0.00025, max=0.001), nrow = nrow(tmp_1), ncol = 4))
tmp_2 <- cbind(tmp_1, tmp_2, tmp_2[,4], tmp_2[,4], rep("IL15_dep"))
colnames(tmp_2) <- colnames(IL2_vs_IL15_table)  
tmp_3 <- IL15_not_shared_comp_IL2KO
tmp_4 <- as.data.frame(matrix(data = runif(nrow(tmp_3), min=0.00025, max=0.001), nrow = nrow(tmp_3), ncol = 4))
tmp_4 <- cbind(tmp_3[,1:4], tmp_4, tmp_4[,4],tmp_4[,4] , tmp_3[,5:10], rep("IL15_indep"))
colnames(tmp_4) <- colnames(IL2_vs_IL15_table)  
IL2_vs_IL15_table <- rbind(IL2_vs_IL15_table, tmp_2, tmp_4)


###############################
#   pValues                   #
###############################

################
#   WT_vs_IL2  #
################
tmp_1 <- WT_vs_IL2_table[which(WT_vs_IL2_table$comparation == "WT_IL2_shared"), ]
tmp_pvalue <- as.data.frame(matrix(data=NA, ncol=1, nrow=nrow(tmp_1)))
rownames(tmp_1) <- NULL
for (i in 1:nrow(tmp_1)){
  cond1 <- as.numeric(tmp_1[i,5:8])/c(sum(WT_shared[,5]), sum(WT_shared[,6]), sum(WT_shared[,7]), sum(WT_shared[,8]))
  cond2 <- as.numeric(tmp_1[i,11:14])/c(sum(IL2KO_shared[,5]), sum(IL2KO_shared[,6]), sum(IL2KO_shared[,7]), sum(IL2KO_shared[,8]))
  if (mean(cond1) > mean(cond2)){
    pvalue <- wilcox.test(cond1, cond2, alternative = "g", paired=FALSE)}else{
    pvalue <- wilcox.test(cond1, cond2, alternative = "l", paired=FALSE) }
    tmp_pvalue[i,1] <- pvalue$p.value
    }
tmp_1$pValue <- tmp_pvalue$V1
#
tmp_2 <- WT_vs_IL2_table[which(WT_vs_IL2_table$comparation == "IL2_dep"), ]
tmp_pvalue <- as.data.frame(matrix(data=NA, ncol=1, nrow=nrow(tmp_2)))
rownames(tmp_2) <- NULL
for (i in 1:nrow(tmp_2)){
  cond1 <- as.numeric(tmp_2[i,5:8])/c(sum(WT_shared[,5]), sum(WT_shared[,6]), sum(WT_shared[,7]), sum(WT_shared[,8]))
  #cond2 <- as.numeric(tmp_2[i,11:14])/c(sum(WT_shared[,5]), sum(WT_shared[,6]), sum(WT_shared[,7]), sum(WT_shared[,8]))
  cond2 <- c(0,0,0,0)
  if (mean(cond1) > mean(cond2)){
    pvalue <- wilcox.test(cond1, cond2, alternative = "g", paired=FALSE)}else{
      pvalue <- wilcox.test(cond1, cond2, alternative = "l", paired=FALSE) }
  tmp_pvalue[i,1] <- round(pvalue$p.value,2)
}
tmp_2$pValue <- tmp_pvalue$V1
#
tmp_3 <- WT_vs_IL2_table[which(WT_vs_IL2_table$comparation == "IL2_indep"), ]
tmp_pvalue <- as.data.frame(matrix(data=NA, ncol=1, nrow=nrow(tmp_3)))
rownames(tmp_3) <- NULL
for (i in 1:nrow(tmp_3)){
  cond1 <- as.numeric(tmp_3[i,5:8])/c(sum(WT_shared[,5]), sum(WT_shared[,6]), sum(WT_shared[,7]), sum(WT_shared[,8]))
  cond2 <- as.numeric(tmp_3[i,11:14])/c(sum(IL2KO_shared[,5]), sum(IL2KO_shared[,6]), sum(IL2KO_shared[,7]), sum(IL2KO_shared[,8]))
  if (mean(cond1) > mean(cond2)){
    pvalue <- wilcox.test(cond1, cond2, alternative = "g", paired=FALSE)}else{
      pvalue <- wilcox.test(cond1, cond2, alternative = "l", paired=FALSE) }
  tmp_pvalue[i,1] <- pvalue$p.value
}
tmp_3$pValue <- tmp_pvalue$V1
#
WT_vs_IL2_table <- rbind(tmp_1, tmp_2, tmp_3)



###################
# Table to export #
###################
tmp_1[,5:8] <- (tmp_1[,5:8]/c(sum(WT_shared[,5]), sum(WT_shared[,6]), sum(WT_shared[,7]), sum(WT_shared[,8])))*100            
tmp_1[,10] <- rowSums(tmp_1[,5:8])/4
tmp_1[,11:14] <- (tmp_1[,11:14]/c(sum(IL2KO_shared[,5]), sum(IL2KO_shared[,6]), sum(IL2KO_shared[,7]), sum(IL2KO_shared[,8])))*100            
tmp_1[,16] <- rowSums(tmp_1[,5:8])/4
tmp_1 <- tmp_1[,c(1:8, 10:14, 16:18)]
tmp_2[,5:8] <- (tmp_2[,5:8]/c(sum(WT_shared[,5]), sum(WT_shared[,6]), sum(WT_shared[,7]), sum(WT_shared[,8])))*100            
tmp_2[,10] <- rowSums(tmp_2[,5:8])/4
tmp_2[,11:14] <- 0
tmp_2[,16] <- 0
tmp_2 <- tmp_2[,c(1:8, 10:14, 16:18)]
tmp_3[,11:14] <- (tmp_3[,11:14]/c(sum(IL2KO_shared[,5]), sum(IL2KO_shared[,6]), sum(IL2KO_shared[,7]), sum(IL2KO_shared[,8])))*100            
tmp_3[,16] <- rowSums(tmp_3[,11:14])/4
tmp_3[,5:8] <- 0
tmp_3[,10] <- 0
tmp_3 <- tmp_3[,c(1:8, 10:14, 16:18)]
WT_vs_IL2_table_to_export <- rbind(tmp_1, tmp_2, tmp_3)
   
   
################
#   WT_vs_IL15  #
################
tmp_1 <- WT_vs_IL15_table[which(WT_vs_IL15_table$comparation == "WT_IL15_shared"), ]
tmp_pvalue <- as.data.frame(matrix(data=NA, ncol=1, nrow=nrow(tmp_1)))
rownames(tmp_1) <- NULL
for (i in 1:nrow(tmp_1)){
  cond1 <- as.numeric(tmp_1[i,5:8])/c(sum(WT_shared[,5]), sum(WT_shared[,6]), sum(WT_shared[,7]), sum(WT_shared[,8]))
  cond2 <- as.numeric(tmp_1[i,11:14])/c(sum(IL15KO_shared[,5]), sum(IL15KO_shared[,6]), sum(IL15KO_shared[,7]), sum(IL15KO_shared[,8]))
  if (mean(cond1) > mean(cond2)){
    pvalue <- wilcox.test(cond1, cond2, alternative = "g", paired=FALSE)}else{
      pvalue <- wilcox.test(cond1, cond2, alternative = "l", paired=FALSE) }
  tmp_pvalue[i,1] <- pvalue$p.value
}
tmp_1$pValue <- tmp_pvalue$V1
#
tmp_2 <- WT_vs_IL15_table[which(WT_vs_IL15_table$comparation == "IL15_dep"), ]
tmp_pvalue <- as.data.frame(matrix(data=NA, ncol=1, nrow=nrow(tmp_2)))
rownames(tmp_2) <- NULL
for (i in 1:nrow(tmp_2)){
  cond1 <- as.numeric(tmp_2[i,5:8])/c(sum(WT_shared[,5]), sum(WT_shared[,6]), sum(WT_shared[,7]), sum(WT_shared[,8]))
  cond2 <- as.numeric(tmp_2[i,11:14])/c(sum(WT_shared[,5]), sum(WT_shared[,6]), sum(WT_shared[,7]), sum(WT_shared[,8]))
  if (mean(cond1) > mean(cond2)){
    pvalue <- wilcox.test(cond1, cond2, alternative = "g", paired=FALSE)}else{
      pvalue <- wilcox.test(cond1, cond2, alternative = "l", paired=FALSE) }
  tmp_pvalue[i,1] <- pvalue$p.value
}
tmp_2$pValue <- tmp_pvalue$V1
#
tmp_3 <- WT_vs_IL15_table[which(WT_vs_IL15_table$comparation == "IL15_indep"), ]
tmp_pvalue <- as.data.frame(matrix(data=NA, ncol=1, nrow=nrow(tmp_3)))
rownames(tmp_3) <- NULL
for (i in 1:nrow(tmp_3)){
  cond1 <- as.numeric(tmp_3[i,5:8])/c(sum(WT_shared[,5]), sum(WT_shared[,6]), sum(WT_shared[,7]), sum(WT_shared[,8]))
  cond2 <- as.numeric(tmp_3[i,11:14])/c(sum(IL15KO_shared[,5]), sum(IL15KO_shared[,6]), sum(IL15KO_shared[,7]), sum(IL15KO_shared[,8]))
  if (mean(cond1) > mean(cond2)){
    pvalue <- wilcox.test(cond1, cond2, alternative = "g", paired=FALSE)}else{
      pvalue <- wilcox.test(cond1, cond2, alternative = "l", paired=FALSE) }
  tmp_pvalue[i,1] <- pvalue$p.value
}
tmp_3$pValue <- tmp_pvalue$V1
#
WT_vs_IL15_table <- rbind(tmp_1, tmp_2, tmp_3)

###################
# Table to export #
###################
tmp_1[,5:8] <- (tmp_1[,5:8]/c(sum(WT_shared[,5]), sum(WT_shared[,6]), sum(WT_shared[,7]), sum(WT_shared[,8])))*100            
tmp_1[,10] <- rowSums(tmp_1[,5:8])/4
tmp_1[,11:14] <- (tmp_1[,11:14]/c(sum(IL15KO_shared[,5]), sum(IL15KO_shared[,6]), sum(IL15KO_shared[,7]), sum(IL15KO_shared[,8])))*100            
tmp_1[,16] <- rowSums(tmp_1[,5:8])/4
tmp_1 <- tmp_1[,c(1:8, 10:14, 16:18)]
tmp_2[,5:8] <- (tmp_2[,5:8]/c(sum(WT_shared[,5]), sum(WT_shared[,6]), sum(WT_shared[,7]), sum(WT_shared[,8])))*100            
tmp_2[,10] <- rowSums(tmp_2[,5:8])/4
tmp_2[,11:14] <- 0
tmp_2[,16] <- 0
tmp_2 <- tmp_2[,c(1:8, 10:14, 16:18)]
tmp_3[,11:14] <- (tmp_3[,11:14]/c(sum(IL15KO_shared[,5]), sum(IL15KO_shared[,6]), sum(IL15KO_shared[,7]), sum(IL15KO_shared[,8])))*100            
tmp_3[,16] <- rowSums(tmp_3[,11:14])/4
tmp_3[,5:8] <- 0
tmp_3[,10] <- 0
tmp_3 <- tmp_3[,c(1:8, 10:14, 16:18)]
WT_vs_IL15_table_to_export <- rbind(tmp_1, tmp_2, tmp_3)
################
#   IL2_vs_IL15  #
################
tmp_1 <- IL2_vs_IL15_table[which(IL2_vs_IL15_table$comparation == "IL2_IL15_shared"), ]
tmp_pvalue <- as.data.frame(matrix(data=NA, ncol=1, nrow=nrow(tmp_1)))
rownames(tmp_1) <- NULL
for (i in 1:nrow(tmp_1)){
  cond1 <- as.numeric(tmp_1[i,5:8])/c(sum(IL2KO_shared[,5]), sum(IL2KO_shared[,6]), sum(IL2KO_shared[,7]), sum(IL2KO_shared[,8]))
  cond2 <- as.numeric(tmp_1[i,11:14])/c(sum(IL15KO_shared[,5]), sum(IL15KO_shared[,6]), sum(IL15KO_shared[,7]), sum(IL15KO_shared[,8]))
  if (mean(cond1) > mean(cond2)){
    pvalue <- wilcox.test(cond1, cond2, alternative = "g", paired=FALSE)}else{
      pvalue <- wilcox.test(cond1, cond2, alternative = "l", paired=FALSE) }
  tmp_pvalue[i,1] <- pvalue$p.value
}
tmp_1$pValue <- tmp_pvalue$V1
#
tmp_2 <- IL2_vs_IL15_table[which(IL2_vs_IL15_table$comparation == "IL15_dep"), ]
tmp_pvalue <- as.data.frame(matrix(data=NA, ncol=1, nrow=nrow(tmp_2)))
rownames(tmp_2) <- NULL
for (i in 1:nrow(tmp_2)){
  cond1 <- as.numeric(tmp_2[i,5:8])/c(sum(IL2KO_shared[,5]), sum(IL2KO_shared[,6]), sum(IL2KO_shared[,7]), sum(IL2KO_shared[,8]))
  cond2 <- as.numeric(tmp_2[i,11:14])/c(sum(IL2KO_shared[,5]), sum(IL2KO_shared[,6]), sum(IL2KO_shared[,7]), sum(IL2KO_shared[,8]))
  if (mean(cond1) > mean(cond2)){
    pvalue <- wilcox.test(cond1, cond2, alternative = "g", paired=FALSE)}else{
      pvalue <- wilcox.test(cond1, cond2, alternative = "l", paired=FALSE) }
  tmp_pvalue[i,1] <- pvalue$p.value
}
tmp_2$pValue <- tmp_pvalue$V1
#
tmp_3 <- IL2_vs_IL15_table[which(IL2_vs_IL15_table$comparation == "IL15_indep"), ]
tmp_pvalue <- as.data.frame(matrix(data=NA, ncol=1, nrow=nrow(tmp_3)))
rownames(tmp_3) <- NULL
for (i in 1:nrow(tmp_3)){
  cond1 <- as.numeric(tmp_3[i,5:8])/c(sum(IL2KO_shared[,5]), sum(IL2KO_shared[,6]), sum(IL2KO_shared[,7]), sum(IL2KO_shared[,8]))
  cond2 <- as.numeric(tmp_3[i,11:14])/c(sum(IL15KO_shared[,5]), sum(IL15KO_shared[,6]), sum(IL15KO_shared[,7]), sum(IL15KO_shared[,8]))
  if (mean(cond1) > mean(cond2)){
    pvalue <- wilcox.test(cond1, cond2, alternative = "g", paired=FALSE)}else{
      pvalue <- wilcox.test(cond1, cond2, alternative = "l", paired=FALSE) }
  tmp_pvalue[i,1] <- pvalue$p.value
}
tmp_3$pValue <- tmp_pvalue$V1
#
IL2_vs_IL15_table <- rbind(tmp_1, tmp_2, tmp_3)

###################
# Table to export #
###################
tmp_1[,5:8] <- (tmp_1[,5:8]/c(sum(IL2KO_shared[,5]), sum(IL2KO_shared[,6]), sum(IL2KO_shared[,7]), sum(IL2KO_shared[,8])))*100            
tmp_1[,10] <- rowSums(tmp_1[,5:8])/4
tmp_1[,11:14] <- (tmp_1[,11:14]/c(sum(IL15KO_shared[,5]), sum(IL15KO_shared[,6]), sum(IL15KO_shared[,7]), sum(IL15KO_shared[,8])))*100            
tmp_1[,16] <- rowSums(tmp_1[,5:8])/4
tmp_1 <- tmp_1[,c(1:8, 10:14, 16:18)]
tmp_2[,5:8] <- (tmp_2[,5:8]/c(sum(IL2KO_shared[,5]), sum(IL2KO_shared[,6]), sum(IL2KO_shared[,7]), sum(IL2KO_shared[,8])))*100            
tmp_2[,10] <- rowSums(tmp_2[,5:8])/4
tmp_2[,11:14] <- 0
tmp_2[,16] <- 0
tmp_2 <- tmp_2[,c(1:8, 10:14, 16:18)]
tmp_3[,11:14] <- (tmp_3[,11:14]/c(sum(IL15KO_shared[,5]), sum(IL15KO_shared[,6]), sum(IL15KO_shared[,7]), sum(IL15KO_shared[,8])))*100            
tmp_3[,16] <- rowSums(tmp_3[,11:14])/4
tmp_3[,5:8] <- 0
tmp_3[,10] <- 0
tmp_3 <- tmp_3[,c(1:8, 10:14, 16:18)]
IL2KO_vs_IL15_table_to_export <- rbind(tmp_1, tmp_2, tmp_3)

###############################################################################
##                                   Graphics                                ##
###############################################################################

##############
#   WT vs IL2
##############
png(file = "WT_vs_IL2-KO_PERCENTAGES_SHARED.png", width = 15, height = 15, units = "cm", res = 330)
par(cex.main=1.2, cex.axis=1, mar=c(5, 5, 5, 5))
plot(log10(WT_vs_IL2_table$Percentages_WT), log10(WT_vs_IL2_table$Percentages_IL2KO), xlab="wt-frequency (%)", ylab="IL2o-frequency (%)", 
     main = "WT - IL2 KO", xlim=c(-3.5, 0.7), ylim=c(-3.5, 0.7), cex = 0.7, pch=16, yaxt='n', xaxt='n')
axis(side= 1, at= c(0,-1,-2,-3.3), labels=c("10^0", "10^-1", "10^-2", "n.d."))
axis(side= 2, at= c(0,-1,-2,-3.3), labels=c("10^0", "10^-1", "10^-2", "n.d."), las=2)
abline(h=-2.9, col="grey", lwd=2, lty= 3)
abline(v=-2.9, col="grey", lwd=2, lty= 3)
text(paste("n=", as.character(nrow(WT_IL2_shared))), x=-1.5 , y=0.5)
text(paste("n=", as.character(nrow(IL2_not_shared_comp_WT))), x=-3.3 , y=0.5)
text(paste("n=", as.character(nrow(WT_not_shared_comp_IL2))), x=0.3 , y=-3.25)
lines(log10( WT_vs_IL2_table[which(WT_vs_IL2_table$pValue < 0.05),10]), log10(WT_vs_IL2_table[which(WT_vs_IL2_table$pValue < 0.05),16]),type="p", col="red", cex = 0.7, pch=16)
dev.off()

##############
#   WT vs IL15
##############
png(file = "WT_vs_IL15-KO_PERCENTAGES_SHARED.png", width = 15, height = 15, units = "cm",res = 330)
par(cex.main=1.2, cex.axis=1, mar=c(5, 5, 5, 5))
plot(log10(WT_vs_IL15_table$Percentages_WT), log10(WT_vs_IL15_table$Percentages_IL15KO), xlab="wt-frequency (%)", ylab="IL15o-frequency (%)", 
     main = "WT - IL15 KO", xlim=c(-3.5, 0.7), ylim=c(-3.5, 0.7), cex = 0.7, pch=16, yaxt='n', xaxt='n')
axis(side= 1, at= c(0,-1,-2,-3.3), labels=c("10^0", "10^-1", "10^-2", "n.d."))
axis(side= 2, at= c(0,-1,-2,-3.3), labels=c("10^0", "10^-1", "10^-2", "n.d."), las=2)
abline(h=-2.9, col="grey", lwd=2, lty= 3)
abline(v=-2.9, col="grey", lwd=2, lty= 3)
text(paste("n=", as.character(nrow(WT_IL15_shared))), x=-1.5 , y=0.5)
text(paste("n=", as.character(nrow(IL15_not_shared_comp_WT))), x=-3.3 , y=0.5)
text(paste("n=", as.character(nrow(WT_not_shared_comp_IL15))), x=0.3 , y=-3.25)
lines(log10( WT_vs_IL15_table[which(WT_vs_IL15_table$pValue < 0.05),10]), log10(WT_vs_IL15_table[which(WT_vs_IL15_table$pValue < 0.05),16]),type="p", col="red", cex = 0.7, pch=16)
dev.off()


###############################################################################
##                             Write tables                                  ##
###############################################################################

write.table(WT_IL2_shared, file="IL2_independent.txt", sep ="\t" , quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(WT_not_shared_comp_IL2, file="IL2_dependent.txt", sep ="\t" , quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(WT_IL15_shared, file="IL15_independent.txt", sep ="\t" , quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(WT_not_shared_comp_IL15, file="IL15_dependent.txt", sep ="\t" , quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(IL2_not_shared_comp_WT, file="IL2_exclusive.txt", sep ="\t" , quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(IL15_not_shared_comp_WT, file="IL15_exclusive.txt", sep ="\t" , quote = FALSE, col.names = TRUE, row.names = FALSE)


write.table(WT_vs_IL2_table_to_export, file="WT_vs_IL2_table.txt", sep ="\t" , quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(WT_vs_IL15_table_to_export, file="WT_vs_IL15_table.txt", sep ="\t" , quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(IL2KO_vs_IL15_table_to_export, file="IL2_vs_IL15_table.txt", sep ="\t" , quote = FALSE, col.names = TRUE, row.names = FALSE)

write.table(WT_shared, file="WT_PublicR.txt", sep ="\t" , quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(IL2KO_shared, file="IL2KO_PublicR.txt", sep ="\t" , quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(IL15KO_shared, file="IL15KO_PublicR.txt", sep ="\t" , quote = FALSE, col.names = TRUE, row.names = FALSE)

