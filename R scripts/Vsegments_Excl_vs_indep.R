#!/usr/bin/env Rscript
##############################################################################
##                           SetUP                                          ##
##############################################################################
#Title: IL2-IL15 Treg dependent: V segments
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
if(!require(matrixStats))
{install.packages("matrixStats")
}
library(matrixStats)
#################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#READ file:
IL2_exclusive <- read.table("IL2_exclusive.txt", sep ="\t", header = TRUE, comment="", quote="")
IL15_exclusive <- read.table("IL15_exclusive.txt", sep ="\t", header = TRUE, comment="", quote="")
IL2_independent <- read.table("IL2_independent.txt", sep ="\t", header = TRUE, comment="", quote="")
IL15_independent <- read.table("IL15_independent.txt", sep ="\t", header = TRUE, comment="", quote="")
V_segments <- read.csv(file = "VAsegments.txt", sep = "\t")


IL2_independent <-IL2_independent[,c(1:4,11:16)]
IL15_independent <-IL15_independent[,c(1:4,11:16)]

###############################################################################
##                          Process Data                                     ##
###############################################################################
###### TRAV present in all samples

#Sums WT
sum_counts_WT <- matrix(data=NA, ncol=2, nrow = 4)
sum_counts_WT[,1] <- c("WT1", "WT2", "WT3", "WT4")
j=1
for (i in 5:8) {
sum_counts_WT[j,2] <- c(sum(IL2_exclusive[,i])+sum(IL2_independent[,i]))
j=j+1
}
colnames(sum_counts_WT) <- c("replicate", "counts")
sum_counts_WT <- as.data.frame(sum_counts_WT)
sum_counts_WT[,1] <- as.character(sum_counts_WT[,1])
sum_counts_WT[,2] <- as.numeric(as.character(sum_counts_WT[,2]))
###############
#     IL2     #
###############
#IL2 dep
condition= IL2_exclusive
tmp <- condition
tmp_1 <- as.data.frame(tmp[,c(2,5)])
tmp_1 <- aggregate(IL2_3_KO ~ v, data=tmp_1, sum)
tmp_1[,2] <- (tmp_1[,2]/sum(tmp_1[,2]))*100
tmp_2 <- as.data.frame(tmp[,c(2,6)])
tmp_2 <- aggregate(IL2_5_KO ~ v, data=tmp_2, sum)
tmp_2[,2] <- (tmp_2[,2]/sum(tmp_2[,2]))*100
tmp_3 <- as.data.frame(tmp[,c(2,7)])
tmp_3 <- aggregate(IL2_6_KO ~ v, data=tmp_3, sum)
tmp_3[,2] <- (tmp_3[,2]/sum(tmp_3[,2]))*100
tmp_4 <- as.data.frame(tmp[,c(2,8)])
tmp_4 <- aggregate(IL2_7_KO ~ v, data=tmp_4, sum)
tmp_4[,2] <- (tmp_4[,2]/sum(tmp_4[,2]))*100
Data_IL2_exclusive <- cbind(tmp_1, tmp_2[,2], tmp_3[,2], tmp_4[,2])
colnames(Data_IL2_exclusive) <- c("V","WT1", "WT2", "WT3", "WT4")

#IL2 indep
condition= IL2_independent
tmp <- condition
tmp_1 <- as.data.frame(tmp[,c(2,5)])
tmp_1 <- aggregate(IL2_3_KO ~ v, data=tmp_1, sum)
tmp_1[,2] <- (tmp_1[,2]/sum(tmp_1[,2]))*100
tmp_2 <- as.data.frame(tmp[,c(2,6)])
tmp_2 <- aggregate(IL2_5_KO ~ v, data=tmp_2, sum)
tmp_2[,2] <- (tmp_2[,2]/sum(tmp_2[,2]))*100
tmp_3 <- as.data.frame(tmp[,c(2,7)])
tmp_3 <- aggregate(IL2_6_KO ~ v, data=tmp_3, sum)
tmp_3[,2] <- (tmp_3[,2]/sum(tmp_3[,2]))*100
tmp_4 <- as.data.frame(tmp[,c(2,8)])
tmp_4 <- aggregate(IL2_7_KO ~ v, data=tmp_4, sum)
tmp_4[,2] <- (tmp_4[,2]/sum(tmp_4[,2]))*100
Data_IL2_independent <- cbind(tmp_1, tmp_2[,2], tmp_3[,2], tmp_4[,2])
colnames(Data_IL2_independent) <- c("V","WT1", "WT2", "WT3", "WT4")

###############
#     IL15     #
###############
#IL15 dep
condition= IL15_exclusive
tmp <- condition
tmp_1 <- as.data.frame(tmp[,c(2,5)])
tmp_1 <- aggregate(IL15_9_KO ~ v, data=tmp_1, sum)
tmp_1[,2] <- (tmp_1[,2]/sum(tmp_1[,2]))*100
tmp_2 <- as.data.frame(tmp[,c(2,6)])
tmp_2 <- aggregate(IL15_10_KO ~ v, data=tmp_2, sum)
tmp_2[,2] <- (tmp_2[,2]/sum(tmp_2[,2]))*100
tmp_3 <- as.data.frame(tmp[,c(2,7)])
tmp_3 <- aggregate(IL15_11_KO ~ v, data=tmp_3, sum)
tmp_3[,2] <- (tmp_3[,2]/sum(tmp_3[,2]))*100
tmp_4 <- as.data.frame(tmp[,c(2,8)])
tmp_4 <- aggregate(IL15_12_KO ~ v, data=tmp_4, sum)
tmp_4[,2] <- (tmp_4[,2]/sum(tmp_4[,2]))*100
Data_IL15_exclusive <- cbind(tmp_1, tmp_2[,2], tmp_3[,2], tmp_4[,2])
colnames(Data_IL15_exclusive) <- c("V","WT1", "WT2", "WT3", "WT4")

#IL15 indep
condition= IL15_independent
tmp <- condition
tmp_1 <- as.data.frame(tmp[,c(2,5)])
tmp_1 <- aggregate(IL15_9_KO ~ v, data=tmp_1, sum)
tmp_1[,2] <- (tmp_1[,2]/sum(tmp_1[,2]))*100
tmp_2 <- as.data.frame(tmp[,c(2,6)])
tmp_2 <- aggregate(IL15_10_KO ~ v, data=tmp_2, sum)
tmp_2[,2] <- (tmp_2[,2]/sum(tmp_2[,2]))*100
tmp_3 <- as.data.frame(tmp[,c(2,7)])
tmp_3 <- aggregate(IL15_11_KO ~ v, data=tmp_3, sum)
tmp_3[,2] <- (tmp_3[,2]/sum(tmp_3[,2]))*100
tmp_4 <- as.data.frame(tmp[,c(2,8)])
tmp_4 <- aggregate(IL15_12_KO ~ v, data=tmp_4, sum)
tmp_4[,2] <- (tmp_4[,2]/sum(tmp_4[,2]))*100
Data_IL15_independent <- cbind(tmp_1, tmp_2[,2], tmp_3[,2], tmp_4[,2])
colnames(Data_IL15_independent) <- c("V","WT1", "WT2", "WT3", "WT4")

###############
# Mean and SD #
###############
Data_IL2_exclusive$mean <- rowMeans(Data_IL2_exclusive[,2:5])
Data_IL2_independent$mean <- rowMeans(Data_IL2_independent[,2:5])
Data_IL15_exclusive$mean <- rowMeans(Data_IL15_exclusive[,2:5])
Data_IL15_independent$mean <- rowMeans(Data_IL15_independent[,2:5])
#
Data_IL2_exclusive$SD <- rowSds(as.matrix(Data_IL2_exclusive[,2:5]))
Data_IL2_independent$SD <- rowSds(as.matrix(Data_IL2_independent[,2:5]))
Data_IL15_exclusive$SD <- rowSds(as.matrix(Data_IL15_exclusive[,2:5]))
Data_IL15_independent$SD <- rowSds(as.matrix(Data_IL15_independent[,2:5]))




################################
# Complete the tables and order#
################################
#######
# IL2 #
#######
tmp_1 <- setdiff(Data_IL2_exclusive[,1], Data_IL2_independent[,1]) 
tmp_2 <- matrix(data = 0, ncol = ncol(Data_IL2_exclusive) , nrow= length(tmp_1))
tmp_2[,1] <- tmp_1
colnames(tmp_2) <- colnames(Data_IL2_exclusive)
Data_IL2_independent <- rbind(Data_IL2_independent, tmp_2)
Data_IL2_independent[,1] <- as.character(Data_IL2_independent[,1])
Data_IL2_independent <- Data_IL2_independent[order(Data_IL2_independent[,1]), ]
#
tmp_1 <- setdiff(Data_IL2_independent[,1], Data_IL2_exclusive[,1]) 
tmp_2 <- matrix(data = 0, ncol = ncol(Data_IL2_independent) , nrow= length(tmp_1))
tmp_2[,1] <- tmp_1
colnames(tmp_2) <- colnames(Data_IL2_independent)
Data_IL2_exclusive <- rbind(Data_IL2_exclusive, tmp_2)
Data_IL2_exclusive[,1] <- as.character(Data_IL2_exclusive[,1])
Data_IL2_exclusive <- Data_IL2_exclusive[order(Data_IL2_exclusive[,1]), ]
for(i in 2:7){Data_IL2_exclusive[,i] <- as.numeric(Data_IL2_exclusive[,i])}
for(i in 2:7){Data_IL2_independent[,i] <- as.numeric(Data_IL2_independent[,i])}
#Order by Independent
Data_IL2_independent <- Data_IL2_independent[order(Data_IL2_independent[,6], decreasing = TRUE ), ]
rownames(Data_IL2_exclusive) <- NULL
rownames(Data_IL2_independent) <- NULL
index <- match(Data_IL2_independent[,1], Data_IL2_exclusive[,1], nomatch = 0)
Data_IL2_exclusive <- Data_IL2_exclusive[index,]
rownames(Data_IL2_exclusive) <- NULL



########
# IL15 #
########
tmp_1 <- setdiff(Data_IL15_exclusive[,1], Data_IL15_independent[,1]) 
tmp_2 <- matrix(data = 0, ncol = ncol(Data_IL15_exclusive) , nrow= length(tmp_1))
tmp_2[,1] <- tmp_1
colnames(tmp_2) <- colnames(Data_IL15_exclusive)
Data_IL15_independent <- rbind(Data_IL15_independent, tmp_2)
Data_IL15_independent[,1] <- as.character(Data_IL15_independent[,1])
Data_IL15_independent <- Data_IL15_independent[order(Data_IL15_independent[,1]), ]
#
tmp_1 <- setdiff(Data_IL15_independent[,1], Data_IL15_exclusive[,1]) 
tmp_2 <- matrix(data = 0, ncol = ncol(Data_IL15_independent) , nrow= length(tmp_1))
tmp_2[,1] <- tmp_1
colnames(tmp_2) <- colnames(Data_IL15_independent)
Data_IL15_exclusive <- rbind(Data_IL15_exclusive, tmp_2)
Data_IL15_exclusive[,1] <- as.character(Data_IL15_exclusive[,1])
Data_IL15_exclusive <- Data_IL15_exclusive[order(Data_IL15_exclusive[,1]), ]
for(i in 2:7){Data_IL15_exclusive[,i] <- as.numeric(Data_IL15_exclusive[,i])}
for(i in 2:7){Data_IL15_independent[,i] <- as.numeric(Data_IL15_independent[,i])}
#Order by Independent
Data_IL15_independent <- Data_IL15_independent[order(Data_IL15_independent[,6], decreasing = TRUE ), ]
rownames(Data_IL15_exclusive) <- NULL
rownames(Data_IL15_independent) <- NULL
index <- match(Data_IL15_independent[,1], Data_IL15_exclusive[,1], nomatch = 0)
Data_IL15_exclusive <- Data_IL15_exclusive[index,]
rownames(Data_IL15_exclusive) <- NULL



########################################
# Transform data frame Graphics tables #
########################################

######
#IL2 #
######
#IL2_exclusive
IL2_exclusive_graphic <- matrix(data = NA, ncol = 6, nrow = nrow(Data_IL2_exclusive))
colnames(IL2_exclusive_graphic) <- c("x_values", "mean", "SD", "V_segment", "condition", "pValue")
IL2_exclusive_graphic[,2] <- c(Data_IL2_exclusive[,6])
IL2_exclusive_graphic[,3] <- Data_IL2_exclusive[,7]
IL2_exclusive_graphic[,4] <- Data_IL2_exclusive[,1]
IL2_exclusive_graphic[,5] <- rep("IL-2 Exclusive",nrow(Data_IL2_exclusive))
#Add X
IL2_exclusive_graphic[,1] <- seq(from=1, to=nrow(Data_IL2_exclusive)*3, by=3)
#pValues
for (i in 1:nrow(Data_IL2_exclusive)) {
  if (as.numeric(Data_IL2_exclusive[i,6]) > as.numeric(Data_IL2_independent[i,6]))
  {pValue <- wilcox.test(as.numeric(Data_IL2_exclusive[i,2:5]),as.numeric(Data_IL2_independent[i,2:5]),  alternative = "g", paired = FALSE)
  }else {pValue <- wilcox.test(as.numeric(Data_IL2_exclusive[i,2:5]),as.numeric(Data_IL2_independent[i,2:5]),  alternative = "l", paired = FALSE)}
  IL2_exclusive_graphic[i,6] <- round(pValue$p.value, digits = 4)}

#IL2_independent
IL2_independent_graphic <- matrix(data = NA, ncol = 6, nrow = nrow(Data_IL2_independent))
colnames(IL2_independent_graphic) <- c("x_values", "mean", "SD", "V_segment", "condition", "pValue")
IL2_independent_graphic[,2] <- c(Data_IL2_independent[,6])
IL2_independent_graphic[,3] <- Data_IL2_independent[,7]
IL2_independent_graphic[,4] <- Data_IL2_independent[,1]
IL2_independent_graphic[,5] <- rep("IL-2 independent",nrow(Data_IL2_independent))
#Add X
IL2_independent_graphic[,1] <- seq(from=2, to=nrow(Data_IL2_independent)*3, by=3)
# Together
IL2_graphic <- as.data.frame(rbind(IL2_independent_graphic, IL2_exclusive_graphic))
for(i in c(1:3,6)){IL2_graphic[,i] <- as.numeric(as.character(IL2_graphic[,i]))}
for(i in c(4,5)){IL2_graphic[,i] <- as.character(IL2_graphic[,i])}


######
#IL15 #
######
#IL15_exclusive
IL15_exclusive_graphic <- matrix(data = NA, ncol = 6, nrow = nrow(Data_IL15_exclusive))
colnames(IL15_exclusive_graphic) <- c("x_values", "mean", "SD", "V_segment", "condition", "pValue")
IL15_exclusive_graphic[,2] <- c(Data_IL15_exclusive[,6])
IL15_exclusive_graphic[,3] <- Data_IL15_exclusive[,7]
IL15_exclusive_graphic[,4] <- Data_IL15_exclusive[,1]
IL15_exclusive_graphic[,5] <- rep("IL-15 Exclusive",nrow(Data_IL15_exclusive))
#Add X
IL15_exclusive_graphic[,1] <- seq(from=1, to=nrow(Data_IL15_exclusive)*3, by=3)
#pValues
for (i in 1:nrow(Data_IL15_exclusive)) {
  if (as.numeric(Data_IL15_exclusive[i,6]) > as.numeric(Data_IL15_independent[i,6]))
  {pValue <- wilcox.test(as.numeric(Data_IL15_exclusive[i,2:5]),as.numeric(Data_IL15_independent[i,2:5]),  alternative = "g", paired = FALSE)
  }else {pValue <- wilcox.test(as.numeric(Data_IL15_exclusive[i,2:5]),as.numeric(Data_IL15_independent[i,2:5]),  alternative = "l", paired = FALSE)}
  IL15_exclusive_graphic[i,6] <- round(pValue$p.value, digits = 4)}

#IL15_independent
IL15_independent_graphic <- matrix(data = NA, ncol = 6, nrow = nrow(Data_IL15_independent))
colnames(IL15_independent_graphic) <- c("x_values", "mean", "SD", "V_segment", "condition", "pValue")
IL15_independent_graphic[,2] <- c(Data_IL15_independent[,6])
IL15_independent_graphic[,3] <- Data_IL15_independent[,7]
IL15_independent_graphic[,4] <- Data_IL15_independent[,1]
IL15_independent_graphic[,5] <- rep("IL-2 independent",nrow(Data_IL15_independent))
#Add X
IL15_independent_graphic[,1] <- seq(from=2, to=nrow(Data_IL15_independent)*3, by=3)
# Together
IL15_graphic <- as.data.frame(rbind(IL15_independent_graphic, IL15_exclusive_graphic))
for(i in c(1:3,6)){IL15_graphic[,i] <- as.numeric(as.character(IL15_graphic[,i]))}
for(i in c(4,5)){IL15_graphic[,i] <- as.character(IL15_graphic[,i])}

###############################################################################
##                          Select specific TRAVs                            ##
###############################################################################

###########
#   IL2   #
###########
index_1 <- which(IL2_graphic[,2]>=5) 
tmp_trav <- IL2_graphic[index_1,]
tmp_trav <- unique(tmp_trav[,4])
tmp <- matrix(data = NA, ncol = ncol(IL2_graphic), nrow = 0)
colnames(tmp) <- colnames(IL2_graphic)
for (i in 1:length(tmp_trav)) {
  tmp <- rbind(tmp, IL2_graphic[which(IL2_graphic[,4] == tmp_trav[i]), ])}
IL2_graphic_selected <- tmp
rownames(IL2_graphic_selected) <- NULL
tmp <- c(seq(from=1, to=nrow(IL2_graphic_selected )*1.5, by=3), seq(from=2, to=nrow(IL2_graphic_selected )*1.5, by=3))
IL2_graphic_selected[,1] <- sort(tmp)

###########
#   IL15   #
###########
index_1 <- which(IL15_graphic[,2]>=5) 
tmp_trav <- IL15_graphic[index_1,]
tmp_trav <- unique(tmp_trav[,4])
tmp <- matrix(data = NA, ncol = ncol(IL15_graphic), nrow = 0)
colnames(tmp) <- colnames(IL15_graphic)
for (i in 1:length(tmp_trav)) {
  tmp <- rbind(tmp, IL15_graphic[which(IL15_graphic[,4] == tmp_trav[i]), ])}
IL15_graphic_selected <- tmp
rownames(IL15_graphic_selected) <- NULL
tmp <- c(seq(from=1, to=nrow(IL15_graphic_selected )*1.5, by=3), seq(from=2, to=nrow(IL15_graphic_selected )*1.5, by=3))
IL15_graphic_selected[,1] <- sort(tmp)


##############################
##           Stars          ##
##############################
# *P < 0.05, **P < 0.01, ***P < 0.001, ****P < 0.0001 
stars_signif_IL2 <- matrix(data = NA, ncol = 5, nrow = (nrow(IL2_graphic_selected)/2))
colnames(stars_signif_IL2) <- c("x_pos", "y_pos", "pValue", "stars", "y_segment")
stars_signif_IL2[,1] <- IL2_graphic_selected[seq(from=1, to=nrow(IL2_graphic_selected), by=2 ) ,1]+0.5
tmp <- as.data.frame(cbind( IL2_graphic_selected[,2], IL2_graphic_selected[,3]))
tmp$sum <- rowSums(tmp)
tmp <- cbind(tmp[seq(from=1, to=nrow(IL2_graphic_selected), by=2 ),3],tmp[seq(from=2, to=nrow(IL2_graphic_selected), by=2 ),3])
for (i in 1:nrow(tmp)) {stars_signif_IL2[i,2] <- max(tmp[i,1],tmp[i,2])}
tmp <- stars_signif_IL2[,2]
stars_signif_IL2[,2] <- tmp+(round(max(tmp) /10,0)) 
stars_signif_IL2[,5] <- tmp+(round(max(tmp) /15,0)) 
stars_signif_IL2[,3] <- IL2_graphic_selected[seq(from=2, to=nrow(IL2_graphic_selected)+1, by=2 ),6]     #1:(nrow(IL2_graphic_selected)/2)

if (length (stars_signif_IL2[which(stars_signif_IL2[,3] > 0.05),]) == 5){
  tmp_ns <- as.data.frame(t(stars_signif_IL2[which(stars_signif_IL2[,3] > 0.05),]))
  tmp_ns[1:nrow(tmp_ns),4] <- "ns"}
if (length (stars_signif_IL2[which(stars_signif_IL2[,3] > 0.05),]) > 5){
  tmp_ns <- as.data.frame((stars_signif_IL2[which(stars_signif_IL2[,3] > 0.05),]))
  tmp_ns[1:nrow(tmp_ns),4] <- "ns"}
stars_signif_IL2 <- stars_signif_IL2[which(stars_signif_IL2[,3] <= 0.05),]
rownames(stars_signif_IL2) <- NULL
stars_signif_IL2[which(stars_signif_IL2[,3] >= 0.01),4] <- "*"
stars_signif_IL2[which(stars_signif_IL2[,3] <= 0.01 &  stars_signif_IL2[,3] >= 0.001), 4] <- "**"
stars_signif_IL2[which(stars_signif_IL2[,3] <= 0.001 &  stars_signif_IL2[,3] >= 0.0001), 4] <- "***"
#stars_signif_IL2[ which(stars_signif_IL2[,3] <= 0.0001), 4] <- "****"
if(exists("tmp_ns")) {stars_signif_IL2 <- rbind(stars_signif_IL2, tmp_ns)
for (i in c(1:3,5)) {stars_signif_IL2[,i] <- as.numeric(stars_signif_IL2[,i])}
stars_signif_IL2 <- stars_signif_IL2[order(stars_signif_IL2[,1]),]}
rm(tmp_ns)



# *P < 0.05, **P < 0.01, ***P < 0.001, ****P < 0.0001 
stars_signif_IL15 <- matrix(data = NA, ncol = 5, nrow = (nrow(IL15_graphic_selected)/2))
colnames(stars_signif_IL15) <- c("x_pos", "y_pos", "pValue", "stars", "y_segment")
stars_signif_IL15[,1] <- IL15_graphic_selected[seq(from=1, to=nrow(IL15_graphic_selected), by=2 ) ,1]+0.5
tmp <- as.data.frame(cbind( IL15_graphic_selected[,2], IL15_graphic_selected[,3]))
tmp$sum <- rowSums(tmp)
tmp <- cbind(tmp[seq(from=1, to=nrow(IL15_graphic_selected), by=2 ),3],tmp[seq(from=2, to=nrow(IL15_graphic_selected), by=2 ),3])
for (i in 1:nrow(tmp)) {stars_signif_IL15[i,2] <- max(tmp[i,1],tmp[i,2])}
tmp <- stars_signif_IL15[,2]
stars_signif_IL15[,2] <- tmp+(round(max(tmp) /10,0)) 
stars_signif_IL15[,5] <- tmp+(round(max(tmp) /15,0)) 
stars_signif_IL15[,3] <- IL15_graphic_selected[seq(from=2, to=nrow(IL15_graphic_selected)+1, by=2 ),6]     #1:(nrow(IL15_graphic_selected)/2)
if (length (stars_signif_IL15[which(stars_signif_IL15[,3] > 0.05),]) == 5){
  tmp_ns <- as.data.frame(t(stars_signif_IL15[which(stars_signif_IL15[,3] > 0.05),]))
  tmp_ns[1:nrow(tmp_ns),4] <- "ns"}
if (length (stars_signif_IL15[which(stars_signif_IL15[,3] > 0.05),]) > 5){
  tmp_ns <- as.data.frame((stars_signif_IL15[which(stars_signif_IL15[,3] > 0.05),]))
  tmp_ns[1:nrow(tmp_ns),4] <- "ns"}
stars_signif_IL15 <- stars_signif_IL15[which(stars_signif_IL15[,3] <= 0.05),]
rownames(stars_signif_IL15) <- NULL
stars_signif_IL15[which(stars_signif_IL15[,3] >= 0.01),4] <- "*"
stars_signif_IL15[which(stars_signif_IL15[,3] <= 0.01 &  stars_signif_IL15[,3] >= 0.001), 4] <- "**"
stars_signif_IL15[which(stars_signif_IL15[,3] <= 0.001 &  stars_signif_IL15[,3] >= 0.0001), 4] <- "***"
#stars_signif_IL15[ which(stars_signif_IL15[,3] <= 0.0001), 4] <- "****"
if(exists("tmp_ns")) {stars_signif_IL15 <- rbind(stars_signif_IL15, tmp_ns)
for (i in c(1:3,5)) {stars_signif_IL15[,i] <- as.numeric(stars_signif_IL15[,i])}
stars_signif_IL15 <- stars_signif_IL15[order(stars_signif_IL15[,1]),]}




###################
# Eliminate "TRAV"#
###################
IL2_graphic_selected[,4] <- gsub("TRAV", "", IL2_graphic_selected[,4])
IL15_graphic_selected[,4] <- gsub("TRAV", "", IL15_graphic_selected[,4])

#####################################################################################
#                               Graphics                                            #
#####################################################################################

size_chr=16
custom_theme <- function () { theme(legend.text=element_text(size=size_chr), plot.title= element_text(size=size_chr, hjust = 0.5), 
                                    axis.text.x = element_text(size=size_chr, colour = "black",angle = 90, hjust=0.95,vjust=0.5), axis.text.y = element_text(size=size_chr, colour = "black"), 
                                    axis.title.x= element_text(size=size_chr+2, hjust = 0.5) , axis.title.y= element_text(size=size_chr+4),
                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.ticks.length=unit(.25, "cm"), axis.ticks = element_line(size = 1.5),
                                    legend.position = "none", axis.line = element_line(size = 1, colour = "black"))}  





ggplot(IL2_graphic_selected, aes(x = x_values, y = mean, fill=condition))+ 
  geom_bar(stat="identity")+
  xlab("") +
  ylab("") +
  scale_y_continuous(expand = c(0, 0.3)) +
  custom_theme()+
  guides(fill = guide_legend(reverse = TRUE))+
  expand_limits(y=c(0, 20))+ 
  scale_fill_manual(values=c("#156cff", "#ff492e"))+
  scale_x_continuous(expand = c(0, 1), breaks=seq(from=1.5, to=IL2_graphic_selected[(nrow(IL2_graphic_selected)),1], by=3), labels = unique(IL2_graphic_selected[,4]))+
  geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.2,position=position_dodge(.9),color="black")+
  annotate("text", x=as.numeric(stars_signif_IL2[,1]), y=as.numeric(stars_signif_IL2[,2]), label=stars_signif_IL2[,4], color="black", size=6)+
  ggsave("Vseg_IL2_excl_indep.pdf", width=12, height=10, units = "cm")


ggplot(IL15_graphic_selected, aes(x = x_values, y = mean, fill=condition))+ 
  geom_bar(stat="identity")+
  xlab("") +
  ylab("") +
  scale_y_continuous(expand = c(0, 0.3)) +
  custom_theme()+
  guides(fill = guide_legend(reverse = TRUE))+
  expand_limits(y=c(0, 25))+ 
  scale_fill_manual(values=c("#156cff", "#ff492e"))+
  scale_x_continuous(expand = c(0, 1), breaks=seq(from=1.5, to=IL15_graphic_selected[(nrow(IL15_graphic_selected)),1], by=3), labels = unique(IL15_graphic_selected[,4]))+
  geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.2,position=position_dodge(.9),color="black")+
  annotate("text", x=as.numeric(stars_signif_IL15[,1]), y=as.numeric(stars_signif_IL15[,2]), label=stars_signif_IL15[,4], color="black", size=6)+
  ggsave("Vseg_IL15_excl_indep.pdf", width=14, height=10, units = "cm")

