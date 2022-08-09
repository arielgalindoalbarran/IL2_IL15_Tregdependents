#!/usr/bin/env Rscript
##############################################################################
##                           SetUP                                          ##
##############################################################################
#Title: IL2-IL15 Treg dependent: Spectratype
#Description: 
#Installing and loading required packages
if(!require(ggplot2))
{install.packages("ggplot2")
}
if(!require(reshape))
{install.packages("reshape")
}
if(!require(RColorBrewer))
{install.packages("RColorBrewer")
}

library(ggplot2)
library(reshape)
library(RColorBrewer)
library(matrixStats)
library(dplyr)
##############################################################################
##                           Initial Setup                                 ##
##############################################################################
#set working directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#####################################
#Fill these:
IL2_dependent <- read.table("IL2_dependent.txt", header = TRUE)
IL2_independent <- read.table("IL2_independent.txt", header = TRUE)
IL15_dependent <- read.table("IL15_dependent.txt", header = TRUE)
IL15_independent <- read.table("IL15_independent.txt", header = TRUE)
for (i in 1:4){IL2_dependent[,i] <- as.character(IL2_dependent[,i])}
for (i in 1:4){IL2_independent[,i] <- as.character(IL2_independent[,i])}
for (i in 1:4){IL15_dependent[,i] <- as.character(IL15_dependent[,i])}
for (i in 1:4){IL15_independent[,i] <- as.character(IL15_independent[,i])}


###############################################################################
##                       Process Data  IL2                                  ##
###############################################################################
Data_IL2 <- matrix(data=NA, nrow=7, ncol = 14)
Data_IL2[,1] <- c(10,11,12,13,14,15,16)
colnames(Data_IL2) <- c("Len", rep("IL2_dep", 4),"mean_IL2_dep", "SD_IL2_dep",rep("IL2_indep", 4), "mean_IL2_indep", "SD_IL2_indep", "pValue_IL2_dep_indep")

j=2
for (i in 5:8){
tmp_1 <- nchar(as.character(IL2_dependent$cdr3aa))
tmp_2 <- IL2_dependent[,i]
tmp_3 <- as.data.frame(cbind(tmp_1, tmp_2))
colnames(tmp_3) <- c("Len", "counts")
tmp_4 <- aggregate(counts ~ Len, data=tmp_3, sum)
tmp_4$counts <- (tmp_4$counts/sum(tmp_4$counts))*100
indx <- match(tmp_4[,1] , Data_IL2[,1], nomatch = 0) 
tmp_4 <- tmp_4[indx,]
Data_IL2[,j] <- tmp_4$counts
j=j+1
}
j=8
for (i in 5:8){
  tmp_1 <- nchar(as.character(IL2_independent$cdr3aa))
  tmp_2 <- IL2_independent[,i]
  tmp_3 <- as.data.frame(cbind(tmp_1, tmp_2))
  colnames(tmp_3) <- c("Len", "counts")
  tmp_4 <- aggregate(counts ~ Len, data=tmp_3, sum)
  tmp_4$counts <- (tmp_4$counts/sum(tmp_4$counts))*100
  indx <- match(tmp_4[,1] , Data_IL2[,1], nomatch = 0) 
  tmp_4 <- tmp_4[indx,]
  Data_IL2[,j] <- tmp_4$counts
  j=j+1
}

#means
Data_IL2[,6] <- rowMeans(Data_IL2[,2:5])
Data_IL2[,12] <- rowMeans(Data_IL2[,8:11])
#SD
Data_IL2[,7] <- rowSds(as.matrix(Data_IL2[,2:5]))
Data_IL2[,13] <- rowSds(as.matrix(Data_IL2[,8:11]))

#pValue
for (i in 1:7) {
  if (Data_IL2[i,6] > Data_IL2[i,12]){pValue <- wilcox.test(as.vector(Data_IL2[i,2:5]),as.vector(Data_IL2[i,8:11]),  alternative = "g", paired = FALSE)
    }else {pValue <- wilcox.test(as.vector(Data_IL2[i,2:5]),as.vector(Data_IL2[i,8:11]),  alternative = "l", paired = FALSE)}
    Data_IL2[i,14] <- round(pValue$p.value, digits = 4)}



###############################################################################
##                       Process Data  IL15                                  ##
###############################################################################
Data_IL15 <- matrix(data=NA, nrow=7, ncol = 14)
Data_IL15[,1] <- c(10,11,12,13,14,15,16)
colnames(Data_IL15) <- c("Len", rep("IL15_dep", 4),"mean_IL15_dep", "SD_IL15_dep",rep("IL15_indep", 4), "mean_IL15_indep", "SD_IL15_indep", "pValue_IL15_dep_indep")

j=2
for (i in 5:8){
  tmp_1 <- nchar(as.character(IL15_dependent$cdr3aa))
  tmp_2 <- IL15_dependent[,i]
  tmp_3 <- as.data.frame(cbind(tmp_1, tmp_2))
  colnames(tmp_3) <- c("Len", "counts")
  tmp_4 <- aggregate(counts ~ Len, data=tmp_3, sum)
  tmp_4$counts <- (tmp_4$counts/sum(tmp_4$counts))*100
  indx <- match(tmp_4[,1] , Data_IL15[,1], nomatch = 0) 
  tmp_4 <- tmp_4[indx,]
  Data_IL15[,j] <- tmp_4$counts
  j=j+1
}
j=8
for (i in 5:8){
  tmp_1 <- nchar(as.character(IL15_independent$cdr3aa))
  tmp_2 <- IL15_independent[,i]
  tmp_3 <- as.data.frame(cbind(tmp_1, tmp_2))
  colnames(tmp_3) <- c("Len", "counts")
  tmp_4 <- aggregate(counts ~ Len, data=tmp_3, sum)
  tmp_4$counts <- (tmp_4$counts/sum(tmp_4$counts))*100
  indx <- match(tmp_4[,1] , Data_IL15[,1], nomatch = 0) 
  tmp_4 <- tmp_4[indx,]
  Data_IL15[,j] <- tmp_4$counts
  j=j+1
}
#means
Data_IL15[,6] <- rowMeans(Data_IL15[,2:5])
Data_IL15[,12] <- rowMeans(Data_IL15[,8:11])
#SD
Data_IL15[,7] <- rowSds(as.matrix(Data_IL15[,2:5]))
Data_IL15[,13] <- rowSds(as.matrix(Data_IL15[,8:11]))

#pValue
for (i in 1:7) {
  if (Data_IL15[i,6] > Data_IL15[i,12]){pValue <- wilcox.test(as.vector(Data_IL15[i,2:5]),as.vector(Data_IL15[i,8:11]),  alternative = "g", paired = FALSE)
  }else {pValue <- wilcox.test(as.vector(Data_IL15[i,2:5]),as.vector(Data_IL15[i,8:11]),  alternative = "l", paired = FALSE)}
  Data_IL15[i,14] <- round(pValue$p.value, digits = 4)}


##########################################
#Process data to graph by strains
##########################################
#IL2
tmp_1 <- Data_IL2[,c(1,6,7,14)]
tmp_1 <- as.data.frame(tmp_1)
tmp_1$strain <- rep("IL2_dependent", nrow(tmp_1))
colnames(tmp_1) <- c("Len", "mean", "SD", "pValue", "strain")
tmp_2 <- Data_IL2[,c(1,12,13,14)]
tmp_2 <- as.data.frame(tmp_2)
tmp_2$strain <- rep("IL2_independent", nrow(tmp_2))
colnames(tmp_2) <- c("Len", "mean", "SD", "pValue", "strain")
# join the samples
Data_IL2_graph <- rbind(tmp_1, tmp_2)
#change X values
Data_IL2_graph[,1] <- c(seq(from=1, to=21, by=3), seq(from=1, to=21, by=3)+1)

#Select stars 
# *P < 0.05, **P < 0.01, ***P < 0.001, ****P < 0.0001 
stars_signif_IL2 <- matrix(data = NA, ncol = 5, nrow = (nrow(Data_IL2_graph)/2))
colnames(stars_signif_IL2) <- c("x_pos", "y_pos", "pValue", "stars", "y_segment")
stars_signif_IL2[,1] <- Data_IL2_graph[1:(nrow(Data_IL2_graph)/2) ,1]+0.5
tmp <- as.data.frame(cbind( Data_IL2_graph[,2], Data_IL2_graph[,3]))
tmp$sum <- rowSums(tmp)
tmp <- cbind(tmp[1:(nrow(tmp)/2),3],tmp[((nrow(tmp)/2)+1):nrow(tmp),3])
for (i in 1:nrow(tmp)) {stars_signif_IL2[i,2] <- max(tmp[i,1],tmp[i,2])}
tmp <- stars_signif_IL2[,2]
stars_signif_IL2[,2] <- tmp+(round(max(tmp) /10,0)) 
stars_signif_IL2[,5] <- tmp+(round(max(tmp) /15,0)) 
stars_signif_IL2[,3] <- Data_IL2_graph[1:(nrow(Data_IL2_graph)/2) ,4]
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

#IL15
tmp_1 <- Data_IL15[,c(1,6,7,14)]
tmp_1 <- as.data.frame(tmp_1)
tmp_1$strain <- rep("IL15_dependent", nrow(tmp_1))
colnames(tmp_1) <- c("Len", "mean", "SD", "pValue", "strain")
tmp_2 <- Data_IL15[,c(1,12,13,14)]
tmp_2 <- as.data.frame(tmp_2)
tmp_2$strain <- rep("IL15_independent", nrow(tmp_2))
colnames(tmp_2) <- c("Len", "mean", "SD", "pValue", "strain")
# join the samples
Data_IL15_graph <- rbind(tmp_1, tmp_2)
#change X values
Data_IL15_graph[,1] <- c(seq(from=1, to=21, by=3), seq(from=1, to=21, by=3)+1)
#Select stars 
# *P < 0.05, **P < 0.01, ***P < 0.001, ****P < 0.0001 
stars_signif_IL15 <- matrix(data = NA, ncol = 5, nrow = (nrow(Data_IL15_graph)/2))
colnames(stars_signif_IL15) <- c("x_pos", "y_pos", "pValue", "stars", "y_segment")
stars_signif_IL15[,1] <- Data_IL15_graph[1:(nrow(Data_IL15_graph)/2) ,1]+0.5
tmp <- as.data.frame(cbind( Data_IL15_graph[,2], Data_IL15_graph[,3]))
tmp$sum <- rowSums(tmp)
tmp <- cbind(tmp[1:(nrow(tmp)/2),3],tmp[((nrow(tmp)/2)+1):nrow(tmp),3])
for (i in 1:nrow(tmp)) {stars_signif_IL15[i,2] <- max(tmp[i,1],tmp[i,2])}
tmp <- stars_signif_IL15[,2]
stars_signif_IL15[,2] <- tmp+(round(max(tmp) /10,0)) 
stars_signif_IL15[,5] <- tmp+(round(max(tmp) /15,0)) 
stars_signif_IL15[,3] <- Data_IL15_graph[1:(nrow(Data_IL15_graph)/2) ,4]
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


#########################################
#      To graph general mean length 
#########################################
#####   IL2   #####
Data_IL2_means <- matrix(data=NA, nrow=2, ncol=5)
colnames(Data_IL2_means) <- c("coordinate", "mean", "SD", "condition", "pValue")
Data_IL2_means[,1] <- c(1,2)
Data_IL2_means[,4] <- c("IL2_dependent", "IL2_independent")
#
tmp_1 <- Data_IL2[,1:5]
tmp_means_1 <- matrix(data=NA, ncol=1, nrow=4)
for (i in 1:4) {tmp_means_1[i,1] <- sum(tmp_1[,1]*tmp_1[,i+1]) / (sum(tmp_1[,1]))} 
Data_IL2_means[1,2]  <- round(mean(tmp_means_1[,1]), 2)
Data_IL2_means[1,3]  <- round(sd(tmp_means_1[,1]), 2)
#
tmp_2 <- Data_IL2[,c(1,8:11)]
tmp_means_2 <- matrix(data=NA, ncol=1, nrow=4)
for (i in 1:4) {tmp_means_2[i,1] <- sum(tmp_2[,1]*tmp_2[,i+1]) / (sum(tmp_2[,1]))} 
Data_IL2_means[2,2]  <- round(mean(tmp_means_2[,1]), 2)
Data_IL2_means[2,3]  <- round(sd(tmp_means_2[,1]), 2)
#pValue
  if (as.numeric(Data_IL2_means[1,2]) > as.numeric(Data_IL2_means[2,2])){pValue <- wilcox.test(tmp_means_1[,1], tmp_means_2[,1],  alternative = "g", paired = FALSE) }else 
      {pValue <- wilcox.test(tmp_means_1[,1], tmp_means_2[,1],  alternative = "l", paired = FALSE)}
pValue <- round(pValue$p.value, digits = 4)
Data_IL2_means[,5] <- c(pValue, pValue)

Data_IL2_means <- as.data.frame(Data_IL2_means)
for(i in c(1:3, 5)){Data_IL2_means[,i] <- as.numeric(as.character(Data_IL2_means[,i]))}
Data_IL2_means[,4] <- as.character(Data_IL2_means[,4])

#####   IL15   #####
Data_IL15_means <- matrix(data=NA, nrow=2, ncol=5)
colnames(Data_IL15_means) <- c("coordinate", "mean", "SD", "condition", "pValue")
Data_IL15_means[,1] <- c(1,2)
Data_IL15_means[,4] <- c("IL15_dependent", "IL15_independent")
#
tmp_1 <- Data_IL15[,1:5]
tmp_means_1 <- matrix(data=NA, ncol=1, nrow=4)
for (i in 1:4) {tmp_means_1[i,1] <- sum(tmp_1[,1]*tmp_1[,i+1]) / (sum(tmp_1[,1]))} 
Data_IL15_means[1,2]  <- round(mean(tmp_means_1[,1]), 2)
Data_IL15_means[1,3]  <- round(sd(tmp_means_1[,1]), 2)
#
tmp_2 <- Data_IL15[,c(1,8:11)]
tmp_means_2 <- matrix(data=NA, ncol=1, nrow=4)
for (i in 1:4) {tmp_means_2[i,1] <- sum(tmp_2[,1]*tmp_2[,i+1]) / (sum(tmp_2[,1]))} 
Data_IL15_means[2,2]  <- round(mean(tmp_means_2[,1]), 2)
Data_IL15_means[2,3]  <- round(sd(tmp_means_2[,1]), 2)
#pValue
if (as.numeric(Data_IL15_means[1,2]) > as.numeric(Data_IL15_means[2,2])){pValue <- wilcox.test(tmp_means_1[,1], tmp_means_2[,1],  alternative = "g", paired = FALSE) }else 
{pValue <- wilcox.test(tmp_means_1[,1], tmp_means_2[,1],  alternative = "l", paired = FALSE)}
pValue <- round(pValue$p.value, digits = 4)
Data_IL15_means[,5] <- c(pValue, pValue)

Data_IL15_means <- as.data.frame(Data_IL15_means)
for(i in c(1:3, 5)){Data_IL15_means[,i] <- as.numeric(as.character(Data_IL15_means[,i]))}
Data_IL15_means[,4] <- as.character(Data_IL15_means[,4])

# *P < 0.05, **P < 0.01, ***P < 0.001, ****P < 0.0001 
stars_signif_means_IL2 <- Data_IL2_means[which(Data_IL2_means[,5] <= 0.05),]
rownames(stars_signif_means_IL2) <- NULL
stars_9 <- stars_signif_means_IL2[which(stars_signif_means_IL2[,5] >= 0.01),]
rownames(stars_9) <- NULL
stars_10 <- stars_signif_means_IL2[which(stars_signif_means_IL2[,5] <= 0.01),]
stars_10 <- stars_10[which(stars_10[,4] >= 0.001),]
rownames(stars_10) <- NULL
stars_11 <- stars_signif_means_IL2[which(stars_signif_means_IL2[,5] <= 0.001),]
stars_11 <- stars_11[which(stars_11[,4] >= 0.0001),]
rownames(stars_11) <- NULL
stars_12 <- stars_signif_means_IL2[which(stars_signif_means_IL2[,5] <= 0.0001),]

# *P < 0.05, **P < 0.01, ***P < 0.001, ****P < 0.0001 
stars_signif_means_IL15 <- Data_IL15_means[which(Data_IL15_means[,5] <= 0.05),]
rownames(stars_signif_means_IL15) <- NULL
stars_13 <- stars_signif_means_IL15[which(stars_signif_means_IL15[,5] >= 0.01),]
rownames(stars_13) <- NULL
stars_14 <- stars_signif_means_IL15[which(stars_signif_means_IL15[,5] <= 0.01),]
stars_14 <- stars_14[which(stars_14[,4] >= 0.001),]
rownames(stars_14) <- NULL
stars_15 <- stars_signif_means_IL15[which(stars_signif_means_IL15[,5] <= 0.001),]
stars_15 <- stars_15[which(stars_15[,4] >= 0.0001),]
rownames(stars_15) <- NULL
stars_16 <- stars_signif_means_IL15[which(stars_signif_means_IL15[,5] <= 0.0001),]


############################################################################## 
#                               Graphics                                     #
##############################################################################
##Get colors from RGB to hexadecimal:
#rgb2hex <- function(r,g,b) sprintf('#%s',paste(as.hexmode(c(r,g,b)),collapse = ''))
#rgb2hex(0,186,56) 
#c("#00ba38", "#f8766d", "#619cff")
size_chr=16
custom_theme <- function () { theme(legend.text=element_text(size=size_chr), plot.title= element_text(size=size_chr, hjust = 0.5), 
                                    axis.text.x = element_text(size=size_chr, colour = "black"), axis.text.y = element_text(size=size_chr, colour = "black"), 
                                    axis.title.x= element_text(size=size_chr+2, hjust = 0.5) , axis.title.y= element_text(size=size_chr+2),
                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.ticks.length=unit(.25, "cm"), axis.ticks = element_line(size = 1.5),
                                    legend.position = "none",
                                    axis.line = element_line(size = 1, colour = "black"))}  

custom_theme_2 <- function () { theme(legend.text=element_text(size=size_chr), plot.title= element_text(size=16, hjust = 0.5), 
                                    axis.text.x = element_text(size=size_chr-4, colour = "black"), axis.text.y = element_text(size=size_chr, colour = "black"), 
                                    axis.title.x= element_text(size=size_chr+2, hjust = 0.5) , axis.title.y= element_text(size=size_chr+2),
                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",
                                    panel.background = element_blank(), axis.ticks.length=unit(.25, "cm"), axis.ticks.x=element_blank(),axis.ticks = element_line(size = 1.5),
                                    axis.line = element_line(size = 1, colour = "black"))}  


ggplot(Data_IL2_graph, aes(x = Len, y = mean, fill = strain))+ 
  geom_bar(width = 0.95, stat = "identity") +
  xlab("") +
  ylab("") +
  labs(fill="", title = "") +
  scale_y_continuous(expand = c(0, 0)) +
  custom_theme() +
  guides(fill = guide_legend(reverse = FALSE))+
  expand_limits(y=c(0, 50))+ 
  scale_fill_manual(values=c("#156cff", "#ff492e"))+
  scale_x_continuous(breaks=seq(from=1.5, to=22.5, by=3), labels = seq(from=10, to=17))+
  geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=0.3, position=position_dodge(.9), size=0.7, alpha=0.6)+
  #stars
  annotate("text", x=as.numeric(stars_signif_IL2[,1]), y=as.numeric(stars_signif_IL2[,2]), label=stars_signif_IL2[,4], color="black", size=6)+
  ggsave("Spectratype_IL2_Dep_vs_indep.pdf", width=10, height=10, units = "cm")


  
ggplot(Data_IL15_graph, aes(x = Len, y = mean, fill = strain))+ 
  geom_bar(width = 0.95, stat = "identity") +
  xlab("") +
  ylab("") +
  labs(fill="", title = "") +
  scale_y_continuous(expand = c(0, 0)) +
  custom_theme() +
  guides(fill = guide_legend(reverse = FALSE))+
  expand_limits(y=c(0, 50))+ 
  scale_fill_manual(values=c("#156cff", "#ff492e"))+
  scale_x_continuous(breaks=seq(from=1.5, to=22.5, by=3), labels = seq(from=10, to=17))+
  geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=0.3, position=position_dodge(.9), size=0.7, alpha=0.6)+
  #stars
  annotate("text", x=as.numeric(stars_signif_IL15[,1]), y=as.numeric(stars_signif_IL15[,2]), label=stars_signif_IL15[,4], color="black", size=6)+
  ggsave("Spectratype_IL15_Dep_vs_indep.pdf", width=10, height=10, units = "cm")




#################
#  Means simple #
#################

Data_IL2_means[,1] <- c(1.2,2.2)
Data_IL15_means[,1] <- c(1.2,2.2)

### Means IL2
ggplot(Data_IL2_means, aes(x = coordinate, y = mean, fill = condition))+ 
  geom_bar(width = 0.5, stat = "identity") +
  xlab("") +
  ggtitle("")+
  ylab("") +
  custom_theme_2() +
  expand_limits(y=c(0, 20), x=c(0.8,2.3))+ 
  scale_fill_manual(values=c("#156cff", "#ff492e"))+
  scale_y_continuous(breaks=seq(from=0, to=20, by=5), expand = c(0, 0))+
  scale_x_continuous(breaks=c(1.2,2.2), label = c("Dependent", "Independent"))+
  geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=0.3, position=position_dodge(1), size=0.8, alpha=0.4)+
  #annotate("segment", x = 1.2, xend = 2.2, y = 18, yend = 18, colour = "grey40")+
  #stars
  annotate(geom="text", x=(stars_9[1,1]+0.7) , y=19, label="*",color="black", size=8)+
  annotate(geom="text", x=(stars_10[1,1]+0.7), y=19, label="**",color="black", size=8)+
  annotate(geom="text", x=(stars_11[1,1]+0.7), y=19, label="***",color="black", size=8)+
  annotate(geom="text", x=(stars_12[1,1]+0.7), y=19, label="****",color="black", size=8)+
  ggsave("Spectratype_mean_IL2_simple_Dep_vs_indep.tiff", width=8, height=11, units = "cm")

### Means IL15
ggplot(Data_IL15_means, aes(x = coordinate, y = mean, fill = condition))+ 
  geom_bar(width = 0.5, stat = "identity") +
  xlab("") +
  ggtitle("")+
  ylab("") +
  custom_theme_2() +
  expand_limits(y=c(0, 20), x=c(0.8,2.3))+ 
  scale_fill_manual(values=c("#156cff", "#ff492e"))+
  scale_y_continuous(breaks=seq(from=0, to=20, by=5), expand = c(0, 0))+
  scale_x_continuous(breaks=c(1.2,2.2), label = c("Dependent", "Independent"))+
  geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=0.2, position=position_dodge(1), size=0.8, alpha=0.4)+
  #annotate("segment", x = 1.2, xend = 2.2, y = 18, yend = 18, colour = "grey40")+
  annotate(geom="text", x=(stars_13[1,1]+0.7) , y=19, label="*",color="black", size=8)+
  annotate(geom="text", x=(stars_14[1,1]+0.7), y=19, label="**",color="black", size=8)+
  annotate(geom="text", x=(stars_15[1,1]+0.7), y=19, label="***",color="black", size=8)+
  annotate(geom="text", x=(stars_16[1,1]+0.7), y=19, label="****",color="black", size=8)+
  ggsave("Spectratype_mean_IL15_simple_Dep_vs_indep.tiff", width=8, height=11, units = "cm")


