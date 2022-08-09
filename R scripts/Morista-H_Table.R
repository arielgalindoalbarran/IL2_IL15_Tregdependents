#!/usr/bin/env Rscript
##############################################################################
##                           SetUP                                          ##
##############################################################################
#Title: IL2-IL15 Treg dependent:Morista-H table
#Description: 
#Installing and loading required packages

if(!require(divo))
{install.packages("divo")
}
if(!require(tidyr))
{install.packages("tidyr")
}

if(!require(abdiv))
{install.packages("abdiv")
}
library(divo)
library(tidyr)
library(abdiv)
#################################
#Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################################################################
##                           Initial Setup                                 ##
##############################################################################

#READ files in the directory:
Files_data <- as.data.frame(list.files(getwd(), pattern = ".txt"))
colnames(Files_data) <- "files"
Files_data$files <- as.character(Files_data$files)

#Do the short name and sort
tmp_1 <- as.data.frame(gsub(".txt", "", Files_data[,1]))
colnames(tmp_1) <- c("cond")
tmp_1 <- tmp_1 %>% separate(cond, c("A", "B"), sep = "-")

Files_data$number <- tmp_1$A
Files_data$strain <- tmp_1$B
Files_data$number <- as.numeric(Files_data$number)
Files_data <- Files_data[order(Files_data$strain, Files_data$number, decreasing = FALSE),]
Files_data$shortname <- paste(Files_data$strain, c(1:4, 1:4, 1:4), sep = "_")
rownames(Files_data) <- NULL
###############################################################################
##                          Process Data ALL SAMPLES                         ##
###############################################################################

#Do the table
for (i in 1:nrow(Files_data)) {
  if (i != 1){ 
  tmp_1 <- read.table(Files_data[i,1], header = TRUE)
  clonotypes <-  c(clonotypes,paste(tmp_1 [,5], tmp_1 [,4], tmp_1 [,7], sep = "_")  ) 
  } else {
    tmp_1 <- read.table(Files_data[i,1], header = TRUE)
    clonotypes <-  paste(tmp_1 [,5], tmp_1 [,4], tmp_1 [,7], sep = "_")}
}
clonotypes <- unique(clonotypes)
table_comp <- matrix(data = 0, nrow = length(clonotypes), ncol=nrow(Files_data)+1, )
colnames(table_comp) <- c("clonotype_VAAJ", Files_data$shortname)
table_comp[,1] <- clonotypes

#Fill the table
for (i in 1:nrow(Files_data)) {
  tmp_1 <- read.table(Files_data[i,1], header = TRUE)
  tmp_1 <- as.data.frame(tmp_1 [,1], paste(tmp_1 [,5], tmp_1 [,4], tmp_1 [,7], sep = "_"))
  tmp_1$VAAJ <- rownames(tmp_1)
  rownames(tmp_1) <- NULL
  colnames(tmp_1) <- c("Counts", "VAAJ")
  tmp_1 <- as.data.frame(aggregate(Counts ~ VAAJ , data = tmp_1, sum))
  indx <- match(tmp_1$VAAJ, table_comp[,1], nomatch = 0) 
  table_comp[indx,i+1] <- tmp_1$Counts
   }

table_comp_Matrix <- as.matrix(table_comp)
table_comp_Matrix <- table_comp_Matrix[,-1]
for (i in 1:ncol(table_comp_Matrix)) {
  table_comp_Matrix[,i] <- as.numeric(table_comp_Matrix[,i])   
}
class(table_comp_Matrix) <- "numeric"

result <- mh(table_comp_Matrix, csv_output = "Morista-H_Table", PlugIn = TRUE)

