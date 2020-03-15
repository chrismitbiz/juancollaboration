##############R ANALYSIS PREPARATION##############
#Set working directory
setwd("~/Documents/Study/LaTrobe/Research/phD/Juan_project") #set your own directory, store all files there#Load packages
library(vegan)          #Community Ecology Package: for metamds, rda (pca), cca (ca), capscale, anosim, adonis, vectors, cluster
library(labdsv)         #Ordination and Multivariate Analysis for Ecology: for pco, pca, ca
library(dplyr)
library(phyloseq)  
#packageVersion('phyloseq')
library(qiime2R)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggfortify)
library(ggpubr)
library(ggrepel)
library(factoextra)
library(Rmisc)
library(reshape2)


theme_set(theme_bw())  



treatments<-read.csv("treatments.csv",header = TRUE)
OTUtable <- read.csv("OTUdata.csv",header = TRUE, row.names = 1)
OTUtable <- t(OTUtable)
OTUtable_abu <- OTUtable / rowSums(OTUtable) #relative abundances
OTUtable_abu_log <- log(OTUtable_abu + 1) 

rowSums(OTUtable_abu)

#NMDS
#Data_abu <- OTUtable_abu[match(treatments$Sample,row.names(OTUtable_abu)),] 

## NMDS using Bray-Curtis dissimilarities
nmds_data <- metaMDS(comm = OTUtable_abu, distance = "bray")
## get x and y coordinates
nmds_data_df <- as.data.frame(nmds_data[["points"]])  

#head(area_nmds_df)
nmds_data_df <- cbind(nmds_data_df,treatments = treatments[["Treatments"]],   #make sure the treatments are added to df
                       Soil = treatments[["Soil.type"]])
#head(nmds_data_df)
#the same for the species scores
#species_nmds_df <- as.data.frame(nmds_data[["species"]])  
#species_nmds_df$substance <- rownames(species_nmds_df)  # create a column of species, from the rownames of species.scores
#head(species_nmds_df)

#ratio - looking at ordination - site scores only
ggplot(data = nmds_data_df,aes(MDS1,MDS2, color = Soil, #label = treatments,
                          size = 3, shape = treatments, palette = "simpsons")) +
  geom_point() + 
  theme_void() +  theme_bw()   
  #geom_text() 


##permanova
disp <- betadisper(vegdist(log1p(ARISAdataFUN)), envdata$Property) #across all samples
disp <- betadisper(vegdist(log1p(ARISAdataBAC_sel)), env_select$Property) #on selected data
plot(disp)  
boxplot(disp)   
anova(disp)   

permutest(disp) #to check if group median distances differ between groups
#if there are differences in multivariate dispersion between the groups i.e. p values are low and visually they separate well
#then PERMANOVA is not reliable (Warton, Wright & Wang 2011)
TukeyHSD(disp)


mod1 <- adonis(log1p(OTUtable_abu) ~ Soil.type*Treatments, data = treatments, method = "bray") 
mod1









