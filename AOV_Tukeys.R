#Analysis of variance and Tukey's HSD

#Loading packages
library(vegan)
library(vegetarian)
library("car")
library("multcompView")
#library(multcomp)
#library(purrr)

#Import Data
setwd("path name")
Cell_Counts<-read.csv("Cell_Counts.csv")
FP_data<-read.csv("FP_data.csv")
TNP_data<-read.csv("TNP_data.csv")
Isotope_xenic<-read.csv("Isotope_data.csv")
axenic_isotopes<-read.csv(file="Isotope_Axenic_Data.csv")

#Organizing Cell_counts dataframe
newdat<-Cell_Counts[,5:9]
Shannon.In<-diversity(newdat,"shannon")
log.cell<-log(Cell_Counts$Total_Density)
Cell_Counts<-cbind(Cell_Counts,Shannon.In,log.cell)
x <- c(1, 2, 3, 4) #define vector
names(x) <- c('Oligotrophic', 'Mesotrophic', 'Eutrophic', 'Hypereutrophic') #provide names
Cell_Counts$Nutrient_Treatment[Cell_Counts$Nutrient_Treatment == "OG"]<-"Oligotrophic" 
Cell_Counts$Nutrient_Treatment[Cell_Counts$Nutrient_Treatment == "MS"]<-"Mesotrophic" 
Cell_Counts$Nutrient_Treatment[Cell_Counts$Nutrient_Treatment == "EU"]<-"Eutrophic" 
Cell_Counts$Nutrient_Treatment[Cell_Counts$Nutrient_Treatment == "HE"]<-"Hypereutrophic" 
Cell_Counts$Nutrient_Treatment<- factor(Cell_Counts$Nutrient_Treatment,levels = c('Oligotrophic', 'Mesotrophic', 'Eutrophic', 'Hypereutrophic'))
Cell_Counts$Microbiome_Treatment <- factor(Cell_Counts$Microbiome_Treatment,levels=c("Low","Medium","High"))
Cell_Counts$Temperature_Treatment<-as.factor(Cell_Counts$Temperature_Treatment)
str(Cell_Counts)

#Organizing FP_data dataframe
FP_data$Nutrient_Treatment[FP_data$Nutrient_Treatment == "OG"]<-"Oligotrophic" 
FP_data$Nutrient_Treatment[FP_data$Nutrient_Treatment == "MS"]<-"Mesotrophic" 
FP_data$Nutrient_Treatment[FP_data$Nutrient_Treatment == "EU"]<-"Eutrophic" 
FP_data$Nutrient_Treatment[FP_data$Nutrient_Treatment == "HE"]<-"Hypereutrophic" 
FP_data$Nutrient_Treatment<- factor(FP_data$Nutrient_Treatment,levels = c('Oligotrophic', 'Mesotrophic', 'Eutrophic', 'Hypereutrophic'))
FP_data$Nutrient_Treatment<- factor(FP_data$Nutrient_Treatment,levels = c('Oligotrophic', 'Mesotrophic', 'Eutrophic', 'Hypereutrophic'))
FP_data$Microbiome_Treatment <- factor(FP_data$Microbiome_Treatment,levels=c("Low","Medium","High"))
FP_data$Temperature_Treatment<-as.factor(FP_data$Temperature_Treatment)
str(FP_data)

#Organizing TNP_data dataframe
names(x) <- c('Oligotrophic', 'Mesotrophic', 'Eutrophic', 'Hypereutrophic')
TNP_data$Nutrient_Treatment[TNP_data$Nutrient_Treatment == "OG"]<-"Oligotrophic" 
TNP_data$Nutrient_Treatment[TNP_data$Nutrient_Treatment == "MS"]<-"Mesotrophic" 
TNP_data$Nutrient_Treatment[TNP_data$Nutrient_Treatment == "EU"]<-"Eutrophic" 
TNP_data$Nutrient_Treatment[TNP_data$Nutrient_Treatment == "HE"]<-"Hypereutrophic" 
TNP_data$Nutrient_Treatment<- factor(TNP_data$Nutrient_Treatment,levels = c('Oligotrophic', 'Mesotrophic', 'Eutrophic', 'Hypereutrophic'))
TNP_data$Microbiome_Treatment <- factor(TNP_data$Microbiome_Treatment,levels=c("Low","Medium","High"))
TNP_data$Temperature_Treatment<-as.factor(TNP_data$Temperature_Treatment)
str(TNP_data)

#µg/L = (Molar concentration (mol/L) x Molecular weight (g/mol) x 1,000,000)
TP_M<-TNP_data$TP*(0.000001)
TP_gpermol<-TP_M*30.97
TP_ugml<-TP_gpermol*1000000

TN_M<-TNP_data$TN*(0.000001)
TN_gpermol<-TN_M*14.01 
TN_ugml<-TN_gpermol*1000000

#attach back to TNP
TNP_data1<-cbind(TNP_data,TP_ugml,TN_ugml)

#Organizing Isotope_xenic dataframe 
names(x) <- c('Oligotrophic', 'Mesotrophic', 'Eutrophic', 'Hypereutrophic')
Isotope_xenic$Nutrient_Treatment[Isotope_xenic$Nutrient_Treatment == "OG"]<-"Oligotrophic" 
Isotope_xenic$Nutrient_Treatment[Isotope_xenic$Nutrient_Treatment == "MS"]<-"Mesotrophic" 
Isotope_xenic$Nutrient_Treatment[Isotope_xenic$Nutrient_Treatment == "EU"]<-"Eutrophic" 
Isotope_xenic$Nutrient_Treatment[Isotope_xenic$Nutrient_Treatment == "HE"]<-"Hypereutrophic" 
Isotope_xenic$Nutrient_Treatment<- factor(Isotope_xenic$Nutrient_Treatment,levels = c('Oligotrophic', 'Mesotrophic', 'Eutrophic', 'Hypereutrophic'))
Isotope_xenic$Microbiome_Treatment <- factor(Isotope_xenic$Microbiome_Treatment,levels=c("Low","Medium","High"))
Isotope_xenic$Temperature_Treatment<-as.factor(Isotope_xenic$Temperature_Treatment)
str(Isotope_xenic)

#AOVs
CD_aov<-aov(log.cell~Nutrient_Treatment*Microbiome_Treatment*Temperature_Treatment, data = Cell_Counts)
summary(CD_aov) #matches supp. info. 

FP_aov<-aov(Sample_Weight~Nutrient_Treatment*Microbiome_Treatment*Temperature_Treatment, data = FP_data)
summary(FP_aov) #matches supp. info. 

TN_aov<-aov(TN_ugml~Nutrient_Treatment*Microbiome_Treatment*Temperature_Treatment, data = TNP_data1)
summary(TN_aov) #matches supp. info. 

TP_aov<-aov(log(TP_ugml)~Nutrient_Treatment*Microbiome_Treatment*Temperature_Treatment, data = TNP_data1)
summary(TP_aov)

CN_ratio_aov<-aov(CN_Ratio~Nutrient_Treatment*Microbiome_Treatment*Temperature_Treatment, data = Isotope_xenic)
summary(CN_ratio_aov)

c13_aov<-aov(Delta_C~Nutrient_Treatment*Microbiome_Treatment*Temperature_Treatment, data = Isotope_xenic)
summary(c13_aov) #matches supp. info. 

n15_aov<-aov(Delta_N~Nutrient_Treatment*Microbiome_Treatment*Temperature_Treatment, data = Isotope_xenic)
summary(n15_aov) #matches supp. info. 

SD_aov<-aov(Shannon.In~Nutrient_Treatment*Microbiome_Treatment*Temperature_Treatment, data = Cell_Counts)
summary(SD_aov) #matches supp. info.

#Tukey's HSD using TukeyHSD
CD_THSD<-TukeyHSD(CD_aov)
multcompLetters4(CD_aov, CD_THSD)

FP_THSD<-TukeyHSD(FP_aov)
multcompLetters4(FP_aov, FP_THSD)

TN_THSD<-TukeyHSD(TN_aov)
multcompLetters4(TN_aov, TN_THSD)

TP_THSD<-TukeyHSD(TP_aov)
multcompLetters4(TP_aov, TP_THSD)

c13_THSD<-TukeyHSD(c13_aov)
multcompLetters4(c13_aov, c13_THSD)

n15_THSD<-TukeyHSD(n15_aov)
multcompLetters4(n15_aov, c13_THSD)

SD_THSD<-TukeyHSD(SD_aov)
multcompLetters4(SD_aov, SD_THSD)

#########Tukey's HSD within nutrient treatment
#log cell counts
Oligo_Cell_Counts<- Cell_Counts[which(Cell_Counts$Nutrient_Treatment=="Oligotrophic"),]
Oligo_CD_aov<-aov(log.cell~Microbiome_Treatment*Temperature_Treatment, data = Oligo_Cell_Counts)
Oligo_CD_THSD<-TukeyHSD(Oligo_CD_aov)
multcompLetters4(Oligo_CD_aov, Oligo_CD_THSD)

Meso_Cell_Counts<- Cell_Counts[which(Cell_Counts$Nutrient_Treatment=="Mesotrophic"),]
Meso_CD_aov<-aov(log.cell~Microbiome_Treatment*Temperature_Treatment, data = Meso_Cell_Counts)
Meso_CD_THSD<-TukeyHSD(Meso_CD_aov)
multcompLetters4(Meso_CD_aov, Meso_CD_THSD)

Eut_Cell_Counts<- Cell_Counts[which(Cell_Counts$Nutrient_Treatment=="Eutrophic"),]
Eut_CD_aov<-aov(log.cell~Microbiome_Treatment*Temperature_Treatment, data = Eut_Cell_Counts)
Eut_CD_THSD<-TukeyHSD(Eut_CD_aov)
multcompLetters4(Eut_CD_aov, Eut_CD_THSD)

Hyper_Cell_Counts<- Cell_Counts[which(Cell_Counts$Nutrient_Treatment=="Hypereutrophic"),]
Hyper_CD_aov<-aov(log.cell~Microbiome_Treatment*Temperature_Treatment, data = Hyper_Cell_Counts)
Hyper_CD_THSD<-TukeyHSD(Hyper_CD_aov)
multcompLetters4(Hyper_CD_aov, Hyper_CD_THSD)

#biomass
Oligo_FP_data<- FP_data[which(FP_data$Nutrient_Treatment=="Oligotrophic"),]
Oligo_FP_aov<-aov(Sample_Weight~Microbiome_Treatment*Temperature_Treatment, data = Oligo_FP_data)
Oligo_FP_THSD<-TukeyHSD(Oligo_FP_aov)
multcompLetters4(Oligo_FP_aov, Oligo_FP_THSD)

Meso_FP_data<- FP_data[which(FP_data$Nutrient_Treatment=="Mesotrophic"),]
Meso_FP_aov<-aov(Sample_Weight~Microbiome_Treatment*Temperature_Treatment, data = Meso_FP_data)
Meso_FP_THSD<-TukeyHSD(Meso_FP_aov)
multcompLetters4(Meso_FP_aov, Meso_FP_THSD)

Eut_FP_data<- FP_data[which(FP_data$Nutrient_Treatment=="Eutrophic"),]
Eut_FP_aov<-aov(Sample_Weight~Microbiome_Treatment*Temperature_Treatment, data = Eut_FP_data)
Eut_FP_THSD<-TukeyHSD(Eut_FP_aov)
multcompLetters4(Eut_FP_aov, Eut_FP_THSD)

Hyper_FP_data<- FP_data[which(FP_data$Nutrient_Treatment=="Hypereutrophic"),]
Hyper_FP_aov<-aov(Sample_Weight~Microbiome_Treatment*Temperature_Treatment, data = Hyper_FP_data)
Hyper_FP_THSD<-TukeyHSD(Hyper_FP_aov)
multcompLetters4(Hyper_FP_aov, Hyper_FP_THSD)

#total N
Oligo_TN_data<- TNP_data1[which(TNP_data1$Nutrient_Treatment=="Oligotrophic"),]
Oligo_TN_aov<-aov(TN_ugml~Microbiome_Treatment*Temperature_Treatment, data = Oligo_TN_data)
Oligo_TN_THSD<-TukeyHSD(Oligo_TN_aov)
multcompLetters4(Oligo_TN_aov, Oligo_TN_THSD)

Meso_TN_data<- TNP_data1[which(TNP_data1$Nutrient_Treatment=="Mesotrophic"),]
Meso_TN_aov<-aov(TN_ugml~Microbiome_Treatment*Temperature_Treatment, data = Meso_TN_data)
Meso_TN_THSD<-TukeyHSD(Meso_TN_aov)
multcompLetters4(Meso_TN_aov, Meso_TN_THSD)

Eut_TN_data<- TNP_data1[which(TNP_data1$Nutrient_Treatment=="Eutrophic"),]
Eut_TN_aov<-aov(TN_ugml~Microbiome_Treatment*Temperature_Treatment, data = Eut_TN_data)
Eut_TN_THSD<-TukeyHSD(Eut_TN_aov)
multcompLetters4(Eut_TN_aov, Eut_TN_THSD)

Hyper_TN_data<- TNP_data1[which(TNP_data1$Nutrient_Treatment=="Hypereutrophic"),]
Hyper_TN_aov<-aov(TN_ugml~Microbiome_Treatment*Temperature_Treatment, data = Hyper_TN_data)
Hyper_TN_THSD<-TukeyHSD(Hyper_TN_aov)
multcompLetters4(Hyper_TN_aov, Hyper_TN_THSD)

#total p
Oligo_TP_data<- TNP_data1[which(TNP_data1$Nutrient_Treatment=="Oligotrophic"),]
Oligo_TP_aov<-aov(TP_ugml~Microbiome_Treatment*Temperature_Treatment, data = Oligo_TP_data)
Oligo_TP_THSD<-TukeyHSD(Oligo_TP_aov)
multcompLetters4(Oligo_TP_aov, Oligo_TP_THSD)

Meso_TP_data<- TNP_data1[which(TNP_data1$Nutrient_Treatment=="Mesotrophic"),]
Meso_TP_aov<-aov(TP_ugml~Microbiome_Treatment*Temperature_Treatment, data = Meso_TP_data)
Meso_TP_THSD<-TukeyHSD(Meso_TP_aov)
multcompLetters4(Meso_TP_aov, Meso_TP_THSD)

Eut_TP_data<- TNP_data1[which(TNP_data1$Nutrient_Treatment=="Eutrophic"),]
Eut_TP_aov<-aov(TP_ugml~Microbiome_Treatment*Temperature_Treatment, data = Eut_TP_data)
Eut_TP_THSD<-TukeyHSD(Eut_TP_aov)
multcompLetters4(Eut_TP_aov, Eut_TP_THSD)

Hyper_TP_data<- TNP_data1[which(TNP_data1$Nutrient_Treatment=="Hypereutrophic"),]
Hyper_TP_aov<-aov(TP_ugml~Microbiome_Treatment*Temperature_Treatment, data = Hyper_TP_data)
Hyper_TP_THSD<-TukeyHSD(Hyper_TP_aov)
multcompLetters4(Hyper_TP_aov, Hyper_TP_THSD)

#delta C
Oligo_Delta_C_data<- Isotope_xenic[which(Isotope_xenic$Nutrient_Treatment=="Oligotrophic"),]
Oligo_Delta_C_aov<-aov(Delta_C~Microbiome_Treatment*Temperature_Treatment, data = Oligo_Delta_C_data)
Oligo_Delta_C_THSD<-TukeyHSD(Oligo_Delta_C_aov)
multcompLetters4(Oligo_Delta_C_aov, Oligo_Delta_C_THSD)

Meso_Delta_C_data<- Isotope_xenic[which(Isotope_xenic$Nutrient_Treatment=="Mesotrophic"),]
Meso_Delta_C_aov<-aov(Delta_C~Microbiome_Treatment*Temperature_Treatment, data = Meso_Delta_C_data)
Meso_Delta_C_THSD<-TukeyHSD(Meso_Delta_C_aov)
multcompLetters4(Meso_Delta_C_aov, Meso_Delta_C_THSD)

Eut_Delta_C_data<- Isotope_xenic[which(Isotope_xenic$Nutrient_Treatment=="Eutrophic"),]
Eut_Delta_C_aov<-aov(Delta_C~Microbiome_Treatment*Temperature_Treatment, data = Eut_Delta_C_data)
Eut_Delta_C_THSD<-TukeyHSD(Eut_Delta_C_aov)
multcompLetters4(Eut_Delta_C_aov, Eut_Delta_C_THSD)

Hyper_Delta_C_data<- Isotope_xenic[which(Isotope_xenic$Nutrient_Treatment=="Hypereutrophic"),]
Hyper_Delta_C_aov<-aov(Delta_C~Microbiome_Treatment*Temperature_Treatment, data = Hyper_Delta_C_data)
Hyper_Delta_C_THSD<-TukeyHSD(Hyper_Delta_C_aov)
multcompLetters4(Hyper_Delta_C_aov, Hyper_Delta_C_THSD)

#delta N
Oligo_Delta_N_data<- Isotope_xenic[which(Isotope_xenic$Nutrient_Treatment=="Oligotrophic"),]
Oligo_Delta_N_aov<-aov(Delta_N~Microbiome_Treatment*Temperature_Treatment, data = Oligo_Delta_N_data)
Oligo_Delta_N_THSD<-TukeyHSD(Oligo_Delta_N_aov)
multcompLetters4(Oligo_Delta_N_aov, Oligo_Delta_N_THSD)

Meso_Delta_N_data<- Isotope_xenic[which(Isotope_xenic$Nutrient_Treatment=="Mesotrophic"),]
Meso_Delta_N_aov<-aov(Delta_N~Microbiome_Treatment*Temperature_Treatment, data = Meso_Delta_N_data)
Meso_Delta_N_THSD<-TukeyHSD(Meso_Delta_N_aov)
multcompLetters4(Meso_Delta_N_aov, Meso_Delta_N_THSD)

Eut_Delta_N_data<- Isotope_xenic[which(Isotope_xenic$Nutrient_Treatment=="Eutrophic"),]
Eut_Delta_N_aov<-aov(Delta_N~Microbiome_Treatment*Temperature_Treatment, data = Eut_Delta_N_data)
Eut_Delta_N_THSD<-TukeyHSD(Eut_Delta_N_aov)
multcompLetters4(Eut_Delta_N_aov, Eut_Delta_N_THSD)

Hyper_Delta_N_data<- Isotope_xenic[which(Isotope_xenic$Nutrient_Treatment=="Hypereutrophic"),]
Hyper_Delta_N_aov<-aov(Delta_N~Microbiome_Treatment*Temperature_Treatment, data = Hyper_Delta_N_data)
Hyper_Delta_N_THSD<-TukeyHSD(Hyper_Delta_N_aov)
multcompLetters4(Hyper_Delta_N_aov, Hyper_Delta_N_THSD)

#CN Ratio
Oligo_CN_Ratio_data<- Isotope_xenic[which(Isotope_xenic$Nutrient_Treatment=="Oligotrophic"),]
Oligo_CN_Ratio_aov<-aov(CN_Ratio~Microbiome_Treatment*Temperature_Treatment, data = Oligo_CN_Ratio_data)
Oligo_CN_Ratio_THSD<-TukeyHSD(Oligo_CN_Ratio_aov)
multcompLetters4(Oligo_CN_Ratio_aov, Oligo_CN_Ratio_THSD)

Meso_CN_Ratio_data<- Isotope_xenic[which(Isotope_xenic$Nutrient_Treatment=="Mesotrophic"),]
Meso_CN_Ratio_aov<-aov(CN_Ratio~Microbiome_Treatment*Temperature_Treatment, data = Meso_CN_Ratio_data)
Meso_CN_Ratio_THSD<-TukeyHSD(Meso_CN_Ratio_aov)
multcompLetters4(Meso_CN_Ratio_aov, Meso_CN_Ratio_THSD)

Eut_CN_Ratio_data<- Isotope_xenic[which(Isotope_xenic$Nutrient_Treatment=="Eutrophic"),]
Eut_CN_Ratio_aov<-aov(CN_Ratio~Microbiome_Treatment*Temperature_Treatment, data = Eut_CN_Ratio_data)
Eut_CN_Ratio_THSD<-TukeyHSD(Eut_CN_Ratio_aov)
multcompLetters4(Eut_CN_Ratio_aov, Eut_CN_Ratio_THSD)

Hyper_CN_Ratio_data<- Isotope_xenic[which(Isotope_xenic$Nutrient_Treatment=="Hypereutrophic"),]
Hyper_CN_Ratio_aov<-aov(CN_Ratio~Microbiome_Treatment*Temperature_Treatment, data = Hyper_CN_Ratio_data)
Hyper_CN_Ratio_THSD<-TukeyHSD(Hyper_CN_Ratio_aov)
multcompLetters4(Hyper_CN_Ratio_aov, Hyper_CN_Ratio_THSD)

#shanonnons 
Oligo_SD_aov<-aov(Shannon.In~Microbiome_Treatment*Temperature_Treatment, data = Oligo_Cell_Counts)
Oligo_SD_THSD<-TukeyHSD(Oligo_SD_aov)
multcompLetters4(Oligo_SD_aov, Oligo_SD_THSD)


Meso_SD_aov<-aov(Shannon.In~Microbiome_Treatment*Temperature_Treatment, data = Meso_Cell_Counts)
Meso_SD_THSD<-TukeyHSD(Meso_SD_aov)
multcompLetters4(Meso_SD_aov, Meso_SD_THSD)


Eut_SD_aov<-aov(Shannon.In~Microbiome_Treatment*Temperature_Treatment, data = Eut_Cell_Counts)
Eut_SD_THSD<-TukeyHSD(Eut_SD_aov)
multcompLetters4(Eut_SD_aov, Eut_SD_THSD)


Hyper_SD_aov<-aov(Shannon.In~Microbiome_Treatment*Temperature_Treatment, data = Hyper_Cell_Counts)
Hyper_SD_THSD<-TukeyHSD(Hyper_SD_aov)
multcompLetters4(Hyper_SD_aov, Hyper_SD_THSD)
