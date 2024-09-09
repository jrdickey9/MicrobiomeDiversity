#PHYCO STRESS AKA NIKKI'S THESIS DATA

#installing packages
#install.packages("microbiome")
#BiocManager::install("microbiome")


#install.packages("ggvenn") # install via CRAN
library(ggvenn)
ggvenn(asv_lists)

#load packages
library(microbiome)
library(vegetarian)
library(phyloseq); packageVersion("phyloseq") #1.26.1
library(ggplot2); packageVersion("ggplot2") #3.4.1
library(gdata)
library(ecodist)
library(vegan)
library("car")
library(dplyr)
library(biomformat)

set.seed(34)

setwd("/Volumes/Dickey_External/PhycoStress/Dickey_files/qiime2_files/dadatable")
ASV_reads <- read_biom("Dada2FeatTab.biom")
ASV_table <- as.data.frame(as.matrix(biom_data(ASV_reads)))
#write.csv(ASV_table,"ASVtable.csv")
#ASV_table<-read.csv("ASVtable.csv")
otu_tab<-t(ASV_table) #36 291

#setwd
setwd("path name")

#reading in biom table that was converted to txt and then to csv.
#otu_tab<-read.csv("otu_tab_wtaxonomy.csv") #1379 x 38
#accession.numbs<-as.vector(otu_tab[,1])
#taxonomy<-as.vector(otu_tab[,38])
#rownames(otu_tab)<-accession.numbs
colnames(otu_tab)
otu_tab[1:36,1:3]
Phyco_OTU<-otu_tab #36 columns aka samps, rownames are accession IDs, colnames are samples. 
phyOTU<-Phyco_OTU
colnames(Phyco_OTU)
dim(phyOTU) #36 x 291

df.phyOTU<-as.data.frame(phyOTU)
df.phyOTU[1:10,1:10]
rownames(df.phyOTU)
colnames(df.phyOTU)

#reading in metadata file
setwd("path name")
metadata<-read.csv("metadata2.csv",header=TRUE) #28 x 5
dim(metadata) #36x5
colnames(metadata) 

#taxonomy
setwd("path name")
taxa3<-read.csv("taxa4.csv")
ph.headers<-c("ID", "Kingdom","Phylum","Class","Order","Family","Genus","species")
colnames(taxa3)<-ph.headers
rownames(taxa3) #just numbers
rownames(taxa3)<-taxa3[,1]
taxa3<-taxa3[,2:8]
colnames(taxa3)

#reading in phylogenetic tree
library("ape")
library(phytools)
library(castor)
library(doParallel)
setwd("path name")
dated.16Stree<-read.tree(file="tree.nwk") #reading in tree
is.rooted(dated.16Stree)

sample_names(dated.16Stree) #NULL
dated.16Stree$tip.label #for example "bf08d62b32cced86e829cba893bdf318" 

#creating phyloseq object
str(taxa3) #data.frame 291 x 7
taxa3.m<-as.matrix(taxa3) #VERY NECESSARY TO DO, DON'T SKIP. 
str(taxa3.m)
colnames(taxa3.m)
str(df.phyOTU) #data.frame 36 x 291
str(metadata) #data.frame 36 x 5

rownames(df.phyOTU)<-as.character(rownames(df.phyOTU))
colnames(df.phyOTU) #accession numbers
rownames(metadata)<-as.character(metadata[,1])
rownames(taxa3.m)<-as.character(rownames(taxa3.m)) #these the accession numbers
samp.names<-as.character(metadata[,1])

#to set up sample names to match (originally marked NULL)
sample_names(df.phyOTU)<-samp.names
sample_names(metadata)<-samp.names
sample_names(taxa3.m)<-samp.names
sample_names(dated.16Stree)<-samp.names

#to set up taxa names to match (originally marked NULL)
taxa_names(df.phyOTU)<-colnames(df.phyOTU)
taxa_names(dated.16Stree)<-colnames(df.phyOTU)
taxa_names(metadata)<-colnames(df.phyOTU)
taxa_names(taxa3.m)<-colnames(df.phyOTU)

#Here is the actual phyloseq object
phylo.phyco<-phyloseq(otu_table(df.phyOTU, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxa3.m), phy_tree(dated.16Stree))

#replace taxa names to short hand
taxa_names(phylo.phyco) <- paste0("ASV", seq(ntaxa(phylo.phyco)))
phylo.phyco@otu_table[1:10,1:10]
dim(phylo.phyco@otu_table) #36 x 1379
rowSums(phylo.phyco@otu_table)
colSums(phylo.phyco@otu_table)

#Making sure that there are no OTUs assigned to Eukaryota, Chloroplasts, or Mitochondria prior to rarefaction and standardization. Usually found in Kingdom, Class, Order, Family. Checkin' some seqs are sneaky.
table(tax_table(phylo.phyco)[, "Kingdom"], exclude = NULL) #6 unassigned, examine reads
table(tax_table(phylo.phyco)[, "Phylum"], exclude = NULL) 
table(tax_table(phylo.phyco)[, "Class"], exclude = NULL) 
table(tax_table(phylo.phyco)[, "Order"], exclude = NULL) #15 chloroplast 
table(tax_table(phylo.phyco)[, "Family"], exclude = NULL) #6 mitochondria
table(tax_table(phylo.phyco)[, "Genus"], exclude = NULL) 
table(tax_table(phylo.phyco)[, "species"], exclude = NULL) #2 unidentified

#######examining unassigned reads at the Kingdom level########
which(tax_table(phylo.phyco)[, "Kingdom"]=="Unassigned") #94, 189, 192, 256, 284, 286
phylo.phyco@tax_table[94,] 
rownames(taxa3.m)[94] #6283478ee4038101f96fba3ec2fe640f Mus musculus domesticus mitochondria
phylo.phyco@tax_table[189,] 
rownames(taxa3.m)[189] #e01ac2816dfa32ef6e1aba4a0f4e5bd0 #uncultured bacterium possibly chloroplast
phylo.phyco@tax_table[192,] 
rownames(taxa3.m)[192] #66463713cb5a51a7c6fcda6081371c02 #unculture bacterium
phylo.phyco@tax_table[256,] 
rownames(taxa3.m)[256] #e907974f73ffa6e286dcb4c5a022fb81 #uncultured bacterium
phylo.phyco@tax_table[284,] 
rownames(taxa3.m)[284] #18091fd090ed0c59e01959a5176e19ca #uncultured bacterium
phylo.phyco@tax_table[286,] 
rownames(taxa3.m)[286] #93532bfd4c9f2406ad783b009fc58735 #uncultured bacterium

p1<-subset_taxa(phylo.phyco,  !Kingdom %in% "Unassigned") #removing unknowns that asv seqs were blasted against NCBI
p2<-subset_taxa(p1,  !Order %in% "Chloroplast") #removing chloros
p3<-subset_taxa(p2,  !Family %in% "Mitochondria") #removing mitos

#p3@sam_data$sample_name
#removing axenic filtrate, pond 1, and pond 2 filtrates. 
PycoStressSamps<-subset_samples(p3, sample_name != "Pond1Filtrate" & sample_name != "Pond2Filtrate" & sample_name != "AxenicFiltrate" & sample_name != "PCRBlank1" & sample_name != "PCRBlank2" & sample_name != "PCRBlank3" & sample_name != "PCRBlank4" & sample_name != "ZymoPCR.01" & sample_name != "ZymoPCR.1" & sample_name != "ZymoPCR1" & sample_name != "ZymoPCR10")

dim(PycoStressSamps@sam_data) #25x5
#Beginning analyses here
phyco.tab<-PycoStressSamps@otu_table #to look at with controls use p3 object @ otu_table
dim(phyco.tab) #25 x 264
phyco.tab[1:10,1:10]
rowSums(phyco.tab)
mean(rowSums(phyco.tab))
median(rowSums(phyco.tab))
min(rowSums(phyco.tab))
max(rowSums(phyco.tab))
colSums(phyco.tab)
phyco.tab1<-phyco.tab[,-(which(colSums(phyco.tab)==0))]
dim(phyco.tab1) #25x151
colSums(phyco.tab1)
subset.obj<-colnames(phyco.tab1)
min(rowSums(phyco.tab1)) #39
#str(phyco.tab.m)
phyco.tab.df<-as.data.frame(phyco.tab1)
min(rowSums(phyco.tab.df)) #39
str(phyco.tab.df)
rdat<-rrarefy(phyco.tab.df,39) 
dim(rdat)

#standardize abundances into proportions. 
rowSums(rdat) #rarefaction worked
std.phyco.tab<-decostand(rdat,"total")
std.phyco.tab[1:10,1:10]
#write.csv(std.phyco.tab,"MicrobiomeDiv_Rarefied_PropAbund_ASVtable.csv")
dim(std.phyco.tab) #25 x 151 !!!!

#create quantitative Jaccard distance - takes into account relative abundance of OTUs vs presence/absence
drdat<-vegdist(std.phyco.tab,"jaccard")

#modeling
met.dat<-PycoStressSamps@sam_data #to look at with controls use p3 object @ sam_data
met.dat$Temperature_Treatment<-as.factor(met.dat$Temperature_Treatment)
met.dat$Nutrient_Treatment<-as.factor(met.dat$Nutrient_Treatment)
met.dat$Microbiome_Treatment<-as.factor(met.dat$Microbiome_Treatment)

mod<-dbrda(drdat~met.dat$Temperature_Treatment+met.dat$Nutrient_Treatment+met.dat$Microbiome_Treatment)

h<-how(nperm=10000)

anova(mod,permutations = h, by="margin") 

summary(mod)

#Weighted UniFrac distance
registerDoParallel(cores=4)

PycoStressSamps1<-prune_taxa(taxa=subset.obj, x=PycoStressSamps@phy_tree)
str(PycoStressSamps1)

Phylo.weight<-phyloseq(otu_table(rdat, taxa_are_rows=FALSE),sample_data(metadata), phy_tree(PycoStressSamps1))
dim(Phylo.weight@otu_table)
length(Phylo.weight@phy_tree$tip.label)

wdistUni<-UniFrac(physeq=Phylo.weight, weighted=TRUE, parallel=TRUE, fast=TRUE) #weighted #warning about randomly selected root for tree resolved. 

mod.Uni<-dbrda(wdistUni~met.dat$Temperature_Treatment+met.dat$Nutrient_Treatment+met.dat$Microbiome_Treatment)
anova(mod.Uni, permutations =h, by="margin")
summary(mod.Uni)

#venn diagram of shared ASVs
#subset samples by microbiome diversity, use non rarefied data
Venn_MD<-phyloseq(otu_table(phyco.tab1, taxa_are_rows=FALSE),sample_data(metadata), phy_tree(PycoStressSamps1))
dim(Venn_MD@otu_table) #25x151
dim(Venn_MD@sam_data) #25x5
colSums(Venn_MD@otu_table) #there are no ASVs that are singletons. There 

low.div<-subset_samples(Venn_MD, Microbiome_Treatment=="Low")
med.div<-subset_samples(Venn_MD, Microbiome_Treatment=="Medium")
high.div<-subset_samples(Venn_MD, Microbiome_Treatment=="High")

low.div.table<-low.div@otu_table
med.div.table<-med.div@otu_table
high.div.table<-high.div@otu_table

dim(Venn_MD@otu_table) #25x151
dim(low.div.table) #8x151
dim(med.div.table) #8x151
dim(high.div.table) #9x151

uniqueASVsLOW<-low.div.table[,-(which(colSums(low.div.table)==0))] #shouldnt change anything because this was removed earlier before creating quantitative jaccard distance. There are no zeros nor singletons. 
dim(uniqueASVsLOW) #8x55
LowASV<-colnames(uniqueASVsLOW)

uniqueASVsMED<-med.div.table[,-(which(colSums(med.div.table)==0))] #shouldnt change anything because this was removed earlier before creating quantitative jaccard distance. There are no zeros nor singletons. 
dim(uniqueASVsMED) #8x76
MedASV<-colnames(uniqueASVsMED)

uniqueASVsHIGH<-high.div.table[,-(which(colSums(high.div.table)==0))] #shouldnt change anything because this was removed earlier before creating quantitative jaccard distance. There are no zeros nor singletons. 
dim(uniqueASVsHIGH) #9x89
HighASV<-colnames(uniqueASVsHIGH)

##############using the objects above I want to do the following. Create a venn diagram that looks at the number of unique and shared taxa among microbiome diversity treatments for bacterial ASVs that are >1% of the total microbiome within each treatment group. I would need to find which are less than >1% and retain them if and only if they are >1% in another treatment group. 

#uniqueASVsLOW<-low.div.table[,-(which(colSums(low.div.table)==0))] #Line ran above.
dim(uniqueASVsLOW) #8 samples x 55 ASVs

#uniqueASVsMED<-med.div.table[,-(which(colSums(med.div.table)==0))] #Line ran above 
dim(uniqueASVsMED) #8 samples x 76 ASVs

#uniqueASVsHIGH<-high.div.table[,-(which(colSums(high.div.table)==0))] #Line ran above. 
dim(uniqueASVsHIGH) #9 samples x 89 ASVs

uniqueASVsLOW[1:8,1:10] #raw read numbers
Lowcols<-colnames(uniqueASVsLOW)
LowDivSumsL<-as.numeric(colSums(uniqueASVsLOW))
LowDiv.df<-as.data.frame(LowDivSumsL)
rownames(LowDiv.df)<-Lowcols
total_read_depthL<-sum(LowDiv.df[,1])
LowASV_PropAbund<-((LowDiv.df[,1])/total_read_depthL)*100
LowDiv.df<-cbind(LowDiv.df,LowASV_PropAbund)
#write.csv(LowDiv.df,file="LowDivUniqueASVs.csv")

uniqueASVsMED[1:8,1:10] #raw read numbers
MEDcols<-colnames(uniqueASVsMED)
MEDDivSumsM<-as.numeric(colSums(uniqueASVsMED))
MEDDiv.df<-as.data.frame(MEDDivSumsM)
rownames(MEDDiv.df)<-MEDcols
total_read_depthM<-sum(MEDDiv.df[,1])
MEDASV_PropAbund<-((MEDDiv.df[,1])/total_read_depthM)*100
MEDDiv.df<-cbind(MEDDiv.df,MEDASV_PropAbund)
#write.csv(MEDDiv.df,file="MEDDivUniqueASVs.csv")

uniqueASVsHIGH[1:9,1:10] #raw read numbers
HIGHcols<-colnames(uniqueASVsHIGH)
HIGHDivSumsH<-as.numeric(colSums(uniqueASVsHIGH))
HIGHDiv.df<-as.data.frame(HIGHDivSumsH)
rownames(HIGHDiv.df)<-HIGHcols
total_read_depthH<-sum(HIGHDiv.df[,1])
HIGHASV_PropAbund<-((HIGHDiv.df[,1])/total_read_depthH)*100
HIGHDiv.df<-cbind(HIGHDiv.df,HIGHASV_PropAbund)
#write.csv(HIGHDiv.df,file="HIGHDivUniqueASVs.csv")

#remove less than 1% and re-import as single csv table
common.asvs.microbiome.div<-read.csv(file="MicroDiv_CommonASVs.csv")
unique.low<-common.asvs.microbiome.div[which(common.asvs.microbiome.div$Treatment=="Low"),]
unique.med<-common.asvs.microbiome.div[which(common.asvs.microbiome.div$Treatment=="Medium"),]
unique.high<-common.asvs.microbiome.div[which(common.asvs.microbiome.div$Treatment=="High"),]

lowasv.common<-unique.low$ASV_Numb
medasv.common<-unique.med$ASV_Numb
highasv.common<-unique.high$ASV_Numb

asv_lists_common <- list("Low" = lowasv.common, "Medium" = medasv.common, "High" = highasv.common)

ggvenn(asv_lists_common, fill_color = c("#03045e","#ffb703", "#219ebc"), fill_alpha=.8, stroke_size = 0.8, set_name_size = 6, text_size = 6, show_percentage = FALSE)

#removing low microbiome samples. #need to use the 25x151 nonrarefied data object
PycoStressSamps_med_high<-subset_samples(Venn_MD, sample_name != "NMM6" & sample_name != "NMM36" & sample_name != "NMM31" & sample_name != "NMM26" & sample_name != "NMM21" & sample_name != "NMM16" & sample_name != "NMM11" & sample_name != "NMM1")

#Beginning analyses here
phyco.tab<-PycoStressSamps_med_high@otu_table
dim(phyco.tab) #17 x 151
phyco.tab[1:10,1:10]
rowSums(phyco.tab)
colSums(phyco.tab)
#remove singleton ASVs and ones 
phyco.tab1<-phyco.tab[,-(which(colSums(phyco.tab)==0))]
dim(phyco.tab1) #25x125
colSums(phyco.tab1)
subset.obj<-colnames(phyco.tab1)

#Rarefaction - to the minimum - which also is a low microbiome treatment. So some (aka the high microbiome samps are being rarefied to the lowest of low ones). #find the minimum using rowsums 
min(rowSums(phyco.tab)) #6123
#str(phyco.tab.m)
phyco.tab.df<-as.data.frame(phyco.tab)
min(rowSums(phyco.tab.df)) #6123
str(phyco.tab.df)
rdat<-rrarefy(phyco.tab.df,6123) 

#standardize abundances into proportions. 
rowSums(rdat) #rarefaction worked
std.phyco.tab.mh<-decostand(rdat,"total")
std.phyco.tab.mh[1:10,1:10]
dim(std.phyco.tab.mh) #17x151
#write.csv(std.phyco.tab,"MicrobiomeDiv_Rarefied_PropAbund_ASVtable.csv")


#create quantitative Jaccard distance - takes into account relative abundance of OTUs vs presence/absence
drdat<-vegdist(std.phyco.tab.mh,"jaccard")

#modeling
met.dat.mh<-PycoStressSamps_med_high@sam_data #to look at with controls use p3 object @ sam_data
met.dat.mh$Temperature_Treatment<-as.factor(met.dat.mh$Temperature_Treatment)
met.dat.mh$Nutrient_Treatment<-as.factor(met.dat.mh$Nutrient_Treatment)
met.dat.mh$Microbiome_Treatment<-as.factor(met.dat.mh$Microbiome_Treatment)


mod<-dbrda(drdat~met.dat.mh$Temperature_Treatment+met.dat.mh$Nutrient_Treatment+met.dat.mh$Microbiome_Treatment)

h<-how(nperm=10000)

anova(mod,permutations = h, by="margin") 

summary(mod)

registerDoParallel(cores=4)
PycoStressSamps1<-prune_taxa(taxa=subset.obj, x=PycoStressSamps_med_high@phy_tree)
str(PycoStressSamps1)

Phylo.weight<-phyloseq(otu_table(rdat, taxa_are_rows=FALSE),sample_data(metadata), phy_tree(PycoStressSamps1))
dim(Phylo.weight@otu_table) #17x125
length(Phylo.weight@phy_tree$tip.label) #125

wdistUni<-UniFrac(physeq=Phylo.weight, weighted=TRUE, parallel=TRUE, fast=TRUE) #weighted #warning about randomly selected root for tree resolved. 

mod.Uni<-dbrda(wdistUni~met.dat.mh$Temperature_Treatment+met.dat.mh$Nutrient_Treatment+met.dat.mh$Microbiome_Treatment)
anova(mod.Uni, permutations =h, by="margin")
summary(mod.Uni)

#####Alpha diversity analyses
alpha0<-apply(std.phyco.tab,1,vegetarian::d,q=0) #rarefied down to 39 reads. #25x151
alpha1<-apply(std.phyco.tab,1,vegetarian::d,q=1)
alpha0.mh<-apply(std.phyco.tab.mh,1,vegetarian::d,q=0) #rarefied down to 6,123 reads. #17x151
alpha1.mh<-apply(std.phyco.tab.mh,1,vegetarian::d,q=1)
nonrare.alpha0<-apply(phyco.tab1,1,vegetarian::d,q=0) #length needs to be 25x151

dim(met.dat) #25x6
dim(met.dat.mh) #17x6 

alpha0.mod<-lm(alpha0~met.dat$Temperature_Treatment+met.dat$Nutrient_Treatment+met.dat$Microbiome_Treatment)
alpha0.mh.mod<-lm(alpha0.mh~met.dat.mh$Temperature_Treatment+met.dat.mh$Nutrient_Treatment+met.dat.mh$Microbiome_Treatment)

anova(alpha0.mod)
anova(alpha0.mh.mod)


