#Code for calculate alpha and beta diversity of bacterial community associatead wit parrotfish feces 

library(phyloseq)
library(devtools)
library(ggplot2)
library(readxl)
library(dplyr)
library(ape)
library(gridExtra)
library(knitr)
library(vegan)
library(tidyverse)
library(scales)
library(grid)
library(reshape2)

#Files from qiime2 for create a phyloseq object
asvpf<- read_excel("D:/Users/Usuario/Documents/PARROT FISH PAPER/IMPORTANT FILES/asvloro.xlsx")
taxpf<- read_excel("D:/Users/Usuario/Documents/PARROT FISH PAPER/IMPORTANT FILES/taxoloro.xlsx")
samplespf <- read_excel("D:/Users/Usuario/Documents/PARROT FISH PAPER/IMPORTANT FILES/METADATAloro.xlsx")
treepf<-read.tree("D:/Users/Usuario/Documents/PARROT FISH PAPER/loro/exportedtreemc17/treemc17.nwk")

asvpf <- data.frame(asvpf, row.names = 1)
taxpf <- data.frame(taxpf, row.names = 1)
samplespf <- data.frame(samplespf, row.names = 1) 
asv_pf <- as.matrix(asvpf)
tax_pf <- as.matrix(taxpf)

#Create a phyloseq object
ASVPF = otu_table(asv_pf, taxa_are_rows = TRUE)
TAXPF = tax_table(tax_pf)
samplespf = sample_data(samplespf)
TREPF = phy_tree(treepf)
phpf <- phyloseq(ASVPF, TAXPF, samplespf, TREPF)
phpf

#Modified Protocol from Callahan et al. 2016 (Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses)

# Create table, number of features for each phyla
table(tax_table(phpf)[, "Phylum"], exclude = NULL)
pf1 <- subset_taxa(phpf, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
pf1
# Compute prevalence of each feature, store as data.frame
prevpf = apply(X = otu_table(pf1),
               MARGIN = ifelse(taxa_are_rows(pf1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevpf = data.frame(Prevalence = prevpf,
                    TotalAbundance = taxa_sums(pf1),
                    tax_table(pf1))
plyr::ddply(prevpf, "Phylum", function(pf1){cbind(mean(pf1$Prevalence),sum(pf1$Prevalence))})
# Define phyla to filter
filterPhylapp = c("D_1__Zixibacteria", "D_1__WS2","D_1__Thaumarchaeota","D_1__Schekmanbacteria","D_1__Omnitrophicaeota","D_1__Nitrospinae","D_1__Nanoarchaeaeota","D_1__Modulibacteria","D_1__Fibrobacteres","D_1__Euryarchaeota","D_1__Elusimicrobia","D_1__Deinococcus-Thermus","D_1__Dadabacteria","D_1__Crenarchaeota","D_1__Altiarchaeota","D_1__Calditrichaeota","D_1__Deferribacteres","D_1__WPS-2")
# Filter entries with unidentified Phylum.
pf2 = subset_taxa(pf1, !Phylum %in% filterPhylapp)
pf2
# Subset to the remaining phyla
prevdf1 = subset(prevpf, Phylum %in% get_taxa_unique(pf2, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(phpf),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.1 * nsamples(phpf)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
pf3 = prune_taxa(keepTaxa, phpf)
pf3

#For create a Rarefaction curve
rarecurve(t(otu_table(pf3)), step=50, cex=0.5)
ps.rarefied = rarefy_even_depth(pf3, rngseed=1, sample.size=0.9*min(sample_sums(pf3)), replace=F)
ps.rarefied

#Alpha diversity and plots

alpha.diversity <- estimate_richness(pf3, measures = c("Observed","Chao1", "Shannon"))
alpha.diversity#table with values

data <- cbind(sample_data(pf3), alpha.diversity)
phy.anova <- aov(Chao1 ~ Type, data)#doing with chao1 and observed
summary(phy.anova)#significant effect of host specie in richness

#Pairwise comparisons between groups and correct for multiple comparisons, we run Tukey's honest significance test of our ANOVA.

TukeyHSD(phy.anova)#For 

##Beta diversity and NMDS ordination 

ps_bray <- phyloseq::distance(pf3, method = "bray")

ps_nmds <- ordinate(
  physeq = pf3, 
  method = "NMDS", 
  distance = "bray"
)
plot_ordination(
  physeq = pf3,
  ordination = ps_nmds,
  color = "Type",
  shape = "Type",
  title = "NMDS bacterial Communities"
) + 
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4","#C77CFF")
  ) +
  geom_point(aes(color = Type), alpha = 1, size = 3, shape = 16) 
 
##PERMANOVA for evaluate differences between the four types of samples using Bray Curtis dissimilarities

metadata <- as(sample_data(pf3), "data.frame")

adonis(distance(pf3, method="bray") ~ Type,
       data = metadata)
