#atttempting to work with sequencing data in order to do analysis for UROP poster
#reading off of Michelle Berry's Phyloseq Tutorial
#install necessary packages
install.packages("vegan")
source("http://bioconductor.org/biocLite.R")
biocLite("phyloseq")
##########################################################################################
#load necessary packages 
library(ggplot2)
library(magrittr)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(vegan)
##########################################################################################
##saved Final for later use
save(Final, file = "Final.RData")

#set-up working directory
setwd("C:/Users/morga/Documents/UROP/Denef/Morgan/")

#import data
sharedfile = "Sed_QM.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared"
taxfile = "Sed_QM.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy"
mapfile = "Metadata_Table.csv"

#import mothur data and metadata
mothur_data <- import_mothur(mothur_shared_file = sharedfile, mothur_constaxonomy_file = taxfile)
metadata <- read.csv(mapfile)
#convert metadata to phyloseq format
metadata <- sample_data(metadata)
rownames(metadata) <- metadata$New_Name

##########################################################################################
#alterations to metadata_table
metadata$Month_Collected <- factor(metadata$Month_Collected, levels = c("April", "May", "June", 
                                                                        "July", "August", "September", 
                                                                        "October"))
############################################################################################

#merge metadata into mothur file
moth.merge = merge_phyloseq(mothur_data, metadata)
#15 samples lost when metadata and mothur_data were merged together. 11 were blanks and controls. 4 were samples
#missing samples: FM45.Ts, FM45.ZM, MLB.B1.S515.2, MLB.B1.S915

#reformat taxonomy file column names and rename
colnames(tax_table(moth.merge))
colnames(tax_table(moth.merge)) <- c("Kingdom", "Phylum", "Class", "Order", 
                                     "Family", "Genus", "Rank7", "Rank8")

#filter out samples that I do not want to use, on the tutorial it says to use sample type
#but for my metadata table I made a column named "Use" that has the same affect
moth.sub <- subset_samples(moth.merge, Use == "Sample")
moth.sub <- prune_taxa(taxa_sums(moth.sub) > 0, moth.sub) #unsure of what this step does
########################################################################################
#make new data.frame by removing any eukaryotic cells from the (tax_table(moth.sub))
Final <-
  moth.sub %>%
  subset_taxa(Kingdom == "Bacteria" &
                Family != "mitochondria" &
                Class != "Chloroplast")

#Scale the reads for "Final" by using code originally from Michelle Berry 

#Begin graphing with set that has removed eukaryotic (should not be different)
theme_set(theme_bw())

ggplot(data.frame(sum = sample_sums(Final)), aes(sum)) + 
  geom_histogram(color = "black", fill = "indianred", stat = "bin", binwidth = 500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") + 
  ylab("")

#Graph with eukaryotic cells included (we don't expect this to be different)
ggplot(data.frame(sum = sample_sums(moth.sub)), aes(sum)) +
  geom_histogram(color = "black", fill = "blue") +
  ggtitle("Distribution of Sample Sequencing Depth with Chloroplasts") +
  ylab("") + xlab("Read Counts")
#graphs did differ, meaning that there were some cholorplast/eukaryotic reads involved
#########################################################################################################################

#########################################################################################################################
#starting into the graphing of a stacked bar plot  of phyla >5%
sed_phylum5 <- Final %>%
  subset_samples(Type == "Sediment" & Cruise == "Transect") %>%
  tax_glom(taxrank = "Phylum") %>%                  
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                        
  filter(Abundance > 0.05) %>%                        
  arrange(Phylum)

#set plylum colors
phylum.colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")

ggplot(sed_phylum5, aes(x = Month_Collected, y = Abundance, fill = Phylum)) + 
 # facet_grid(Station~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum.colors) +
  #scale_x_discrete(
    #breaks = c("7/8", "8/4", "9/2", "10/6"),
    #labels = c("Jul", "Aug", "Sep", "Oct"), 
    #drop = FALSE
#  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 5%) \n") +
  ggtitle("Phylum Composition of Lake Michigan Sediment \n Bacterial Communities by Month Collected") 
#######################################################################################
#stacked bar blot of sediment phyla >2%
sed_phylum2 <- Final %>%
  subset_samples(Type == "Sediment" & Cruise == "Transect") %>%
  tax_glom(taxrank = "Phylum") %>%                  
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                        
  filter(Abundance > 0.02) %>%                        
  arrange(Phylum)

ggplot(sed_phylum2, aes(x = Month_Collected, y = Abundance, fill = Phylum)) + 
  # facet_grid(Station~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum.colors) +
  #scale_x_discrete(
  #breaks = c("7/8", "8/4", "9/2", "10/6"),
  #labels = c("Jul", "Aug", "Sep", "Oct"), 
  #drop = FALSE
  #  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Lake Michigan Sediment \n Bacterial Communities by Month Collected")
######################################################################################
#stacked bar plot of mussel samples by 2% pyla
mus_phylum2 <- Final %>%
  subset_samples(Sed_or_Mus == "Mussel" & Cruise == "Transect") %>%
  tax_glom(taxrank = "Phylum") %>%                  
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                        
  filter(Abundance > 0.02) %>%                        
  arrange(Phylum)

ggplot(mus_phylum2, aes(x = Type, y = Abundance, fill = Phylum)) + 
  # facet_grid(Station~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum.colors) +
  #scale_x_discrete(
  #breaks = c("7/8", "8/4", "9/2", "10/6"),
  #labels = c("Jul", "Aug", "Sep", "Oct"), 
  #drop = FALSE
  #  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Lake Michigan Quagga Mussel \n Bacterial Communities by Mussel Part")
#########################################################################################################

#########################################################################################################  
#Making a plot based on increasing read number to determine minimum to help normalize sample reads ()
sums_FWDB <- data.frame(colSums(otu_table(Final)))
colnames(sums_FWDB) <- "Sample_TotalSeqs"
sums_FWDB$sample <- row.names(sums_FWDB)
sums_FWDB <- arrange(sums_FWDB, Sample_TotalSeqs)

ggplot(sums_FWDB, aes(x=reorder(sample, Sample_TotalSeqs), y = Sample_TotalSeqs)) + 
  ylab("Number of Sequences per Sample") +
  geom_bar(stat = "identity", colour="black",fill="cornflowerblue")  + xlab("Sample Name") + 
  ggtitle("Total Number of Sequences per Sample") + 
  theme(axis.text.x = element_text(colour = "black", size=6, angle=45, hjust = 1, vjust = 1))
########################################################################################################

#################################################################################################################
#Plotting Ordinations
#scale reads to even depth/// This hhas already been done in the Formula section
Full.scale.format <-
  scale_Final %>%
  subset_samples(Use == "Sample") %>%
  scale_reads(round = "round")

Full.pcoa <- ordinate(physeq = scale_Final, 
                      method = "PCoA", 
                      distance = "bray")

#making plot
plot_ordination(
  physeq = scale_Final,
  ordination = Full.pcoa,
  color = "Type",
  shape = "Cruise",
  title = "PCoA of Lake Michigan and Muskegon Lake Samples"
) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta")
  ) +
  geom_point(aes(color = Type), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

###################################################PCoA of Sediment only
SedimentPCoA <-
  scale_Final %>%
  subset_samples(Type == "Sediment") %>%
  scale_reads(round = "round")

Full_Sed_PCoA <- ordinate(physeq = SedimentPCoA, 
                      method = "PCoA", 
                      distance = "bray")

#making plot
plot_ordination(
  physeq = SedimentPCoA,
  ordination = Full_Sed_PCoA,
  color = "Lake",
  shape = "Cruise",
  title = "PCoA of Lake Michigan and Muskegon Lake Sediment Samples"
) + 
  scale_color_manual(values = c("darkorchid3", "red", "#ffae19",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta")
  ) +
  geom_point(aes(color = Lake), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5)
######################################################################
#Mussel PCoA
MusselPCoA <-
  scale_Final %>%
  subset_samples(Sed_or_Mus == "Mussel") %>%
  scale_reads(round = "round")

Full_Mus_PCoA <- ordinate(physeq = MusselPCoA, 
                          method = "PCoA", 
                          distance = "bray")

#making plot
plot_ordination(
  physeq = MusselPCoA,
  ordination = Full_Mus_PCoA,
  color = "Type",
 # shape = "Mussel.Classification",
  title = "PCoA of Lake Michigan Mussel Samples"
) + 
  scale_color_manual(values = c("darkorchid3", "red", "#ffae19",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta")
  ) +
  geom_point(aes(color = Type), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5)
###############################################################################################################

###############################################################################################################
#Permanova
set.seed(1)

# Calculate bray curtis distance matrix
final_bray <- phyloseq::distance(scale_Final, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(Final))

# Adonis test
adonis(final_bray ~ Type, data = sampledf)

beta <- betadisper(final_bray, sampledf$Type)
permutest(beta)
###############################################################################################################

########FUNCTIONS####################################################################################################
matround <- function(x){trunc(x+0.5)}
scale_Final <-scale_reads(physeq = Final, n = 5000, round = "matround")

  scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "floor") {
    
    # transform counts to n
    physeq.scale <- transform_sample_counts(physeq, 
                                            function(x) {(n * x/sum(x))}
    )
    
    # Pick the rounding functions
    if (round == "floor"){
      otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
    } else if (round == "round"){
      otu_table(physeq.scale) <- round(otu_table(physeq.scale))
    } else if (round == "matround"){
      otu_table(physeq.scale) <- matround(otu_table(physeq.scale))
    }
    
    # Prune taxa and return new phyloseq object
    physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
    return(physeq.scale)
  }  
