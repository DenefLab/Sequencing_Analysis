#atttempting to work with sequencing data in order to do analysis for UROP poster

#reading off of Michelle Berry's Phyloseq Tutorial
#install necessary packages
install.packages("vegan")
source("http://bioconductor.org/biocLite.R")
biocLite("phyloseq")
##########################################################################################
#load necessary packages 
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(vegan)
library(tidyverse)
##########################################################################################

#import data
sharedfile = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.shared"
taxfile = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.0.03.cons.taxonomy"
mapfile = "Metadata_Table.csv"

#import mothur data and metadata
mothur_data <- import_mothur(mothur_shared_file = sharedfile, mothur_constaxonomy_file = taxfile)
metadata <- read.csv(mapfile) %>%
  select(-Read_Number)
#Rename data so that mothur data matches metadata
metadata$New_Name <- str_replace_all(string = metadata$New_Name,
                                     pattern = "\\.", 
                                     replacement = "_")
metadata$Original_Name <- str_replace_all(string = metadata$Original_Name,
                                     pattern = "\\.", 
                                     replacement = "_")
#convert metadata to phyloseq format
metadata <- sample_data(metadata)
rownames(metadata) <- metadata$Original_Name

##########################################################################################
#alterations to metadata_table
metadata$Month_Collected <- factor(metadata$Month_Collected, levels = c("April", "May", "June", 
                                                                        "July", "August", "September", 
                                                                        "October"))
metadata$Mussel_Density <- factor(metadata$Mussel_Density, levels = c("Very Few", "Few", "Few-Moderate", 
"Moderate", "Moderate-high", "Many", "High"))
metadata$Mussel_Density <- as.factor(metadata$Mussel_Density)
metadata$Depth_in_Meters <- as.factor(metadata$Depth_in_Meters)
############################################################################################




#############################################################################################
################################ Merge sample data with taxonomy data ######################
#merge metadata into mothur file
moth.merge = merge_phyloseq(mothur_data, metadata)
#15 samples lost when metadata and mothur_data were merged together. 11 were blanks and controls. 4 were samples
#missing samples: NFW_Blank and MLBB1_S916, one of the samples says CONTAMINATED, so assuming it is this one.

#reformat taxonomy file column names and rename
colnames(tax_table(moth.merge))
colnames(tax_table(moth.merge)) <- c("Domain", "Phylum", "Class", "Order", 
                                     "Family", "Genus", "Species")

#filter out samples that I do not want to use, on the tutorial it says to use sample type
#but for my metadata table I made a column named "Use" that has the same affect
moth.sub <- subset_samples(moth.merge, Use == "Sample")
moth.sub <- prune_taxa(taxa_sums(moth.sub) > 0, moth.sub) #unsure of what this step does

#make new data.frame by removing any eukaryotic cells from the (tax_table(moth.sub))
Final <-
  moth.sub %>%
  subset_taxa(Domain == "Bacteria" &
                Family != "mitochondria" &
                Class != "Chloroplast")

#Scale the reads for "Final" by using code originally from Michelle Berry 
#functions for "matround" and "scale_reads" found at end of document
scale_Final <-scale_reads(physeq = Final,n=5000, round = "matround")

###########################################################################################################################



############################################################################################################################
########################### Graph Sample Reads #############################################################################
#Begin graphing with set that has removed eukaryotic "Final" (should not be different)

ggplot(data.frame(sum = sample_sums(Final)), aes(sum)) + 
  geom_histogram(color = "black", fill = "indianred", stat = "bin", binwidth = 500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") + 
  ylab("")
# mean, max and min of sample read counts
smin <- min(sample_sums(Final))
smean <- mean(sample_sums(Final))
smax <- max(sample_sums(Final))

#Graph with eukaryotic cells included "moth.sub" (we don't expect this to be different)
ggplot(data.frame(sum = sample_sums(moth.sub)), aes(sum)) +
  geom_histogram(color = "black", fill = "blue") +
  ggtitle("Distribution of Sample Sequencing Depth with Chloroplasts") +
  ylab("") + xlab("Read Counts")
#graphs did differ, meaning that there were some cholorplast/eukaryotic reads involved
#########################################################################################################################


#########################################################################################################################
###################################Stacked Bar plots of Relative Abundance###############################################
#starting into the graphing of a stacked bar plot of phyla >5% in sediment 
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
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "darkorchid3", "navy")

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
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 5%) \n") +
  ggtitle("Phylum Composition of Lake Michigan Sediment \n Bacterial Communities by Month Collected") 
#################################################################################################################3
#stacked bar blot of sediment phyla >2% total abundance 
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
################################################################################################################
#stacked bar plot of mussel samples with phyla >2% relative abundance
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
####################################################################################################
#stacked bar plot of phyla >3% of relative abundance in all viable samples
All_phylum3 <- scale_Final %>%
  subset_samples(Use == "Sample") %>%
  tax_glom(taxrank = "Phylum") %>%                  
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                        
  filter(Abundance > 0.03) %>%                        
  arrange(Phylum)

ggplot(All_phylum3, aes(x = Sed_or_Mus, y = Abundance, fill = Phylum)) + 
  #facet_grid(~.) +
  geom_bar(stat = "identity", position = "fill") +
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
  ylab("Relative Abundance (Phyla > 3%) \n") +
  ggtitle("Phylum Composition of Transect Sample \n Bacterial Communities by Sample Type")
####################
#Stacked barplot of Clostridium genera in all samples
c_bot_screen <- scale_Final %>%
  subset_samples(Use == "Sample") %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                        
  filter(Abundance > 0.03) %>%                        
  arrange(Family) 
c_bot_screen <- subset_taxa(Final, Family == "Clostridiaceae_1")
c_bot_genus <- subset_taxa(scale_Final, Genus == "Clostridium_sensu_stricto_1")

plot_bar(c_bot_screen, x="Sed_or_Mus", fill="Genus") +   
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
plot_bar(c_bot_screen, x="Type", fill="Genus") +   
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
plot_bar(c_bot_screen, x="Sediment_Type", fill="Genus") +   
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
plot_bar(c_bot_screen, x="Month_Collected", fill="Genus") +   
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
c_bot_genus %>%
  filter(abundance > 0) %>%
plot_bar(., x="New_Name", fill="Genus") +   
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") 

ggplot(c_bot_screen, aes(x = Sediment_Type, y = Abundance, fill = Genus)) + 
  #facet_grid(~.) +
  geom_bar(stat = "identity", position = "stack") +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 3%) \n") +
  ggtitle("Phylum Composition of Transect Sample \n Bacterial Communities by Sample Type")
##########################################################################################
#Stacked bar plot by plyha >3% in the mussel samples distinguished by size
MusselSizePhy3 <- scale_Final %>%
  subset_samples(Sed_or_Mus == "Mussel") %>%
  tax_glom(taxrank = "Phylum") %>%                  
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                        
  filter(Abundance > 0.03) %>%                        
  arrange(Phylum)

ggplot(MusselSizePhy3, aes(x = Mussel.Classification, y = Abundance, fill = Phylum)) + 
  # facet_grid(Station~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum.colors) +
  #scale_x_discrete(
  #breaks = c("7/8", "8/4", "9/2", "10/6"),
  #labels = c("Jul", "Aug", "Sep", "Oct"), 
  #drop = FALSE
  #  ) +
  # Remove x axis title
  xlab("Size of Mussel") + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 3%) \n") +
  ggtitle("Phylum Composition of Mussel \n Bacterial Communities by Mussel Size")
#########################################################################################################



######################################################################################################### 
########################## Graph of increasing number of Reads to determine scale factor ################
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
############################################## PCoA #############################################################
#scale reads to even depth/// This hhas already been done in the Formula section
Full.scale.format <-
  scale_Final %>%
  subset_samples(Use == "Sample") %>%
  scale_reads(round = "round")

#PCoA of all Viable Samples
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
##########################################################################
#PCoA of Sediment only
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
#####################################################################################
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
##########################################################################################
#PCoA by mussel density on collection site
DensityPCoA <-
  scale_Final %>%
  subset_samples(Mussel_Density %in%c("Very Few", "Few", "Few-Moderate", "Moderate", "Moderate-high", "Many", "High")) %>%
  scale_reads(round = "round")

Full_PCoA <- ordinate(physeq = DensityPCoA, 
                      method = "PCoA", 
                      distance = "bray")

#making plot
plot_ordination(
  physeq = DensityPCoA,
  ordination = Full_PCoA,
  color = "Mussel_Density",
  shape = "Sed_or_Mus",
  title = "PCoA of Lake Michigan Samples by Mussel Density"
) + 
  scale_color_manual(values = c("darkorchid3", "red", "#ffae19",
                                "#4daf4a", "#1919ff", "navy", "magenta"),
                     name = "Mussel Density"
  ) +
  scale_shape(name = "Sample Type") +
  geom_point(aes(color = Mussel_Density), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5)

##########################################################################################
#NMDS
set.seed(1)

# Ordinate
Final_nmds <- ordinate(
  physeq = scale_Final, 
  method = "NMDS", 
  distance = "bray"
)

plot_ordination(
  physeq = scale_Final,
  ordination = Final_nmds,
  color = "Type",
 # shape = "",
  title = "NMDS of Lake Micigan bacterial Communities"
) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta")
  ) +
  geom_point(aes(color = Type), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 
###############################################################################################################



###############################################################################################################
#################################Constrained Ordination #######################################################
#All Samples using scaled number of reads
Final_not_na <- scale_Final %>%
  subset_samples(!is.na(Type) &
                   !is.na(Lake))

bray_not_na <- phyloseq::distance(physeq = Final_not_na, method = "bray")


# CAP ordinate
cap_ord <- ordinate(
  physeq = Final_not_na, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~ Type + Lake)

# CAP plot
cap_plot <- plot_ordination(
  physeq = Final_not_na, 
  ordination = cap_ord, 
  color = "Type", 
  axes = c(1,2)
) + 
  aes(shape = Lake) + 
  geom_point(aes(colour = Type), alpha = 0.4, size = 4) + 
  geom_point(colour = "grey90", size = 1.5) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", 
                                "#1919ff", "darkorchid3", "magenta")
  )


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )+
  ggtitle("Constrained Ordination of Effect\n of Sample Type and Location on Dissimilarity")
##############################################################################################################



##############################################################################################################
##################################### Alpha Diversity ########################################################
#subset for depth of sediment duing transect cruise
transect_sed <- subset_samples(Final, Cruise == "Transect" & Type == "Sediment")

Final_sum_reps <- merge_samples(x = transect_sed, group = "Duplicate", fun = "sum") # Sum between replicate samples


Final_sum_reps <- prune_taxa(taxa_sums(Final_sum_reps) > 0, Final_sum_reps)
min(sample_sums(Final_sum_reps))
min(sample_sums(scale_Final))


sums_FWDB <- data.frame(rowSums(otu_table(Final_sum_reps)))
colnames(sums_FWDB) <- "Sample_TotalSeqs"
sums_FWDB$sample <- row.names(sums_FWDB)
sums_FWDB <- arrange(sums_FWDB, Sample_TotalSeqs)

ggplot(sums_FWDB, aes(x=reorder(sample, Sample_TotalSeqs), y = Sample_TotalSeqs)) + 
  ylab("Number of Sequences per Sample") +
  geom_bar(stat = "identity", colour="black",fill="cornflowerblue")  + xlab("Sample Name") + 
  ggtitle("Total Number of Sequences per Sample") + 
  theme(axis.text.x = element_text(colour = "black", size=6, angle=45, hjust = 1, vjust = 1))

View(sums_FWDB)

Final_sum_reps
rarefy_num <-  min(sample_sums(Final_sum_reps))-1
rarefy_num

nsamp = nsamples(Final_sum_reps)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(Final_sum_reps)
str(richness)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(Final_sum_reps)
str(evenness)


# It is always important to set a seed when you subsample so your result is replicable 
set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(Final_sum_reps, sample.size = rarefy_num, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even}


# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

alpha <- rbind(rich_stats, even_stats)
View(alpha)
#create new row called Duplicate that has the same info as SampleID so we can merge it with the metadata table
alpha$Duplicate <- alpha$SampleID

#merge the meta data with the new alpha data frame since factors were converted to numbers
alphadiv <- left_join(alpha, metadata, by = "Duplicate")
#disregard warning message 
View(alphadiv)

#plot the diversity
ggplot(alphadiv, aes(x = Month_Collected, y = mean, color = Depth_in_Meters, group = Depth_in_Meters)) +
  #geom_point(size = 2) + 
  geom_line(size = 0.8) +
  facet_wrap(~measure, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("#E96446", "#302F3D", "#87CEFA", "red")) +
  scale_x_discrete("Month_Collected", 
                   breaks = c("April", "May", "June", "July", "August", "September", "October"),
                   labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct"), 
                   drop = FALSE) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank())+
  ggtitle("Diversity of Transect Sediments \n By Depth and Month")
###################################################################################################################



####################################################################################################################
########FUNCTIONS####################################################################################################
matround <- function(x){trunc(x+0.5)}

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
