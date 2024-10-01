
#### Package setup ####

Sys.setenv(language = "EN")
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(phyloseqCompanion)
theme_set(theme_bw())

#### End ####

#### Set Directory ####
list.files("Data")

# Create directory to save results of the analysis and subdirectory for quality control

res.dir <- "Results_check"
qc.dir <- file.path(res.dir, "1.QC")
r.dir <- file.path(res.dir, "RData")
bar.dir <- file.path(res.dir, "6.Taxonomy/Bar_plots")
stats.dir <- file.path(res.dir, "6.Taxonomy/Stats")
lefse.dir <- file.path(res.dir, "7.Lefse")
#load(file.path(r.dir,"4.physeq.bracken.filtered.RData"))
dir.create(bar.dir, recursive = T)
dir.create(stats.dir, recursive = T)
dir.create(lefse.dir, recursive = T)

#### End ####

#### Import kraken-biom file into phyloseq object ####

biomfilename = "Data/S-table.biom"

data <- import_biom(biomfilename, parseFunction=parse_taxonomy_default)

# Change sample names
sample_names(data)
sample_names(data) = gsub("_\\.S|_lib.*|_merged.S", "", sample_names(data)) %>% gsub("NG-28767_", "P",.) %>%
  gsub("1_BA", "T1", .) %>% gsub("2_BA", "T2", .) %>% gsub("3_BA", "T3", .)
sample_names(data)
physeq = data

# Change rank names
rank_names(physeq)
colnames(tax_table(physeq)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "_Species")
rank_names(physeq)

#' Add column with Genus and Species to tax table called "Species",
#' change non-assigned taxa to NA and delete prefixes for easier plotting
tax <- as.data.frame(tax_table(physeq))
tax[tax$`_Species` == "s__", "_Species"] <- NA
tax[tax$Genus == "g__", "Genus"] <- NA
tax[tax$Family == "f__", "Family"] <- NA
tax[tax$Order == "o__", "Order"] <- NA
tax[tax$Class == "c__", "Class"] <- NA
tax[tax$Phylum == "p__", "Phylum"] <- NA
tax[tax$Kingdom == "k__", "Kingdom"] <- NA
tax <- dplyr::mutate(tax, Species = paste(tax[,"Genus"], tax[,"_Species"]))
tax[,"Species"] <- tax[,"Species"] %>% gsub("s__", "", .) %>% gsub("g__", "", .) %>% gsub(" NA", " sp.", .) %>% gsub("NA sp.", NA, .)
tax[,"Genus"] <- tax[,"Genus"] %>% gsub("g__", "", .)
tax[,"Family"] <- tax[,"Family"] %>% gsub("f__", "", .)
tax[,"Order"] <- tax[,"Order"] %>% gsub("o__", "", .)
tax[,"Class"] <- tax[,"Class"] %>% gsub("c__", "", .)
tax[,"Phylum"] <- tax[,"Phylum"] %>% gsub("p__", "", .)
tax$`_Species` <- NULL
tax <- as.matrix.data.frame(tax)
tax <- phyloseq::tax_table(tax)
phyloseq::tax_table(physeq) <- tax; rm(tax)

head(tax_table(physeq))

physeq

#### End ####

#### Save original count table and taxa table ####

count_tab = as(otu_table(physeq), "matrix")
write.table(count_tab, file.path(res.dir, "count_tab_bracken.tsv"), sep="\t", quote=F, col.names=NA)

tax_tab = as(tax_table(physeq), "matrix")
write.table(tax_tab, file.path(res.dir, "tax_tab_bracken.tsv"), sep="\t", quote=F, col.names=NA)

#### End ####

#### Positive control analysis ####

#### Remove positive control ####

to_remove <- c("PK1", "PK2")

physeq_patients <- prune_samples(!(sample_names(physeq) %in% to_remove), physeq)

# Remove ASVs that only appear in positive controls
physeq_patients = filter_taxa(physeq_patients, function(x) sum(x > 1) >= 1, TRUE)

physeq
physeq_patients
physeq=physeq_patients

#### End ####

#### Add samples metadata and match to count table####

sample_info_tab <- read.table("Data/zahn_metadata.txt",header=T, 
                              row.names = 1, check.names=F)

rownames(sample_info_tab) = gsub("_\\.kreport2|_merged|_lib.*", "", rownames(sample_info_tab)) %>% gsub("NG-28767_", "P",.) %>%
  gsub("1_BA", "T1", .) %>% gsub("2_BA", "T2", .) %>% gsub("3_BA", "T3", .)

# Examine consistancy in order between cts col names and coldata rownames 

sample_names(physeq)

all(rownames(sample_info_tab) %in% sample_names(physeq))
all(sample_names(physeq) %in% rownames(sample_info_tab))

# Add metadata
physeq_met <- merge_phyloseq(physeq, sample_data(sample_info_tab)) 
physeq = physeq_met

#### End ####

#### Remove Kingdoms other than Bacteria and NA Phyla ####

rank_names(physeq)

# Keep only Bacteria

table(tax_table(physeq)[, "Kingdom"], exclude = NULL)

physeq_B <- subset_taxa(physeq, !is.na(Kingdom) & Kingdom == "k__Bacteria")

table(tax_table(physeq_B)[, "Kingdom"], exclude = NULL)

physeq = physeq_B

table(tax_table(physeq)[, "Kingdom"], exclude = NULL)

# Remove NA Phyla

table(tax_table(physeq)[, "Phylum"], exclude = NULL)

physeq_oNA <- subset_taxa(physeq, !is.na(Phylum) & !Phylum %in% c("", "p__"))

physeq = physeq_oNA

table(tax_table(physeq)[, "Phylum"], exclude = NULL)

#### End ####

#### Inspect library sizes of bacteria ####

sample_data(physeq)$LibrarySizeBacteria <- sample_sums(physeq_B)

df <- as.data.frame(sample_data(physeq)) 
df <- df[order(df$LibrarySizeBacteria),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySizeBacteria, color=tooth, shape = implant)) + geom_point(size = 2)

rm(df)

#### End ####

#### Save count table and meta data ####

count_tab = as(otu_table(physeq), "matrix")
write.table(count_tab, "Results/count_tab_bracken_bacteria.tsv", sep="\t", quote=F, col.names=NA)

tax_tab = as(tax_table(physeq), "matrix")
write.table(tax_tab, "Results/tax_tab_bracken_bacteria.tsv", sep="\t", quote=F, col.names=NA)

meta_tab = as(sample_data(physeq), "matrix")
write.table(meta_tab, "Results/meta_tab_bracken.tsv", sep="\t", quote=F, col.names=NA)

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

save.image(file.path(r.dir,"3.physeq.bracken.original.RData"))

#### End ####

## Filter low Abundance Taxa and count table normalization ##

#### Define prevalence of each taxa # (in how many samples did each taxa appear at least once) ####
prev0 = apply(X = otu_table(physeq),
              MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(physeq),
                    tax_table(physeq))

# Save ASV Prevalence and Abundance table before filtering
write.table(prevdf, "Results/asv_prevdf_bracken.tsv", sep="\t", quote=F, col.names=NA)

# Plot Taxa prevalence v. total counts. Each point is a different taxa. 
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(physeq),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#### End ####

#### Remove taxa not seen more than 3 times in at least 5% of the samples #### 
# This protects against an OTU with small mean & trivially large C.V.
# Setting filter parameters :

countperphyla = 3
Samplepercentage = 0.05

physeq_filtered = filter_taxa(physeq, function(x) sum(x > countperphyla) > (Samplepercentage*length(x)), TRUE)
physeq
physeq_filtered
physeq = physeq_filtered

#### End ####

#### Normalize number of reads in each sample using median sequencing depth.####

total = median(sample_sums(physeq))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq_mednorm = transform_sample_counts(physeq, standf)

# Transform to relative abundance. Save as new object.
physeq_re = transform_sample_counts(physeq_mednorm, function(x){x / sum(x)})

#### End ####

#### Exploratory plots after filtering and normalization ####
# Check individual phylum Abundance
# Abundance value transformation function

plot_abundance = function(physeq, ylabn = "",
                          Facet = "Phylum",
                          Color = "Class", 
                          x = "implant"){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq,
         mapping = aes_string(x = x, y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + ylab(ylabn) +
    scale_y_log10()
}


#plot the abundance values before and after transformation

pl_ab_original  = plot_abundance(physeq,"Original Abundances", Color="Phylum")
pl_ab_original_norm  =plot_abundance(physeq_mednorm,"Normalized to sequencing depth Abundances", Color = "Phylum")
pl_ab_original_norm_re  =plot_abundance(physeq_re,"Normalized Relative Abundances", Color = "Phylum")

grid.arrange(pl_ab_original, pl_ab_original_norm, pl_ab_original_norm_re)

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

save.image(file.path(r.dir,"4.physeq.bracken.filtered.RData"))

#### End ####

################################################################################
########################### Taxonomy abundance plots ###########################
################################################################################

#### Create variable to merge samples from the same condition ####

sample_data(physeq)$group <- paste(sample_data(physeq)$implant, sample_data(physeq)$time, sample_data(physeq)$tooth, sep = "_")
sample_data(physeq_re)$group <- paste(sample_data(physeq_re)$implant, sample_data(physeq_re)$time, sample_data(physeq_re)$tooth, sep = "_")

#### Plot abundance in Fibula or DCIA samples ####

# Function to create phyloseq object according to condition

phy_cond <- function(phy, condition){
  
  # Merging samples according to condition
  ps.condition <- merge_samples(phy, condition)  
  ps.condition_per <- transform_sample_counts(ps.condition, function(OTU) OTU/sum(OTU)) #2nd transformation to make it again in percentage
  sample_data(ps.condition_per)[,condition] <- rownames(sample_data(ps.condition_per))
  
  return(ps.condition_per)
  
}

# Function to plot the top n taxa in each type of implant

top <- function(phy, taxa, ntop, implant){
  sum = tapply(taxa_sums(phy), tax_table(phy)[, taxa], sum, na.rm=TRUE)
  top = names(sort(sum, TRUE))[1:ntop]

  phy_top = prune_taxa((tax_table(phy)[, taxa] %in% top), phy)
  
  # New facet label names for tooth variable
  tooth.labs <- c("Ceramic", "Tooth")
  names(tooth.labs) <- c("ceramic", "tooth")
  
  plot <- plot_bar(phy_top, taxa, fill=taxa,
           title = paste("Relative abundance from top", ntop, taxa, "in", implant, sep = " ")) +
    geom_bar(stat = "Identity", position = "stack") +
    facet_grid(time~tooth, labeller = labeller(tooth = tooth.labs)) +
    ylab("Relative abundance") + xlab(NULL) +
    theme(text = element_text(size = 22, color = "black", face = "bold"),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(angle = -30, size = 18, face = "bold"))
  return(plot)
}

### Fibula
phy_Fib <- prune_samples((sample_data(physeq_re)$implant == "Fibula"),physeq_re)
ps.group.Fib <- phy_cond(phy_Fib, "group")

# As metadata values have been converted to NA, extract metadata from the group name
sample_data(ps.group.Fib)$tooth <- gsub(".*_", "", sample_data(ps.group.Fib)$group)
sample_data(ps.group.Fib)$time <- gsub("DCIA_", "", sample_data(ps.group.Fib)$group) %>% gsub("Fibula_", "", .) %>%
  gsub("_tooth", "", .) %>% gsub("_ceramic", "", .)

# Plot top 10 taxa in fibula implants
top(ps.group.Fib, "Phylum", 5, "Fibula")

F5B <- top(ps.group.Fib, "Genus", 10, "Fibula")

ggsave(file.path(bar.dir,"Figure_5B_Top_Genus_Fibula.tiff"), dpi = 300, height = 12, width = 14)
ggsave(file.path(bar.dir,"Figure_5B_Top_Genus_Fibula.svg"), dpi = 300, height = 12, width = 14)
ggsave(file.path(bar.dir,"Figure_5B_Top_Genus_Fibula.png"), dpi = 600, height = 12, width = 14)

F6B <- top(ps.group.Fib, "Species", 10, "Fibula")

ggsave(file.path(bar.dir,"Figure_6B_Top_Species_Fibula.tiff"), dpi = 300, height = 12, width = 14)
ggsave(file.path(bar.dir,"Figure_6B_Top_Species_Fibula.svg"), dpi = 300, height = 12, width = 14)
ggsave(file.path(bar.dir,"Figure_6B_Top_Species_Fibula.png"), dpi = 600, height = 12, width = 14)

### DCIA
phy_DCIA <- prune_samples((sample_data(physeq_re)$implant == "DCIA"),physeq_re)
ps.group.DCIA <- phy_cond(phy_DCIA, "group")

# As metadata values have been converted to NA, extract metadata from the group name
sample_data(ps.group.DCIA)$tooth <- gsub(".*_", "", sample_data(ps.group.DCIA)$group)
sample_data(ps.group.DCIA)$time <- gsub("DCIA_", "", sample_data(ps.group.DCIA)$group) %>% gsub("Fibula_", "", .) %>%
  gsub("_tooth", "", .) %>% gsub("_ceramic", "", .)

# Plot top 10 taxa in fibula implants
top(ps.group.DCIA, "Phylum", 5, "DCIA")

F5A <- top(ps.group.DCIA, "Genus", 10, "DCIA")

ggsave(file.path(bar.dir,"Figure_5A_Top_Genus_DCIA.tiff"), dpi = 300, height = 12, width = 14)
ggsave(file.path(bar.dir,"Figure_5A_Top_Genus_DCIA.svg"), dpi = 300, height = 12, width = 14)
ggsave(file.path(bar.dir,"Figure_5A_Top_Genus_DCIA.png"), dpi = 600, height = 12, width = 14)

F6A <- top(ps.group.DCIA, "Species", 10, "DCIA")

ggsave(file.path(bar.dir,"Figure_6A_Top_Species_DCIA.tiff"), dpi = 300, height = 12, width = 14)
ggsave(file.path(bar.dir,"Figure_6A_Top_Species_DCIA.svg"), dpi = 300, height = 12, width = 14)
ggsave(file.path(bar.dir,"Figure_6A_Top_Species_DCIA.png"), dpi = 600, height = 12, width = 14)

# Combine plots
library(cowplot)

plot_grid(F5A,F5B, labels= "AUTO", label_size = 22)

ggsave(file.path(bar.dir,"Figure_5.tiff"), dpi = 300, height = 12, width = 26)
ggsave(file.path(bar.dir,"Figure_5.svg"), dpi = 300, height = 12, width = 28)
ggsave(file.path(bar.dir,"Figure_5.png"), dpi = 600, height = 12, width = 28)

plot_grid(F6A,F6B, labels= "AUTO", label_size = 22)

ggsave(file.path(bar.dir,"Figure_6.tiff"), dpi = 300, height = 12, width = 26)
ggsave(file.path(bar.dir,"Figure_6.svg"), dpi = 300, height = 12, width = 28)
ggsave(file.path(bar.dir,"Figure_6.png"), dpi = 600, height = 12, width = 28)

#### Bubble plot of specific taxa ####

# Select specific taxa to plot

taxa <- c("Treponema medium", "Treponema denticola", "Tannerella forsythia",
        "Selenomonas sputigena", "Prevotella intermedia", "Porphyromonas gingivalis",
        "Porphyromonas endodontalis", "Parvimonas micra", "Fusobacterium nucleatum",
        "Fretibacterium fastidiosum", "Filifactor alocis", "Eikenella corrodens",
        "Aggregatibacter actinomycetemcomitans", "Actinomyces oris", "Actinomyces naeslundii")

# Create sample variable containing all variables to plot
phy_bubble_group <- phy_cond(physeq_re, "group")

# Function for creating long data frames for plotting
long_df <- function(taxon){
  
  # Get taxa to plot that are present in the phyloseq object
  sel_tax <- unique(tax_table(phy_bubble_group)[,taxon])
  
  # Agglomerate to taxonomy level
  if(taxon != "Species"){
  phy_bubble_glom <- tax_glom(phy_bubble_group, taxrank = taxon)
  }else{phy_bubble_glom <- phy_bubble_group}

  # Get relative abundance matrix for each group and specific taxon level
  taxa_names(phy_bubble_glom) <- tax_table(phy_bubble_glom)[,taxon]
  tax_sums <- as.data.frame(t(otu_table(phy_bubble_glom)))

  # Make long data frame with metadata
  longdf <- tax_sums %>% tibble::rownames_to_column("tax") %>%
    tidyr::pivot_longer(cols = -tax, names_to = "groups", values_to = "abund") %>%
    mutate(
      implant = str_replace(groups, "_.*", ""),
      tooth = str_replace(groups, ".*_", ""),
      time = str_replace(str_replace(groups, "DCIA_|Fibula_", ""), "_.*", ""))

  return(longdf)
}

# Create long data frame for Phylum
longdf_ph <- long_df("Phylum")

# Create long data frame for Species
longdf_sp <- long_df("Species")
longdf_sps <- longdf_sp[longdf_sp$tax %in% taxa,]

# Bubble plots 

bubble_plot <- function(longdf){
  ggplot(longdf,aes(time,tax)) +
    geom_point(aes(size=abund, fill=tax, color=tax),shape=21) +
    facet_wrap(implant~tooth, nrow = 1, scales = "free_x", strip.position = "bottom") +
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          panel.spacing = unit(0, "lines"),
          text = element_text(size = 22, color = "black", face = "bold"),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.text.y = element_text(size = 15, color = "black"),
          axis.text.x = element_text(size = 18, color = "black"))+
    ylab("") +
    xlab(NULL) +
    guides(fill = "none", color = "none") +
    scale_size(name = "Relative\nabundance")
}

#bubble_plot(longdf_ph)
bubble_plot(longdf_sps)

ggsave(file.path(bar.dir,"Figure_9.tiff"), dpi = 300, height = 8, width = 14)
ggsave(file.path(bar.dir,"Figure_9.svg"), dpi = 300, height = 8, width = 14)
ggsave(file.path(bar.dir,"Figure_9.png"), dpi = 600, height = 8, width = 14)


################################################################################
############################## Create LEfSe input ##############################
################################################################################

#### Save as input for LEfSe #### 

# Select phyloseq object to save
phy2lefse <- physeq_re         

# Changing taxa names so that they are easier to recognize in LEfSe
tax <- as.data.frame(tax_table(phy2lefse))
tax <- tax[,c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")]
tax$Species <- paste0("s__", gsub(" ", "_", tax$Species))
tax$Genus <- paste0("g__", tax$Genus)
tax$Family <- paste0("f__", tax$Family)
tax$Order <- paste0("o__", tax$Order)
tax$Class <- paste0("c__", tax$Class)
tax$Phylum <- paste0("p__", tax$Phylum)
tax <- as.matrix.data.frame(tax)
tax <- phyloseq::tax_table(tax)
phyloseq::tax_table(phy2lefse) <- tax; rm(tax)

head(tax_table(phy2lefse))

# Save as input for LEfSe
phyloseq2lefse(phy2lefse, c("implant", "tooth"), file.name = file.path(lefse.dir, "physeq_re_lfse_imp_too_p-s.txt"),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)

# Save as input for LEfSe only with DCIA flaps
phy2lefse_DCIA <- subset_samples(phy2lefse, implant == "DCIA")
phy2lefse_DCIA <- prune_taxa(taxa_sums(phy2lefse_DCIA) > 0, phy2lefse_DCIA)
phyloseq2lefse(phy2lefse_DCIA, c("tooth"), file.name = file.path(lefse.dir, "physeq_re_lfse_too_p-s_DCIA_reordered.txt"),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)

# Save as input for LEfSe only with Fibula flaps
phy2lefse_Fibula <- subset_samples(phy2lefse, implant == "Fibula")
phy2lefse_Fibula <- prune_taxa(taxa_sums(phy2lefse_Fibula) > 0, phy2lefse_Fibula)
phyloseq2lefse(phy2lefse_Fibula, c("tooth"), file.name = file.path(lefse.dir, "physeq_re_lfse_too_p-s_Fibula_reordered.txt"),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)

#### End ####