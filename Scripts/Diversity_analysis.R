
#### Package setup ####

Sys.setenv(language = "EN")
library(vegan)
library(phyloseq)
library(dplyr)
library(openxlsx)
library(vegan) 
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(qqplotr)
library(ggpattern)
library(DescTools)
theme_set(theme_bw())

#### End ####

#### Set Directory ####
list.files("Data")

# Create directory to save results of the analysis and subdirectory for quality control

res.dir <- "Results_check"
qc.dir <- file.path(res.dir, "1.QC")
r.dir <- file.path(res.dir, "RData")
a.stats <- file.path(res.dir, "2.Alpha_stats")
a.plots <- file.path(res.dir, "3.Alpha_plots")
b.stats <- file.path(res.dir, "4.Beta_stats")
b.plots <- file.path(res.dir, "5.Beta_plots")

dir.create(qc.dir, recursive = T)
dir.create(r.dir, recursive = T)
dir.create(a.stats, recursive = T)
dir.create(a.plots, recursive = T)
dir.create(b.stats, recursive = T)
dir.create(b.plots, recursive = T)

#### End ####

#### Import kraken-biom file into phyloseq object ####

biomfilename = "Data/A-table-kraken.biom"

data <- import_biom(biomfilename, parseFunction=parse_taxonomy_default)

# Merge 2_1_K # Delete from final script!!!
otu_table(data)[,"NG-28767_2_1_BA_K_lib545751_7734_2_.kreport2"] <- otu_table(data)[,"NG-28767_2_1_BA_K_lib545751_7734_2_.kreport2"] +
  otu_table(data)[,"NG-28767_2_1_BA_K_lib545751_7766_1_.kreport2"]

data <- prune_samples((sample_names(data) != "NG-28767_2_1_BA_K_lib545751_7766_1_.kreport2"), data)

# Change sample names
sample_names(data)
sample_names(data) = gsub("_\\.kreport2|_lib.*", "", sample_names(data)) %>% gsub("NG-28767_", "P",.) %>%
  gsub("1_BA", "T1", .) %>% gsub("2_BA", "T2", .) %>% gsub("3_BA", "T3", .)
sample_names(data)
physeq = data

# Change rank names
rank_names(physeq)
colnames(tax_table(physeq)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "_Species")
rank_names(physeq)

#### End ####

#### Save original count table and taxa table ####

count_tab = as(otu_table(physeq), "matrix")
write.table(count_tab, file.path(res.dir, "count_tab_kraken.tsv"), sep="\t", quote=F, col.names=NA)

tax_tab = as(tax_table(physeq), "matrix")
write.table(tax_tab, file.path(res.dir, "tax_tab_kraken.tsv"), sep="\t", quote=F, col.names=NA)

#### End ####

#### Rarefaction curves ####
count_tab <- as.data.frame(otu_table(data))

col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")

raremax <- min(rowSums(t(count_tab)))
raremax

# Without Mock Community
count_tab_nomock <- count_tab[, -grep("PK", colnames(count_tab))]
raremax_nomock <- min(rowSums(t(count_tab_nomock)))
raremax_nomock

# Option 1
rarecurve(t(count_tab_nomock), step = 100,  col = col, lwd=2, lty = lty, ylab = "ASVs", label = F)
abline(v=(raremax_nomock))

# Option 2
rarecurve(t(raremax_nomock), step = 100, sample = raremax, col = col, lwd=2, lty = lty, ylab = "ASVs", label = F)

#### End ####

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
physeq_met <- merge_phyloseq(physeq_patients, sample_data(sample_info_tab)) 
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
write.table(count_tab, "Results/count_tab_kraken_bacteria.tsv", sep="\t", quote=F, col.names=NA)

tax_tab = as(tax_table(physeq), "matrix")
write.table(tax_tab, "Results/tax_tab_kraken_bacteria.tsv", sep="\t", quote=F, col.names=NA)

meta_tab = as(sample_data(physeq), "matrix")
write.table(meta_tab, "Results/meta_tab.tsv", sep="\t", quote=F, col.names=NA)

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

save.image(file.path(r.dir,"1.physeq.kraken.original.RData"))

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
write.table(prevdf, "Results/asv_prevdf_kraken.tsv", sep="\t", quote=F, col.names=NA)

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

save.image(file.path(r.dir,"2.physeq.kraken.filtered.RData"))

#### End ####

################################################################################
############################### Alpha diversity ################################
################################################################################

#### Create richness matrix ####
# Has to be counts before low abundance filtering, so we take physeq_oNA
rich = estimate_richness(physeq_oNA, measures = c("Observed", "Shannon","InvSimpson"))

# Merge with sample_data for plotting
rich <- merge(rich, meta_tab, by = "row.names")
rich <- rich %>% mutate(tooth=factor(tooth, levels=c("tooth","ceramic")))

# Save table with alpha diversity

write.xlsx(rich, file.path(a.stats, "richness.xlsx"))

#### End ####

#### Check distribution of alpha diversity indices ####

# Produce descriptive statistics
alpha <- list()
for(i in c("Observed", "Shannon", "InvSimpson")){
  alpha[[i]] <- cbind("alpha.index" = i,
                      "n" = length(rich[,i]), 
                      "mean" = mean(rich[,i], na.rm = TRUE), 
                      "sd "= sd(rich[,i], na.rm = TRUE),
                      "stderr" = sd(rich[,i], na.rm = TRUE)/sqrt(length(rich[,i])),
                      "LCL" = mean(rich[,i], na.rm = TRUE) - qt(1 - (0.05 / 2), length(rich[,i]) - 1) * sd(rich[,i], na.rm = TRUE)/sqrt(length(rich[,i])),
                      "UCL" = mean(rich[,i], na.rm = TRUE) + qt(1 - (0.05 / 2), length(rich[,i]) - 1) * sd(rich[,i], na.rm = TRUE)/sqrt(length(rich[,i])),
                      "median" = median(rich[,i], na.rm = TRUE),
                      "min" = min(rich[,i], na.rm = TRUE), 
                      "max" = max(rich[,i], na.rm = TRUE),
                      "IQR" = IQR(rich[,i], na.rm = TRUE),
                      "LCLmed" = MedianCI(rich[,i], na.rm=TRUE)[2],
                      "UCLmed" = MedianCI(rich[,i], na.rm=TRUE)[3],
                      "shapiro.W.Stat" = shapiro.test(rich[,i])$statistic,
                      "shapiro.p.value" = shapiro.test(rich[,i])$p.value)
  
}
alpha <- as.data.frame(do.call(rbind, alpha))

# View result
alpha

# View only Shapiro test
alpha[,c("alpha.index", "shapiro.W.Stat", "shapiro.p.value")] 

# Save
write.xlsx(alpha, file.path(a.stats, "descriptive_stats_shapiro.xlsx"))

# Create column for facets
rich$Sh <- "Shannon"
rich$Ob <- "Observed"
rich$Si <- "InvSimpson"

# Create histograms 
H_Ob <- ggplot(data = rich, mapping = aes(x = Observed)) +
  geom_histogram(aes(y=..density..)) +
  geom_density(alpha=.2, fill="#FF6666") +
  facet_grid(. ~ Ob) 

H_Sh <- ggplot(data = rich, mapping = aes(x = Shannon)) +
  geom_histogram(aes(y=..density..)) +
  geom_density(alpha=.2, fill="#FF6666") +
  facet_grid(. ~ Sh) 

H_Si <- ggplot(data = rich, mapping = aes(x = InvSimpson)) +
  geom_histogram(aes(y=..density..)) +
  geom_density(alpha=.2, fill="#FF6666") +
  facet_grid(. ~ Si) 

grid.arrange(H_Ob, H_Sh, H_Si, nrow = 1)

# Create QQ plots
QQ_Ob <- ggplot(data = rich, mapping = aes(sample = Observed)) +
  stat_qq_band(alpha=0.5, conf=0.95, qtype=1, bandType = "boot") +
  stat_qq_line(identity=TRUE) +
  stat_qq_point(col="black") +
  facet_grid(. ~ Ob) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") + theme_bw()

QQ_Sh <- ggplot(data = rich, mapping = aes(sample = Shannon)) +
  stat_qq_band(alpha=0.5, conf=0.95, qtype=1, bandType = "boot") +
  stat_qq_line(identity=TRUE) +
  stat_qq_point(col="black") +
  facet_grid(. ~ Sh) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") + theme_bw()

QQ_Si <- ggplot(data = rich, mapping = aes(sample = InvSimpson)) +
  stat_qq_band(alpha=0.5, conf=0.95, qtype=1, bandType = "boot") +
  stat_qq_line(identity=TRUE) +
  stat_qq_point(col="black") +
  facet_grid(. ~ Si) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") + theme_bw()

grid.arrange(QQ_Ob, QQ_Sh, QQ_Si, nrow = 1)

#' Observed and Shannon are normally distributed and InvSimpson has a nonparametric distribution
#' as indicated by the Shapiro test, histograms and qq plots

#### End ####

#### Alpha diversity statistics ####
#' Kruskal-Wallis and Mann-Whitney are used for the non-normally distributed variables
#' One Way Anova and t-test are used for the normally distributed variables
#### Krustal-Wallis ####

## Comparison of all timepoints including all samples ##
kruskal <- kruskal.test(rich$InvSimpson ~ rich$time)
eta_squared <- (kruskal$statistic - 3 + 1)/(nrow(rich) - 3)
f <- sqrt(eta_squared/(1-eta_squared))

kw <- c(kruskal$method, "InvSimpson", "time", kruskal$statistic, kruskal$p.value, f)
names(kw) <- c("Test", "Variable1", "Variable2", "Chi2", "p-value", "Effect size")

write.xlsx(t(as.data.frame(kw)), file = file.path(a.stats, "Kruskal.wallis_time_all_samples.xlsx"))

## Comparison of all timepoints per group of implant and tooth ##
kw <- list()

for (implant in c("Fibula", "DCIA")) {
  for (tooth in c("tooth", "ceramic")) {
    dat <- rich[rich$tooth == tooth & rich$implant == implant,]
    kruskal <- kruskal.test(dat$InvSimpson ~ dat$time)
    eta_squared <- (kruskal$statistic - 3 + 1)/(nrow(dat) - 3)
    f <- sqrt(eta_squared/(1-eta_squared))
    
    tab <- c(kruskal$method, implant, tooth, "InvSimpson", "time", kruskal$statistic, kruskal$p.value, f)
    kw[[paste0(implant, "_", tooth)]] <- tab
    
  }}


kw <- as.data.frame(do.call(rbind, kw))
colnames(kw) <- c("Test", "Implant", "Tooth", "Variable1", "Variable2", "Chi2", "p-value", "Effect size")

write.xlsx(kw, file = file.path(a.stats, "Kruskal.wallis_time_per_group.xlsx"))

#### Mann-Whitney-U-Test #### 

## Comparison of tooth and implant including all samples ##
mwut <- list()
for (i in c("tooth" , "implant")){
  if(nrow(rich) <51){exact = T}else{exact = F}
  U_test <- wilcox.test(rich[,"InvSimpson"] ~ rich[,i], exact = exact)
  z <- abs(qnorm(U_test$p.value/2))
  r <- z/sqrt(nrow(rich))           #Effect size calculation
  
  tab <- c(U_test$method, "InvSimpson", i , U_test$statistic, U_test$p.value, r)
  mwut[[i]] <- tab
}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

write.xlsx(mwut, file = file.path(a.stats, "Mann.Whitney_tooth_implant_all_samples.xlsx"))


## Comparison of tooth and implant per timepoint ##
mwut <- list()

for (t in c("T1", "T2", "T3")) {
  for (i in c("implant", "tooth")) {
    U_test <- wilcox.test(rich[rich$time == t, "InvSimpson"] ~ rich[rich$time == t, i], data = rich, exact = T)
    z <- abs(qnorm(U_test$p.value/2))
    r <- z/sqrt(nrow(rich))
    tab <- c(U_test$method, t, "InvSimpson", i, U_test$statistic, U_test$p.value, r)
    mwut[[paste0(t, "_", i)]] <- tab
    
  }
}
mwut <- do.call(rbind, mwut)
colnames(mwut) <- c("Test", "Time point", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

write.xlsx(mwut, file = file.path(a.stats, "Mann.Whitney_tooth_implant_per_timepoint.xlsx"))


## Comparison of tooth and implant per group of the other varible (tooth or implant) and per timepoint ##
mwut <- list()

for (t in c("T1", "T2", "T3")) {
  for (i in c("implant", "tooth")) {
    if(i == "implant"){v1 = "ceramic" ; v2 = "tooth" ; i2 = "tooth"}
    else if(i == "tooth"){v1 = "Fibula" ; v2 = "DCIA" ; i2 = "implant"}
    for(v in c(v1, v2)){
      dat <- rich[rich$time == t & rich[, i2] == v,]
      if(nrow(dat) <51){exact = T}else{exact = F}
      U_test <- wilcox.test(dat[, "InvSimpson"] ~ dat[, i], data = rich, exact = exact)
      z <- abs(qnorm(U_test$p.value/2))
      r <- z/sqrt(nrow(dat))
      
      tab <- c(U_test$method, t, v, "InvSimpson", i, U_test$statistic, U_test$p.value, r)
      mwut[[paste0(t, "_", i, "_", v)]] <- tab
      
    }}}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Time point", "Group", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")
mwut <- mwut[order(mwut$Group),]
mwut$p.adj.BH <- p.adjust(mwut$`p-value`, method = "BH")

write.xlsx(mwut, file = file.path(a.stats, "Mann.Whitney_tooth_implant_per_group_and_timepoint.xlsx"))


## Comparison of two timepoints including all samples ##
mwut <- list()
time <- levels(as.factor(rich$time))

for (t in time) {
  dat <- rich[rich$time != t,]
  if(nrow(dat) <51){exact = T}else{exact = F}
  U_test <- wilcox.test(dat[, "InvSimpson"] ~ dat[, "time"], data = dat, exact = exact)
  z <- abs(qnorm(U_test$p.value/2))
  r <- z/sqrt(nrow(dat))    #Effect size calculation
  
  tab <- c(U_test$method, time[time != t][1], time[time != t][2], "InvSimpson", "time" , U_test$statistic, U_test$p.value, r)
  mwut[[t]] <- tab
  
}
mwut2 <- as.data.frame(do.call(rbind, mwut))
colnames(mwut2) <- c("Test", "Timepoint A", "Timepoint B", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")
mwut2$p.adj <- p.adjust(mwut2$`p-value`, method = "BH")

write.xlsx(kw, file = file.path(a.stats, "Mann.Whitney_time_all_samples.xlsx"))


## Comparison of two timepoints per group of implant and tooth ##
mwut <- list()

time <- levels(as.factor(rich$time))
for (t in time) {
  for (implant in c("Fibula", "DCIA")) {
    for (tooth in c("tooth", "ceramic")) {
      dat <- rich[rich$time != t & rich$tooth == tooth & rich$implant == implant,]
      if(nrow(dat) <51){exact = T}else{exact = F}
      U_test <- wilcox.test(dat[, "InvSimpson"] ~ dat[, "time"], data = dat, exact = exact)
      z <- abs(qnorm(U_test$p.value/2))
      r <- z/sqrt(nrow(dat))    #Effect size calculation
      
      tab <- c(U_test$method, time[time != t], implant, tooth, "InvSimpson", "time", U_test$statistic, U_test$p.value, r)
      mwut[[paste0(t, "_", implant, "_", tooth)]] <- tab
      
    }}}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Time1", "Time2", "Implant", "Tooth", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")
mwut <- mwut %>% arrange(Time1) %>% arrange(Time2) %>% as.data.frame(.)
mwut$p.adj.BH <- p.adjust(mwut$`p-value`, method = "BH")

write.xlsx(mwut, file = file.path(a.stats, "Mann.Whitney_time_per_group.xlsx"))


#### One Way Anova test ####

## Comparison of all timepoints including all samples ##
aov <- list()
for(n in c("Shannon", "Observed")){
  aov_test <- aov(rich[,n] ~ rich[,"time"], data = rich)
  aov_sum <- summary(aov_test)
  tab <- c(n, "time", aov_sum[[1]]$`Sum Sq`, aov_sum[[1]]$`Mean Sq`, aov_sum[[1]]$`F value`, aov_sum[[1]]$`Pr(>F)`)
  aov[[n]] <- tab
}
aov <- as.data.frame(do.call(rbind, aov)) ; aov <- aov[ , colSums(is.na(aov)) == 0]
colnames(aov) <- c("Variable1", "Variable2", "Sum Sq [i]", "Sum Sq Residuals", "Mean Sq [i]", "Mean Sq Residuals", "F value", "Pr(>F)")

write.xlsx(as.data.frame(aov), file = file.path(a.stats, "One.Way.Anova_time_all_samples.xlsx"))


## Comparison of all timepoints per group of implant and tooth ##
aov <- list()

for (implant in c("Fibula", "DCIA")){
  for (tooth in c("tooth", "ceramic")){
    for(n in c("Shannon", "Observed")){
      dat <- rich[rich$tooth == tooth & rich$implant == implant,]
      aov_test <- aov(dat[,n] ~ time, data = dat)
      aov_sum <- summary(aov_test)
      tab <- c(n , "time", implant, tooth, aov_sum[[1]]$`Sum Sq`, aov_sum[[1]]$`Mean Sq`, aov_sum[[1]]$`F value`, aov_sum[[1]]$`Pr(>F)`)
      aov[[paste0(implant, "_", tooth, "_", n)]] <- tab
      
    }}}

aov <- as.data.frame(do.call(rbind, aov)) ; aov <- aov[ , colSums(is.na(aov)) == 0]
colnames(aov) <- c("Variable1", "Variable2", "Implant", "Tooth", "Sum Sq [i]", "Sum Sq Residuals", "Mean Sq [i]", "Mean Sq Residuals", "F value", "Pr(>F)")

write.xlsx(as.data.frame(aov), file = file.path(a.stats, "One.Way.Anova_time_per_group.xlsx"))


#### t-test ####

## Comparison of tooth and implant including all samples ##
tt <- list()
for(i in c("implant", "tooth")){
  for(n in c("Shannon", "Observed")){
    t_test <- t.test(rich[,n] ~ rich[,i])
    tab <- c(t_test$method, n, i, t_test$statistic, t_test$parameter, t_test$estimate, t_test$p.value)
    tt[[paste0(i, "_", n)]] <- tab
  }}
tt <- as.data.frame(do.call(rbind, tt)) 
colnames(tt) <- c("Method", "Variable1", "Variable2", "t", "df", "mean in group 1", "mean in group 2", "p-value")

write.xlsx(tt, file = file.path(a.stats, "t-test_tooth_implant_all_samples.xlsx"))

## Comparison of tooth and implant per timepoint ##
tt <- list()
for (t in c("T1", "T2", "T3")) {
  for (i in c("implant", "tooth")) {
    for(n in c("Shannon", "Observed")){
      t_test <- t.test(rich[rich$time == t, "InvSimpson"] ~ rich[rich$time == t, i])
      tab <- c(t_test$method, t, n, i, t_test$statistic, t_test$parameter, t_test$estimate, t_test$p.value)
      tt[[paste0(t, "_", i, "_", n)]] <- tab
      
    }
  }
}
tt <- as.data.frame(do.call(rbind, tt)) 
colnames(tt) <- c("Method", "Time point", "Variable1", "Variable2", "t", "df", "mean in group 1", "mean in group 2", "p-value")
write.xlsx(tt, file = file.path(a.stats, "t-test_tooth_implant_per_timepoint.xlsx"))

## Comparison of tooth and implant per group of the other varible (tooth or implant) and per timepoint ##
tt <- list()
for (t in c("T1", "T2", "T3")) {
  for (i in c("implant", "tooth")) {
    if(i == "implant"){v1 = "ceramic" ; v2 = "tooth" ; i2 = "tooth"}
    else if(i == "tooth"){v1 = "Fibula" ; v2 = "DCIA" ; i2 = "implant"}
    for(v in c(v1, v2)){
      for(n in c("Shannon", "Observed")){
        dat <- rich[rich$time == t & rich[, i2] == v,]
        t_test <- t.test(dat[,n] ~ dat[,i])
        tab <- c(t_test$method, t, v, n, i, t_test$statistic, t_test$parameter, t_test$estimate, t_test$p.value)
        tt[[paste0(n, "_", t, "_", i, "_", v)]] <- tab
        
      }}}}

tt <- as.data.frame(do.call(rbind, tt)) 
colnames(tt) <- c("Method", "Time point", "Group", "Variable1", "Variable2", "t", "df", "mean in group 1", "mean in group 2", "p-value")
write.xlsx(tt, file = file.path(a.stats, "t-test_tooth_implant_per_group_and_timepoint.xlsx"))


#### End ####

#### Alpha diversity plot ####

# Create plot with the right statistics

# Observed diversity
P_Ob <- ggplot(rich, aes(x = time, y = Observed, shape = implant, color = tooth, Ob)) 

P_Ob.1 <- P_Ob + geom_boxplot_pattern(aes(pattern = implant, color = tooth, fill = tooth, pattern_fill = "white", alpha = 0.3),
                                      pattern_spacing = 0.02, outlier.shape=NA) +
  scale_pattern_manual(values = c("stripe","none"), name= "Implant") +
  guides(pattern_fill = F, alpha = F,
         color=guide_legend(title = "Tooth"),
         fill=guide_legend(title = "Tooth"),
         shape=guide_legend(title = "Implant")) +
  geom_point(position=position_jitterdodge(jitter.width = 0),alpha=0.6, size = 3) +
  facet_grid(. ~ Ob) +
  ylab(NULL) + xlab("") +
  theme(text = element_text(size = 22, color = "black", face = "bold"))

# Shannon diversity
P_Sh <- ggplot(rich, aes(x = time, y = Shannon, shape = implant, color = tooth, Sh)) 

P_Sh.1 <- P_Sh + geom_boxplot_pattern(aes(pattern = implant, color = tooth, fill = tooth, pattern_fill = "white", alpha = 0.3),
                                      pattern_spacing = 0.02, outlier.shape=NA) +
  scale_pattern_manual(values = c("stripe","none"), name= "Implant") +
  guides(pattern_fill = F, alpha = F,
         color=guide_legend(title = "Tooth"),
         fill=guide_legend(title = "Tooth"),
         shape=guide_legend(title = "Implant")) +
  geom_point(position=position_jitterdodge(jitter.width = 0),alpha=0.6, size = 3) +
  facet_grid(. ~ Sh) +
  ylab(NULL) + xlab("Time") +
  theme(text = element_text(size = 22, color = "black", face = "bold"))

# Inversed Simpson
P_Si <- ggplot(rich, aes(x = time, y = InvSimpson, shape = implant, color = tooth, Si)) 

P_Si.1 <- P_Si + geom_boxplot_pattern(aes(pattern = implant, color = tooth, fill = tooth, pattern_fill = "white", alpha = 0.3),
                                      pattern_spacing = 0.02,outlier.shape=NA) +
  scale_pattern_manual(values = c("stripe","none"), name= "Implant") +
  guides(pattern_fill = F, alpha = F,
         color=guide_legend(title = "Tooth"),
         fill=guide_legend(title = "Tooth"),
         shape=guide_legend(title = "Implant")) +
  geom_point(position=position_jitterdodge(jitter.width = 0),alpha=0.6, size = 3) +
  facet_grid(. ~ Si) +
  ylab(NULL) + xlab("") +
  theme(text = element_text(size = 22, color = "black", face = "bold")) 

# Combine all plots
P_alpha <- grid.arrange(nrow = 1, P_Ob.1 + theme(legend.position = "none"),
                   P_Sh.1 + theme(legend.position = "none"),
                   P_Si.1 + theme(legend.position = "none"), 
                   top = textGrob("Alpha diversity",gp=gpar(fontsize=25)))

Fig_3 <- grid.arrange(nrow = 2, P_alpha, get_plot_component(P_Sh.1 + theme(legend.position = "bottom"), pattern = "guide-box-bottom"), heights= c(6,1))

ggsave(plot = Fig_3, file.path(a.plots,"Figure_3.tiff"), dpi = 300, height = 8, width = 12)
ggsave(plot = Fig_3, file.path(a.plots,"Figure_3.svg"), dpi = 300, height = 8, width = 12)
ggsave(plot = Fig_3, file.path(a.plots,"Figure_3.png"), dpi = 600, height = 8, width = 12)

#### End ####

################################################################################
################################ Beta diversity ################################
################################################################################

#### Ordination ####
# Select phyloseq object with relative abundance to calculate beta diversity
physeq_beta=physeq_re
sampdf= data.frame(sample_data(physeq_beta))

# Ordinate with method "NMDS" and "PCoA" on distance "Bray-Curtis" 

ord.nmds.bray <- ordinate(physeq_beta, method="NMDS", distance="bray")
ord.pcoa.bray <- ordinate(physeq_beta, method="PCoA", distance="bray")

#### End ####

#### Plot ordination with spider and hull plots plot_ordination (NMDS or PCoA) ####

# Select text format
textt <- theme(plot.title = element_text(hjust = 0.5),
              text = element_text(size = 22, color = "black", face = "plain"),
              axis.text = element_text(color = "black",face = "bold"),
              axis.title = element_text(color = "black",face = "bold"))

#' Function to create the input for the plot. It returns:
#' - ord_s: data frame used as input for ggplot2
#' - ord_h: mecessary for creating the hull plot

ordi <- function(ords, col.var, var2, var3){
  ordiplot <- plot_ordination(physeq_beta, ords, color="col.var", justDF =T)
  colnames(ordiplot) <- c("x.axis", "y.axis", colnames(ordiplot)[-c(1:2)])
  ordiplot$shape <- paste(ordiplot[[var2]], ordiplot[[var3]])
  ord_s <- data.frame(cluster=factor(ordiplot[,col.var]), x=ordiplot$x.axis, y=ordiplot$y.axis)
  centroids <- aggregate(cbind(x,y)~cluster,data=ord_s,mean)
  ord_s <- merge(ord_s,centroids,by="cluster",suffixes=c("",".centroid"))
  ord_s <<- merge(ord_s, ordiplot, by.x = "x", by.y = "x.axis", sort = F)
  ordiplot$hull <- ordiplot[, col.var]
  ord_h <<- group_by(ordiplot, hull) %>% slice(chull(x.axis, y.axis))
}

#' Function to create single plots. As we want to plot three variables with only two aesthetics,
#' we create a combined variable for shape, where the shape form will represent one variable
#' and the filling another one.
#' Finally, when we plot them together we will create a new legend.

spihull <- function(data, method, col.var, shape.var, color.name, shape.name, text.size = t_size, title = NULL){
  if(method =="NMDS"){x.lab="NMDS 1" ; y.lab="NMDS 2"}else if(method=="PCoA"){x.lab="PCoA 1" ; y.lab="PCoA 2"}
  ggplot(data, aes(x = x, y = y, color = {{col.var}}, shape = {{shape.var}})) + geom_point(size=2, alpha=0.5) + 
    geom_point(size=7, alpha=0.5) +
    geom_segment(aes(yend=y.centroid,xend=x.centroid,group={{col.var}}), size=1, alpha=0.5) +
    geom_polygon(data=ord_h,aes(x=x.axis,y=y.axis,fill={{col.var}},group={{col.var}}),alpha=0.1) +
    scale_shape_manual(values = c(0, 1, 2, 15, 16, 17), name = shape.name) +
    scale_colour_brewer(type="qual", palette="Set1", name = color.name) +
    guides(fill = F) +
    xlab(x.lab) + ylab(y.lab) +
    # stat_ellipse(mapping = aes(group = {{col.var}}), level = 0.95) +
    ggtitle(title) + textt
    # theme(plot.title = element_text(hjust = 0.5),
    #       text = element_text(size = text.size, color = "black", face = "plain"),
    #       axis.text = element_text(color = "black",face = "bold"),
    #       axis.title = element_text(color = "black",face = "bold")) 
}

#' Function to create new legend

new_legend <- function(plot, shape.var, shape.name, fill.var, fill.name){
  
  shape_legend <- ggplot(ord_s, aes(shape={{shape.var}})) + geom_point(size=7, aes(x=x, y=y, shape={{shape.var}})) + 
    scale_shape_manual(values = shape, name = shape.name) + textt
  
  shape_fill_legend <- ggplot(ord_s, aes(shape={{fill.var}})) + geom_point(size=7, aes(x=x, y= y, shape={{fill.var}})) + 
    scale_shape_manual(values = shape_fill, name = fill.name) + textt
  
  plot_grid(plot + theme(legend.position="none"), plot_grid(
    get_plot_component(shape_legend, pattern = "guide-box-right"),
    get_plot_component(shape_fill_legend, pattern = "guide-box-right"),
    get_plot_component(plot + guides(shape = F), pattern = "guide-box-right"), nrow = 1),
    nrow = 2, rel_heights = c(8,2))
}

#### Tooth ####

# For PCoA
ordi(ord.pcoa.bray, "tooth", var2 = "implant", var3 = "time")
PCoA_tooth <- spihull(data = ord_s, method="PCoA", col.var = tooth, shape.var = shape, color.name = "Tooth", shape.name = "Implant - Time")

# For NMDS
ordi(ord.nmds.bray, "tooth", var2 = "implant", var3 = "time")
NMDS_tooth <- spihull(data = ord_s, method="NMDS", col.var = tooth, shape.var = shape, color.name = "Tooth", shape.name = "Implant - Time")

# Change legend
shape <- c(15, 16, 17) ; names(shape) <- levels(as.factor(sampdf$time))
shape_fill <- c(0,15) ; names(shape_fill) <- levels(as.factor(sampdf$implant))

PCoA_tooth.1 <- new_legend(PCoA_tooth, time, "Time", implant, "Implant")
PCoA_tooth.1

NMDS_tooth.1 <- new_legend(NMDS_tooth, time, "Time", implant, "Implant")
NMDS_tooth.1


#### Implant ####

ordi(ord.pcoa.bray, "implant", var2 = "tooth", var3 = "time")
PCoA_implant <- spihull(data = ord_s, method="PCoA", col.var = implant, shape.var = shape, color.name = "Implant", shape.name = "Tooth - Time")

ordi(ord.nmds.bray, "implant", var2 = "tooth", var3 = "time")
NMDS_implant <- spihull(data = ord_s, method="NMDS", col.var = implant, shape.var = shape, color.name = "Implant", shape.name = "Tooth - Time")

#Change legend
shape <- c(15, 16, 17) ; names(shape) <- levels(as.factor(sampdf$time))
shape_fill <- c(0,15) ; names(shape_fill) <- levels(as.factor(sampdf$tooth))

PCoA_implant.1 <- new_legend(PCoA_tooth, time, "Time", tooth, "Tooth")
PCoA_implant.1

NMDS_implant.1 <- new_legend(NMDS_tooth, time, "Time", tooth, "Tooth")
NMDS_implant.1

#### Time ####

ordi(ord.pcoa.bray, "time", var2 = "implant", var3 = "tooth")
PCoA_time <- spihull(data = ord_s, method="PCoA", col.var = time, shape.var = shape, color.name = "Time", shape.name = "Implant - Tooth")

ordi(ord.nmds.bray, "time", var2 = "implant", var3 = "tooth")
NMDS_time <- spihull(data = ord_s, method="NMDS", col.var = time, shape.var = shape, color.name = "Time", shape.name = "Implant - Tooth")

#Change legend
shape <- c(16, 17) ; names(shape) <- levels(as.factor(sampdf$tooth))
shape_fill <- c(1,16) ; names(shape_fill) <- levels(as.factor(sampdf$implant))

PCoA_time.1 <- new_legend(PCoA_time, tooth, "Tooth", implant, "Implant")
PCoA_time.1

NMDS_time.1 <- new_legend(NMDS_time, tooth, "Tooth", implant, "Implant")
NMDS_time.1

#### Combine variables ####

PCoA <- grid.arrange(plot_grid(PCoA_tooth.1, PCoA_implant.1, PCoA_time.1,
                               nrow = 1, labels = "AUTO"),
                     top = textGrob("PCoA ordination using Bray–Curtis dissimilarity index",
                                    gp=gpar(fontsize=23)))
PCoA

NMDS <- grid.arrange(plot_grid(NMDS_tooth.1, NMDS_implant.1, NMDS_time.1,
                               nrow = 1, labels = "AUTO"),
                     top = textGrob("NMDS ordination using Bray–Curtis dissimilarity index",
                                    gp=gpar(fontsize=23, fontface = "bold")))
NMDS

ggsave(plot = NMDS, file.path(a.plots,"Figure_4.tiff"), dpi = 300, height = 7, width = 16)
ggsave(plot = NMDS, file.path(a.plots,"Figure_4.svg"), dpi = 300, height = 7, width = 16)
ggsave(plot = NMDS, file.path(a.plots,"Figure_4.png"), dpi = 600, height = 7, width = 16)

grid.arrange(plot_grid(NMDS_tooth.1, NMDS_implant.1, NMDS_time.1,
                       PCoA_tooth.1, PCoA_implant.1, PCoA_time.1,
                       nrow = 2, labels = "AUTO"),
             top = textGrob("NMDS and PCoA ordination using Bray–Curtis dissimilarity index",
                            gp=gpar(fontsize=23)))

#### End ####

#### Non-parametric permutation based MANOVA (PERMANOVA) ####

dist = phyloseq::distance(physeq_beta, method="bray")

set.seed(1594)
test.adonis <- adonis2(dist ~ implant+tooth+time, data = sampdf, by="margin")
test.adonis

write.xlsx(test.adonis, file = file.path(b.stats,"test.adonis.xlsx"), rowNames = T)

#### End ####