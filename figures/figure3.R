#--------------
# Packages
#--------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(maftools))

#--------------
# Load data
#--------------

analysis_dataset <- readRDS("~/Desktop/projects/mtDNA/rdata/analysis_dataset.RDS")
all_variants_maf <- readRDS("~/Desktop/projects/mtDNA/rdata/all_variants_maf.RDS")
rna_analysis_dataset <- readRDS("~/Desktop/projects/mtDNA/rdata/rna_analysis_dataset.RDS")

#-------------------
# Labels and colors
#--------------------

#Primary and 1st Recurrecnce
timepoint_labs <- c("primary" = "Primary", "recurrent" = "Recurrent")
timepoint_colors <- c("primary" = "#A6611A", "recurrent" = "#018571")

#IDH subtypes- no codel
idh_subtype_labs <- c("IDHmut" = "IDH mutant", "IDHwt" = "IDH wild-type")
idh_colors <- c("IDHmut" = "#FAD000", "IDHwt" = "#158FAD")

#TCGA subtype
tcga_colors <- c("Classical" = "blue", "Mesenchymal" = "red", "Proneural" = "green")
tcga_comp <- list( c("Classical", "Mesenchymal"), c("Mesenchymal", "Proneural"), c("Classical", "Proneural") )

#----------------------------------------------------------------------
# A- mtDNA copynumber subtype  difference at primary/recurrent timepoint
#----------------------------------------------------------------------

analysis_dataset %>% 
  filter(timepoint %in% c(1,2)) %>%
  mutate(timepoint = ifelse(timepoint == 1, "primary", "recurrent")) %>%
  ggplot(., aes(x=idh_subtype, y=mtDNA_CN)) +
  geom_violin(aes(fill = idh_subtype)) +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  facet_grid(~ timepoint, labeller = labeller(timepoint = timepoint_labs) ) +
  scale_x_discrete(labels = idh_subtype_labs) +
  scale_y_continuous(breaks = sequence(3000, from=0,by=300)) +
  labs(title= NULL,x = NULL, y = "mtDNA copy number",fill = "IDH subytpe") +
  scale_fill_manual(values= idh_colors, labels = idh_subtype_labs) +  
  stat_compare_means(method = "wilcox") +
  theme(legend.position = "bottom")

#--------------------------------------------------------------------
# B- mtDNA copynumber at primary/recurrent timepoint by subtype
#--------------------------------------------------------------------

analysis_dataset %>% 
  filter(timepoint %in% c(1,2)) %>%
  mutate(timepoint = ifelse(timepoint == 1, "primary", "recurrent")) %>%
  ggplot(., aes(x=timepoint, y=mtDNA_CN)) +
  geom_violin(aes(fill = timepoint)) +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  facet_wrap( ~ idh_subtype, labeller = labeller(idh_subtype = idh_subtype_labs)) +
  scale_x_discrete(labels = c("primary" = "Primary", "recurrent" = "Recurrent")) +
  labs(title= NULL,x = NULL, y = "mtDNA copy number",fill = "Timepoint") +
  scale_fill_manual(values = timepoint_colors, labels = timepoint_labs ) + 
  scale_y_continuous(breaks = sequence(3000, from=0,by=300)) +
  stat_compare_means(method = "wilcox") +
  theme(legend.position = "bottom")

#-----------------------------------------------------
# C- mtDNA SNVs vs CN
#-----------------------------------------------------

analysis_dataset %>% 
  filter(timepoint %in% c(1,2)) %>%
  mutate(timepoint = ifelse(timepoint == 1, "primary", "recurrent")) %>%
  ggplot(.,aes(y = mtDNA_CN, x= mtDNA_SNV )) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  facet_wrap( ~ idh_subtype + timepoint, 
              labeller = labeller(idh_subtype = idh_subtype_labs, timepoint = timepoint_labs)) +
  labs(title = NULL, x = "mtDNA variants",y = "mtDNA copy number") +
  scale_x_continuous(limits = c(0,300),breaks = sequence(18, from=0,by=2)) +
  scale_y_continuous(breaks = sequence(2700, from=0,by=300)) +
  coord_cartesian(xlim=c(0,18), ylim=c(0,2700)) +
  stat_cor(method = "spearman")  

#-----------------------------------------------------
# D- mtDNA CN and truncating mutations at timepoint
#-----------------------------------------------------

all_variants_maf@data %>% 
  rowwise() %>%
  mutate(truncating_mutation = ifelse(across(all_effects, ~ grepl('stop_gained|frameshift_variant', .)), TRUE, FALSE)) %>%
  select(aliquot_barcode = Tumor_Sample_Barcode, truncating_mutation) %>%
  group_by(aliquot_barcode) %>%
  summarise(truncating_mutation = ifelse(TRUE %in% truncating_mutation, TRUE, FALSE)) %>%
  plyr::join(., analysis_dataset, type = "full") %>%
  filter(timepoint %in% c(1,2)) %>%
  mutate(timepoint = ifelse(timepoint == 1, "primary", "recurrent"),
         truncating_mutation = ifelse(is.na(truncating_mutation), FALSE, truncating_mutation),
         truncating_mutation = factor(truncating_mutation, levels = c(TRUE, FALSE))) %>%
  ggplot(.,aes(x = idh_subtype, y = mtDNA_CN, fill = truncating_mutation)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap( ~ timepoint, labeller = labeller(timepoint = timepoint_labs)) +
  labs(title = NULL,x = NULL, y = "mtDNA copy number") +
  scale_x_discrete(labels = idh_subtype_labs) +
  scale_fill_manual(values = c("TRUE" = "#F7F7F7", "FALSE" = "#555555"), name = "Truncating mutation") +
  stat_compare_means(method = "wilcox", size = 3.3) +
  scale_y_continuous(breaks = seq(3000, from=0,by=300))+
  theme(legend.position = "bottom")

#----------------------------------------------------------------------
# E- mtDNA copynumber vs mtdna gene enrichment
#----------------------------------------------------------------------

rna_analysis_dataset %>%
  filter(timepoint %in% c(1,2)) %>%
  mutate(timepoint = ifelse(timepoint == 1, "primary","recurrent")) %>%
  ggplot(., aes(y = mtDNA_genes, x= mtDNA_CN )) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_continuous(limits = c(-2,6.5), breaks = seq(6, from = -2, by = 1) ) +
  scale_x_continuous(limits = c(0,2600), breaks = seq(2500, from = 0, by = 500)) +
  facet_wrap( ~ idh_subtype + timepoint ,
              labeller = labeller(idh_subtype = idh_subtype_labs, timepoint = timepoint_labs)) +
  labs(title =NULL, y = "mtDNA gene enrichment", x = "mtDNA copy number") +
  theme_classic() +
  stat_cor(method = "spearman")

#----------------------------------------------------------------------
# F- TCGA subtypes mtdna gene enrichment
#----------------------------------------------------------------------

rna_analysis_dataset %>%
  filter(idh_subtype == "IDHwt" & timepoint %in% c(1,2)) %>%
  mutate(timepoint = ifelse(timepoint == 1, "primary","recurrent")) %>%
  ggplot(., aes(x = tcga_subtype, y = mtDNA_genes)) +
  geom_violin(aes(fill = tcga_subtype)) +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  facet_grid(~ timepoint, labeller = labeller(timepoint = timepoint_labs) ) +
  scale_y_continuous(limits = c(-2,8), breaks = seq(8, from = -2, by = 1) ) +
  scale_fill_manual(values = tcga_colors) +
  labs(title= NULL,x = NULL, y = "mtDNA enrichment score",fill = "TCGA Subtype") +  
  stat_compare_means(comparisons = tcga_comp) +
  theme(legend.position = "bottom") 
