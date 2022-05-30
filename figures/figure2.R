#--------------
# Packages
#--------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(maftools))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scales))

#--------------
# Load data
#--------------

analysis_dataset <- readRDS("~/Desktop/projects/mtDNA/rdata/analysis_dataset.RDS")
all_variants_maf <- readRDS("~/Desktop/projects/mtDNA/rdata/all_variants_maf.RDS")
analysis_dataset_turnover <- readRDS("~/Desktop/projects/mtDNA/rdata/analysis_dataset_turnover.RDS")

#--------------
# Labels and colors
#--------------

#Primary and 1st Recurrecnce
timepoint_labs <- c("primary" = "Primary", "recurrent" = "Recurrent")
timepoint_colors <- c("primary" = "#A6611A", "recurrent" = "#018571")

#IDH subtypes- no codel
idh_subtype_labs <- c("IDHmut" = "IDH mutant", "IDHwt" = "IDH wild-type")
idh_colors <- c("IDHmut" = "#FAD000", "IDHwt" = "#158FAD")

#Mutation substitutions
mut_patterns_colors <- brewer.pal(n = 12, name = "Set3")
names(mut_patterns_colors) <- c("C>A","C>G","C>T","G>A","G>C","G>T","A>C","A>G","A>T","T>A","T>C","T>G")

mut_patterns_labs <- c("A>C" = "T>G (H)", "A>G" = "T>C (H)", "A>T" = "T>A (H)",
                               "C>A" = "C>A (L)", "C>G" = "C>G (L)", "C>T" = "C>T (L)",
                               "G>A" = "C>T (H)", "G>C" = "C>G (H)", "G>T" = "C>A (H)",
                               "T>A" = "T>A (L)", "T>C" = "T>C (L)","T>G" = "T>G (L)")

#--------------------------------------------
# A- Initial, recurrent, and shared mutations
#--------------------------------------------

all_variants_maf@data %>%
  dplyr::select(gene = Hugo_Symbol, type = Variant_Type,
                start = Start_Position, end = End_Position,
                aliquot_barcode = Tumor_Sample_Barcode) %>%
  plyr::join(., analysis_dataset[,c("case_barcode","aliquot_barcode","timepoint", "idh_subtype")]) %>%
  filter(timepoint %in% c(1,2)) %>%
  mutate(timepoint = ifelse( timepoint == 1, "primary","recurrent")) %>%
  group_by(case_barcode,gene,start,end) %>%
  mutate(mutation_timepoint = ifelse(length(unique(timepoint)) > 1, "shared", 
                                      ifelse( "primary" %in% timepoint,"primary", "recurrent"))) %>% 
  ungroup() %>%
  group_by(case_barcode,idh_subtype,mutation_timepoint) %>%
  summarise(total = n()) %>%
  mutate(mutation_timepoint = factor(mutation_timepoint,levels = rev(c("primary","shared","recurrent" )))) %>%
  ggplot(., aes(x = case_barcode, y = total, fill = mutation_timepoint)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic( ) +
  theme(axis.ticks.x = element_blank(),axis.ticks.length = unit(0, "mm")) +
  scale_x_discrete(breaks = order(analysis_dataset$case_barcode),labels = NULL, expand = c(0,0)) +
  labs(title = NULL, x = NULL, y = "Proportion of mutations per patient", fill = "Mutation timepoint")+
  scale_fill_manual(values = list("primary" = "#A6611A", "shared" ="darkmagenta", "recurrent" = "#018571")) +
  facet_wrap(~ idh_subtype, scales = "free_x", labeller = labeller(idh_subtype = idh_subtype_labs)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "bottom")

#----------------
# B- Shared VAFs
#-----------------

all_variants_maf@data %>% 
  select(aliquot_barcode = Tumor_Sample_Barcode, 
         ref = Reference_Allele, alt = Tumor_Seq_Allele2,
         start = Start_Position, end = End_Position, gene = Hugo_Symbol,
         t_alt_count, t_ref_count) %>% 
  plyr::join(., analysis_dataset[,c("case_barcode","aliquot_barcode","timepoint", "idh_subtype")]) %>%
  filter(timepoint %in% c(1,2)) %>%
  mutate(timepoint = ifelse( timepoint == 1, "primary","recurrent"),
          vaf = t_alt_count / (t_alt_count + t_ref_count)) %>%
  group_by(case_barcode,start,ref,alt) %>% 
  filter(n() > 1) %>%
  arrange(case_barcode,start,timepoint) %>% 
  mutate(gene = factor(gene, levels = c("MT-ND1", "MT-CO1","MT-CYB","MT-ATP6", "MT-ND3","MT-ND5","MT-ND4","MT-ND4L"))) %>% 
  ggplot(., aes(x=timepoint, y = vaf)) +
  geom_line(aes(group = interaction(case_barcode,start), 
                color = idh_subtype), 
                linetype = 1, orientation = "x") +
  geom_point(key_glyph = "point") +
  facet_wrap(~ gene, nrow = 2) +
  labs(title = NULL, x = "Timepoint", y = "Variant Allele Frequency") +
  scale_x_discrete(labels = timepoint_labs)+
  scale_y_continuous(limits = c(0,1), breaks = seq(1,from = 0,by =.1)) +
  scale_color_manual(values = idh_colors, labels = idh_subtype_labs, name = "IDH Subtype") +
  theme_classic() +
  theme(legend.position = "bottom")

#----------------
# C- Nuclear vs MT2
#-----------------

ggplot(analysis_dataset_turnover,aes(y = mtDNA_turnover, x= mitocarta_coverage_adj_mut_freq)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap( ~ idh_subtype, labeller = labeller(idh_subtype = idh_subtype_labs)) +
  labs(title =NULL, y = "mtDNA SNV turnover (log2)",
       x = "Mitochondrial-associated nucDNA SNV burden (log10)") +
  theme_classic() +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 2),
                     breaks = trans_breaks(function(x) asinh(x / (2 * 1)) / log(2), 
                                           function(x) 2 * 1 * sinh(x * log(2)), n = 10)(c(0,260)),
                     labels = trans_format(function(x) asinh(x / (2 * 1)) / log(2), 
                                           format = label_number(accuracy = 1))) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = trans_breaks(function(x) asinh(x / (2 * 1)) / log(10), 
                                           function(x) 2 * 1 * sinh(x * log(10)), n = 6)(c(0,14600)),
                     labels = trans_format(function(x) asinh(x / (2 * 1)) / log(10), 
                                           format = label_number(accuracy = .5))) +
  stat_cor(method = "spearman")


#----------------------------
# D- Mutational patterns
#----------------------------

all_variants_maf@data %>%
  dplyr::select(gene = Hugo_Symbol, type = Variant_Type,
                start = Start_Position, end = End_Position,
                ref = Reference_Allele, alt = Tumor_Seq_Allele2,
                aliquot_barcode = Tumor_Sample_Barcode) %>%
  plyr::join(., analysis_dataset[,c("case_barcode","aliquot_barcode","timepoint", "idh_subtype")]) %>%
  filter(timepoint %in% c(1,2) & type == "SNP") %>%
  mutate(timepoint = ifelse( timepoint == 1, "primary","recurrent")) %>%
  rowwise() %>%
  mutate(substitution = paste(ref, ">",alt, sep = "")) %>%
  group_by(aliquot_barcode) %>%
  mutate(total_sample_mutations = n()) %>% 
  group_by(aliquot_barcode, timepoint, idh_subtype, substitution, total_sample_mutations) %>%
  summarise(count = n(),
            weight = count / total_sample_mutations) %>%
  distinct() %>%
  pivot_wider(names_from = substitution, values_from = c("count","weight")) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(c(everything(), -c(aliquot_barcode,timepoint,total_sample_mutations,idh_subtype)), 
               names_to = c("count_weight","substitution"), names_sep = "_") %>%
  pivot_wider(names_from = count_weight, values_from = value) %>%
  dplyr::select(aliquot_barcode,idh_subtype,timepoint,total_sample_mutations,substitution,count,weight) %>% 
  filter(count != 0) %>%
  group_by(substitution,idh_subtype,timepoint) %>%
  mutate(mean_weight_subtype_timepoint = mean(weight),
         unweighted_count_subtype_timepoint = sum(count)) %>% 
  ungroup() %>%
  mutate(weighted_sample_count = mean_weight_subtype_timepoint / weight) %>% 
  group_by(substitution,idh_subtype,timepoint) %>%
  summarise(total_weighted_count = sum(weighted_sample_count)) %>%
  mutate(substitution = factor(substitution,
                               levels = c("A>C","A>G","A>T",
                                          "C>A","C>G","C>T",
                                          "G>A","G>C","G>T",
                                          "T>A","T>C","T>G"))) %>%
  ggplot(., aes(x = timepoint, y = total_weighted_count, fill = substitution)) +
  geom_bar(stat = "identity", position = position_fill()) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(title = NULL, x = "Timepoint", y = "Proportion of substitutions", fill = NULL)+
  scale_fill_manual(values = mut_patterns_colors, labels = mut_patterns_labs) +
  facet_wrap(~ idh_subtype, labeller = labeller(idh_subtype = idh_subtype_labs)) +
  guides(fill = guide_legend(nrow =2)) +
  scale_y_continuous(expand = c(0,0) ,breaks = seq(1,from = 0,by =.1)) +
  scale_x_discrete(expand = c(0,0), labels = timepoint_labs)

#---------------------------------
# E- mtDNA hypermutant case VAFs
#---------------------------------

all_variants_maf@data %>% 
  select(aliquot_barcode = Tumor_Sample_Barcode, 
         ref = Reference_Allele, alt = Tumor_Seq_Allele2,
         start = Start_Position, end = End_Position, gene = Hugo_Symbol,
         t_alt_count, t_ref_count) %>% 
  filter(aliquot_barcode == "GLSS-SN-0016-R1-01D-WGS-5ZOY6H") %>%
  mutate(vaf = t_alt_count / (t_alt_count + t_ref_count)) %>%
  ggplot(., mapping = aes(x=start, y = vaf, color = ref)) +
  geom_point(size = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0,16569),breaks = seq(16000, from = 0, by=1000), expand = c(0, 0),
                       labels = c("Start","1kb","2kb","3kb","4kb","5kb","6kb","7kb","8kb","9kb",
                                  "10kb","11kb","12kb","13kb","14kb","15kb","16kb")) +
  scale_y_continuous(limits = c(0,0.1),breaks = seq(0.1, from = 0, by=0.01),expand = c(0, 0)) +
  scale_color_manual(values = c("G" = "#FDB462", "C" = "#8DD3C7"), 
                     labels = c("G" = "C>A (H)","C" = "C>A (L)"), 
                     name = "Substitution") +
  labs(y = "Variant Allele Frequency", x = "Mitochondrial genome coordinates") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


