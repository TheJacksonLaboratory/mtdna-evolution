#--------------
# Packages
#--------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(NMF))

#--------------
# Load data
#--------------

# Analysis datasets
analysis_dataset <- readRDS("~/Desktop/projects/mtDNA/rdata/analysis_dataset.RDS")
analysis_dataset_coverage <- readRDS("~/Desktop/projects/mtDNA/rdata/analysis_dataset_coverage.RDS")
analysis_dataset_clinical <- readRDS("~/Desktop/projects/mtDNA/rdata/analysis_dataset_clinical.RDS")
analysis_dataset_mitocarta <- readRDS("~/Desktop/projects/mtDNA/rdata/analysis_dataset_mitocarta.RDS")
analysis_dataset_turnover <- readRDS("~/Desktop/projects/mtDNA/rdata/analysis_dataset_turnover.RDS")

rna_analysis_dataset <- readRDS("~/Desktop/projects/mtDNA/rdata/rna_analysis_dataset.RDS")

# MAFs
all_variants_maf <- readRDS("~/Desktop/projects/mtDNA/rdata/all_variants_maf.RDS")

# MitoCarta 3.0 mutations
acquired_repair_mutations <- readRDS("~/Desktop/projects/mtDNA/rdata/acquired_repair_mutations.RDS")

# Survival data
load("~/Desktop/projects/mtDNA/rdata/IDHmutant_survival_vs_mtDNA_CN.RData")
load("~/Desktop/projects/mtDNA/rdata/IDHwildtype_survival_vs_mtDNA_CN.RData")
load("~/Desktop/projects/mtDNA/rdata/IDHmutant_survival_vs_mtDNA_CN_recurrent.RData")
load("~/Desktop/projects/mtDNA/rdata/IDHwildtype_survival_vs_mtDNA_CN_recurrent.RData")

#-------------------
# Labels and colors
#--------------------

#Primary and 1st Recurrecnce
timepoint_labs <- c("primary" = "Primary", "recurrent" = "Recurrent")
timepoint_colors <- c("primary" = "#A6611A", "recurrent" = "#018571")

#IDH subtypes- no codel
idh_subtype_labs <- c("IDHmut" = "IDH mutant", "IDHwt" = "IDH wild-type")
idh_colors <- c("IDHmut" = "#FAD000", "IDHwt" = "#158FAD")

#OXPHOS complexes
oxphos_complex_labs <- c("Complex_I" = "Complex I", "Complex_II" = "Complex II", "Complex_III" = "Complex III",
                         "Complex_IV" = "Complex IV", "Complex_V" = "Complex V")

#---------------------------------------------------
# 1- Coverage histogram w/ vertical line at median
#---------------------------------------------------

# Data: analysis_dataset_coverage

analysis_dataset_coverage %>%
  ggplot(., mapping = aes(x = mtDNA_coverage)) +
  geom_histogram(bins=2^5, colour = 1, fill = "white") +
  geom_vline(xintercept = ceiling(median(analysis_dataset_coverage$mtDNA_coverage)),
             linetype=2, color = "red", size=1.5) +
  theme_classic() +
  scale_x_continuous(limits = c(0,50000),
                     breaks = seq(50000, from = 0, by=5000), 
                     expand = c(0,0), 
                     labels = comma_format()) +
  scale_y_continuous(limits = c(0,20),
                     breaks = seq(20, from = 0, by=2),
                     expand = c(0, 0)) +
  labs(y = "Count", x = "Mitochondrial genome coverage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#--------------------------------------------
# 2- VAFs per subtype/timepoint
#--------------------------------------------

# Data: all_variants_maf, analysis_dataset

# Make plot
all_variants_maf@data %>%
  dplyr::select(gene = Hugo_Symbol, type = Variant_Type,
                start = Start_Position, end = End_Position,
                aliquot_barcode = Tumor_Sample_Barcode, 
                t_alt_count ,t_ref_count) %>%
  plyr::join(., analysis_dataset) %>%
  filter(timepoint %in% c(1,2)) %>%
  mutate(timepoint = ifelse(timepoint == 1, "primary", "recurrent"),
         vaf = t_alt_count / (t_alt_count + t_ref_count)) %>%
  ggplot(., aes(x = vaf)) +
  geom_histogram(breaks = seq(1, from = 0, by=0.05), colour = 1, fill = "white") + 
  theme_classic() +
  scale_x_continuous(limits = c(0,1), breaks = seq(1, from = 0, by=0.1), expand = c(0,0)) +
  facet_wrap(~ idh_subtype + timepoint, 
             scales = "free_y",
             labeller = labeller(idh_subtype = idh_subtype_labs,
                                 timepoint = timepoint_labs)) +
  scale_y_continuous(breaks = pretty_breaks(n = 8)) +
  labs(y = "Count", x = "Variant Allele Frequency (VAF)")

#--------------------------------------------
# 3- MT-CO1 mutations
#--------------------------------------------

# Data: all_variants_maf

# Subset MAFs
all_variants_primary <- subsetMaf(all_variants_maf,
                                            tsb = analysis_dataset$aliquot_barcode[which(analysis_dataset$timepoint == 1)])
all_variants_recurrent_wo_sn16 <- subsetMaf(all_variants_maf,
                                              tsb = analysis_dataset$aliquot_barcode[which(analysis_dataset$timepoint == 2 &
                                                                                             analysis_dataset$aliquot_barcode != "GLSS-SN-0016-R1-01D-WGS-5ZOY6H")])                                         
all_variants_sn16 <- subsetMaf(all_variants_maf,
                                            tsb = "GLSS-SN-0016-R1-01D-WGS-5ZOY6H")

# Lolliplots
lollipopPlot2(all_variants_primary,
              all_variants_recurrent_wo_sn16,
              gene = "MT-CO1",
              m1_name = "Primary",
              m2_name = "Recurrent",
              legendTxtSize = 1)

lollipopPlot(all_variants_sn16, gene = "MT-CO1",
             legendTxtSize = 1, showMutationRate = FALSE)

#--------------------------------------------
# 4- Recurrent mtDNA signature analysis
#--------------------------------------------

# Data: all_variants_maf

# Subset MAF
all_variants_recurrent_wo_sn16 <- subsetMaf(all_variants_maf,
                                              tsb = analysis_dataset$aliquot_barcode[which(analysis_dataset$timepoint != 1 &
                                                                                             analysis_dataset$aliquot_barcode != "GLSS-SN-0016-R1-01D-WGS-5ZOY6H")])                                         

#Create TNM and estimate signatures
tnm_mtDNA_recurrent <- trinucleotideMatrix(maf = all_variants_recurrent,
                           ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
est_mtDNA_recurrent <- estimateSignatures(mat = tnm_mtDNA_recurrent, pConstant = 0.1)

# Using estimate (see Cophenetic plot), extract signatures & compare to SBS
sigs_mtDNA_recurrent <- extractSignatures(mat = tnm_mtDNA_recurrent, n = 3, pConstant = 0.1)
comp_mtDNA_recurrent <- compareSignatures(sigs_mtDNA_recurrent, "SBS")

#Plot signatures
maftools::plotSignatures(nmfRes = sigs_mtDNA_recurrent, title_size = 1.2, sig_db = "SBS")

#--------------------------------------------
# 5- mtDNA SNV vs mitocarta
#--------------------------------------------

# Data : analysis_dataset_coverage, analysis_dataset_mitocarta, analysis_dataset_turnover

a5 <- analysis_dataset_mitocarta %>% 
  plyr::join(.,analysis_dataset_coverage) %>%
  filter(timepoint == 1) %>%
  mutate(mtDNA_coverage_adj_mut_freq = (mtDNA_SNV / mtDNA_coverage )*10^6) %>%
  ggplot(.,aes(y = mtDNA_coverage_adj_mut_freq, x= mitocarta_coverage_adj_mut_freq)) +
  geom_point() +
  geom_smooth(method = "lm",fullrange=TRUE) +
  theme_classic() +
  facet_wrap( ~ idh_subtype, labeller = labeller(idh_subtype = idh_subtype_labs)) +
  labs(title =NULL, x = NULL, y = "mtDNA SNV burden",) +
  scale_x_continuous(limits = c(0,2400),breaks = sequence(2200, from=0,by=400)) +
  scale_y_continuous(limits = c(-600,1500),breaks = sequence(15, from=-600,by=150)) +
  coord_cartesian(ylim=c(0,1500), xlim=c(0,2400)) +
  stat_cor(method = "spearman")

b5 <- analysis_dataset_turnover %>%
  filter(nuc_coverage_adj_mut_freq < 8) %>%
  ggplot(.,aes(y = mtDNA_turnover, x= mitocarta_coverage_adj_mut_freq)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap( ~ idh_subtype, labeller = labeller(idh_subtype = idh_subtype_labs)) +
  labs(title =NULL, y = "mtDNA SNV turnover",
       x = "Mitochondrial-associated nucDNA SNV burden") +
  theme_classic() +
  scale_y_continuous(limits = c(-2,14),breaks = sequence(14, from=-2,by=2)) +
  scale_x_continuous(limits = c(0,400),breaks = sequence(400, from=0,by=50)) +
  coord_cartesian(ylim=c(0,14), xlim=c(0,350)) +
  stat_cor(method = "spearman")

plot_grid(a5,b5, nrow = 2,ncol = 1,align = "v", labels="AUTO")

#--------------------------------------------
# 6- Turnover and repair mutations
#--------------------------------------------

# A

acquired_repair_mutations <- acquired_repair_mutations %>% filter(gene_symbol != "PRIMPOL")

analysis_dataset_turnover %>% 
  mutate(repair_mut = ifelse(aliquot_barcode %in% unique(acquired_repair_mutations$aliquot_barcode),
                             TRUE, FALSE)) %>%
  ggplot(.,aes(y = turnover, x = repair_mut)) +
  geom_violin(aes(fill = repair_mut)) +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  scale_y_continuous(limits = c(0,22), breaks = seq(22,from = 0,by = 2)) +
  facet_wrap( ~ idh_subtype , labeller = labeller(idh_subtype = idh_subtype_labs)) +
  scale_fill_manual(values = c("#555555","#F7F7F7"), name = "Gained mtDNA repair mutation") +
  labs(x = NULL, y = "mtDNA SNV turnover") +
  geom_text(data=data.frame(x=1.5, y=21, idh_subtype=c("IDHmut", "IDHwt"),
                            label = c("Wilcoxon, p = 0.1","Wilcoxon, p = 0.022")), 
            aes(x,y,label=label), inherit.aes=FALSE) +
  theme(legend.position = "bottom")


# B
GLASS_III_mitocarta_signatures <- readRDS("~/Desktop/projects/mtDNA/GLASS_III_mitocarta_signatures.RDS")
maftools::plotSignatures(nmfRes = GLASS_III_mitocarta_signatures, title_size = 1.2, sig_db = "SBS")

#--------------------------------------------
# 7- mtDNA SNV and CN
#--------------------------------------------  

# Data: analysis_dataset, analysis_dataset_clinical

a7 <- analysis_dataset %>% 
  filter(timepoint == 2 & nuc_coverage_adj_mut_freq < 8) %>%
  ggplot(.,aes(x = mtDNA_SNV, y= mtDNA_CN)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap( ~ idh_subtype, labeller = labeller(idh_subtype = idh_subtype_labs)) +
  labs(title =NULL, x = "mtDNA variants",y = "mtDNA copy number") +
  theme_classic() +
  scale_y_continuous(breaks = sequence(2700, from=0,by=300)) +
  scale_x_continuous(limits = c(0,12), breaks = seq(12,from = 0,by = 2)) +
  coord_cartesian(xlim=c(0,12)) +
  stat_cor(method = "spearman")

b7 <- analysis_dataset_clinical %>% 
  filter(timepoint == 2 & idh_subtype == "IDHmut") %>%
  filter(!(is.na(grade))) %>%
  mutate(grade = ifelse(grade == "IV", "Grade IV", "Low-grade (II/III)")) %>%
  ggplot(.,aes(x = mtDNA_SNV, y= mtDNA_CN)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap( ~ grade, labeller = labeller(idh_subtype = idh_subtype_labs)) +
  labs(title =NULL, x = "mtDNA variants",y = "mtDNA copy number") +
  theme_classic() +
  scale_y_continuous(breaks = sequence(2700, from=0,by=300)) +
  scale_x_continuous(limits = c(0,18), breaks = seq(18,from = 0,by = 2)) +
  coord_cartesian(ylim=c(0,2700)) +
  stat_cor(method = "spearman")

plot_grid(a7,b7, nrow = 2,ncol = 1, labels =  "AUTO")

#--------------------------------------------
# 8- VAFs per subtype of truncating mutations
#--------------------------------------------

# Data: all_variants_maf, analysis_dataset

all_variants_maf@data %>% 
  rowwise() %>%
  mutate(truncating_mutation = ifelse(across(all_effects, ~ grepl('stop_gained|frameshift_variant', .)), TRUE, FALSE)) %>%
  select(aliquot_barcode = Tumor_Sample_Barcode, 
         ref = Reference_Allele, alt = Tumor_Seq_Allele2,
         start = Start_Position, end = End_Position, gene = Hugo_Symbol,
         t_alt_count, t_ref_count, truncating_mutation) %>%
  filter(truncating_mutation == TRUE) %>%
  plyr::join(., analysis_dataset) %>%
  filter(timepoint %in% c(1,2)) %>%
  mutate(timepoint = ifelse(timepoint == 1, "primary", "recurrent"),
         vaf = t_alt_count / (t_alt_count + t_ref_count),
         vaf_adj = vaf * purity) %>%
  group_by(case_barcode,ref,alt,start) %>% 
  mutate(timepoint = ifelse(n() > 1, "shared", timepoint)) %>%
  ungroup() %>%
  ggplot(., mapping = aes(x=start, y = vaf_adj, color = timepoint)) +
  geom_point(size = 1) +
  theme_classic() +
  facet_wrap(~ idh_subtype, labeller = labeller(idh_subtype = idh_subtype_labs)) +
  scale_x_continuous(limits = c(0,16569), breaks = seq(16000, from = 0, by=1000), expand = c(0, 0),
                     labels = c("Start","1kb","2kb","3kb","4kb","5kb","6kb","7kb","8kb","9kb",
                                "10kb","11kb","12kb","13kb","14kb","15kb","16kb")) +
  scale_y_continuous(limits = c(0,0.7),
                     breaks = seq(0.7, from = 0, by=0.1),
                     expand = c(0, 0)) +
  scale_color_manual(values = c("primary" = "#A6611A", "shared" ="darkmagenta","recurrent" = "#018571"),
                     labels = c("primary" = "Primary", "shared" ="Shared","recurrent" = "Recurrent")) +
  labs(y = "Purity Adj. VAF", x = "Mitochondrial genome coordinates", color = "Timepoint") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "bottom")

#--------------------------------------------
# 9 - OXPHOS enrichment
#--------------------------------------------

# Data: rna_analysis_dataset

# IDH mutant plot
oxphos_complex_comp_IDHmut <- rna_analysis_dataset %>% 
  select(case_barcode:idh_subtype, Complex_I:Complex_V) %>%
  filter(timepoint %in% c(1,2) & idh_subtype == "IDHmut") %>%
  mutate(timepoint = ifelse(timepoint == 1, "primary","recurrent")) %>%
  group_by(case_barcode) %>% 
  filter(n() > 1) %>% 
  ungroup() %>%
  pivot_longer(cols = Complex_I:Complex_V, names_to = "pathway", values_to = "enrichment_score") %>%
  ggplot(., aes(x = timepoint, y = enrichment_score)) +
  geom_violin(aes(fill = timepoint)) +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  facet_wrap(~ pathway, labeller = labeller(pathway = oxphos_complex_labs)) +
  scale_fill_manual(values = timepoint_colors, labels = timepoint_labs) +
  scale_y_continuous(limits = c(-4,4),breaks = seq(4, from = -4, by = 1)) +
  scale_x_discrete(labels = timepoint_labs) +
  labs(title= NULL,x = NULL, y = "OXPHOS enrichment score",fill = "Timepoint") +  
  stat_compare_means(method = "wilcox", paired = T,label.y = 3.5) +
  theme(legend.position = "none")

# IDH wild-type plot
oxphos_complex_comp_IDHwt <- rna_analysis_dataset %>% 
  select(case_barcode:idh_subtype, Complex_I:Complex_V) %>%
  filter(timepoint %in% c(1,2) & idh_subtype == "IDHwt") %>%
  mutate(timepoint = ifelse(timepoint == 1, "primary","recurrent")) %>%
  group_by(case_barcode) %>% 
  filter(n() > 1) %>% 
  ungroup() %>%
  pivot_longer(cols = Complex_I:Complex_V, names_to = "pathway", values_to = "enrichment_score") %>%
  ggplot(., aes(x = timepoint, y = enrichment_score)) +
  geom_violin(aes(fill = timepoint)) +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  facet_wrap(~ pathway, labeller = labeller(pathway = oxphos_complex_labs)) +
  scale_fill_manual(values = timepoint_colors, labels = timepoint_labs) +
  scale_y_continuous(limits = c(-4,4),breaks = seq(4, from = -4, by = 1)) +
  scale_x_discrete(labels = timepoint_labs) +
  labs(title= NULL,x = NULL, y = "OXPHOS enrichment score",fill = "Timepoint") +  
  stat_compare_means(method = "wilcox", paired = T,label.y = 3.5) +
  theme(legend.position = "bottom")

#Get shared legend
timepoint_legend <- get_legend(oxphos_complex_comp_IDHwt)
oxphos_complex_comp_IDHwt <- oxphos_complex_comp_IDHwt + theme(legend.position = "none")

#Plot both subtypes with shared legend
plot_grid(plot_grid(oxphos_complex_comp_IDHmut,
                    oxphos_complex_comp_IDHwt, 
                    nrow = 1,ncol = 2, labels =  "AUTO"),
          timepoint_legend, 
          nrow = 2,ncol = 1, rel_heights = c(1,.2))

#--------------------------------------------
# 10 - Survival CN primary
#--------------------------------------------

# Data: IDHmutant_survival_vs_mtDNA_CN.RData, IDHwildtype_survival_vs_mtDNA_CN.RData

#IDH mutant plots
survplot_IDHmut_CN_primary <- ggsurvplot(survfit_IDHmut_CN, 
                                       data = coxdf_IDHmut_CN, 
                                       conf.int = FALSE, 
                                       xlim = c(0,150),
                                       break.x.by = 25,
                                       legend.labs=c("Low mtDNA CN", "High mtDNA CN"),
                                       xlab = "Overall Survival (months)",
                                       ggtheme = theme_classic())

hazard_ratio_IDHmut_CN_primary <- ggforest(fit_cox_IDHmut_CN, 
                                         data = survdata_IDHmut_CN,
                                         noDigits = 3, 
                                         fontsize = 0.5,
                                         main = NULL)

#IDH wild-type plots
survplot_IDHwt_CN_primary <- ggsurvplot(survfit_IDHwt_CN, 
                                      data = coxdf_IDHwt_CN,
                                      conf.int = FALSE,
                                      xlim = c(0,60),
                                      break.x.by = 10,
                                      xlab = "Overall Survival (months)",
                                      ggtheme = theme_classic(),
                                      legend = "none")

hazard_ratio_IDHwt_CN_primary <- ggforest(fit_cox_IDHwt_CN, 
                                        data = survdata_IDHwt_CN,
                                        noDigits = 3, 
                                        fontsize = 0.5,
                                        main = NULL)

# Get shared legend
surv_CN_legend <- get_legend(survplot_IDHmut_CN_primary$plot)

survplot_IDHmut_CN_primary$plot$theme$legend.position <- "none"

# Plot legend and survival plots separately
plot_grid( surv_CN_legend,
           plot_grid(survplot_IDHmut_CN_primary$plot, survplot_IDHwt_CN_primary$plot,
                     hazard_ratio_IDHmut_CN_primary, hazard_ratio_IDHwt_CN_primary,
                     nrow = 2,ncol = 2, rel_heights = c(15,10)),
           nrow = 2,ncol = 1, rel_heights = c(1,10))

#--------------------------------------------
# 11 - Survival CN recurrent
#--------------------------------------------

# Data: IDHmutant_survival_vs_mtDNA_CN_recurrent.RData, IDHwildtype_survival_vs_mtDNA_CN_recurrent.RData

survplot_IDHmut_CN_recurrent <- ggsurvplot(survfit_IDHmut_CN_recurrent, 
                                       data = coxdf_IDHmut_CN_recurrent, 
                                       conf.int = FALSE, 
                                       xlim = c(0,150),
                                       break.x.by = 25,
                                       xlab = "Overall Survival (months)",
                                       ggtheme = theme_classic(),
                                       legend = "none")

hazard_ratio_IDHmut_CN_recurrent <- ggforest(fit_cox_IDHmut_CN_recurrent, 
                                         data = survdata_IDHmut_CN_recurrent,
                                         noDigits = 3, 
                                         fontsize = 0.5,
                                         main = NULL)

#IDH wild-type plots
survplot_IDHwt_CN_recurrent <- ggsurvplot(survfit_IDHwt_CN_recurrent, 
                                      data = coxdf_IDHwt_CN_recurrent,
                                      conf.int = FALSE,
                                      xlim = c(0,60),
                                      break.x.by = 10,
                                      xlab = "Overall Survival (months)",
                                      ggtheme = theme_classic(),
                                      legend = "none")

hazard_ratio_IDHwt_CN_recurrent <- ggforest(fit_cox_IDHwt_CN_recurrent, 
                                        data = survdata_IDHwt_CN_recurrent,
                                        noDigits = 3, 
                                        fontsize = 0.5,
                                        main = NULL)

# Plot legend and survival plots separately - reuse legend from Supp. 10
plot_grid( surv_CN_legend,
           plot_grid(survplot_IDHmut_CN_recurrent$plot, survplot_IDHwt_CN_recurrent$plot,
                     hazard_ratio_IDHmut_CN_recurrent, hazard_ratio_IDHwt_CN_recurrent,
                     nrow = 2,ncol = 2, rel_heights = c(15,10)),
           nrow = 2,ncol = 1, rel_heights = c(1,10))