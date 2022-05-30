#--------------
# Packages
#--------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(maftools))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggvenn))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(survival))

#--------------
# Load data
#--------------

analysis_dataset <- readRDS("~/Desktop/projects/mtDNA/rdata/analysis_dataset.RDS")
all_variants_maf <- readRDS("~/Desktop/projects/mtDNA/rdata/all_variants_maf.RDS")

load("~/Desktop/projects/mtDNA/rdata/IDHmutant_survival_vs_mtDNA_turnover.RData")
load("~/Desktop/projects/mtDNA/rdata/IDHwildtype_survival_vs_mtDNA_turnover.RData")

#-------------------
# Labels and colors
#--------------------

#Multiple timepoints colors
multi_time_labs <- c("1" = "Primary","2" = "1st recurrence","3" = "2nd recurrence")
multi_time_colors <- c("1" = "#A6611A", "2" = "#018571", "3" = "#011484")

#IDH subtypes- no codel
idh_subtype_labs <- c("IDHmut" = "IDH mutant", "IDHwt" = "IDH wild-type")
idh_colors <- c("IDHmut" = "#FAD000", "IDHwt" = "#158FAD")

#Multiple timepoints conserved variants 
multi_conserved_labs <- c("MT-CO1" = "MT-CO1", "MT-ND3" = "MT-ND3", "MT-ND1" = "MT-RNR2")

#--------------------------------------------
# A- Multiple timepoints n variants and copy number
#--------------------------------------------
analysis_dataset %>%
  group_by(case_barcode) %>% 
  filter(n() > 2) %>%
  dplyr::mutate(timepoint = as.character(timepoint)) %>%
  ggplot(., aes(x = timepoint, y = mtDNA_SNV, group = case_barcode)) +
  geom_line() +
  geom_point(aes(size = mtDNA_CN, color = timepoint)) +
  facet_wrap(~ idh_subtype, scales = "free_x", labeller = labeller(idh_subtype = idh_subtype_labs)) +
  scale_color_manual(name = "Timepoint", values = multi_time_colors, guide = "none") +
  scale_size(name = "mtDNA copy number", limits = c(50,1700)) +
  scale_x_discrete(name = NULL, labels = multi_time_labs) +
  scale_y_continuous(name ="mtDNA variants", limits = c(0,32),
                     breaks = sequence(16, from = 0, by = 2), expand = c(0,0)) +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  

#--------------------------------------------
# B- Multiple timepoints shared mutations
#--------------------------------------------

all_variants_maf@data %>%
  dplyr::select(gene = Hugo_Symbol, type = Variant_Type,
                start = Start_Position, end = End_Position,
                aliquot_barcode = Tumor_Sample_Barcode) %>%
  plyr::join(.,analysis_dataset) %>%
  filter(case_barcode %in% c(analysis_dataset %>% 
                               group_by(case_barcode) %>% 
                               filter(n() > 2) %>% 
                               select(case_barcode))$case_barcode) %>%
  rowwise() %>%
  mutate(variant_id = paste(gene,start,end, sep = "-",collapse  = "-"),
         value = TRUE) %>%
  select(variant_id, timepoint, value) %>%
  pivot_wider(names_from = timepoint, names_prefix = "timepoint", names_sort = TRUE,
              values_from  = value, values_fill = FALSE) %>%
  rename("Primary" = "timepoint1","1st Recurrence" = "timepoint2", "2nd Recurrence" = "timepoint3") %>%
  ggvenn(., 
         fill_color = c("#A6611A", "#018571", "#011484"),
         fill_alpha = 0.85, text_color = "white",
         stroke_size = 0.5, set_name_size = 4, text_size = 3)

#--------------------------------------------
# C- Multiple timepoints conserved ladder plot
#--------------------------------------------

all_variants_maf@data %>%
  dplyr::select(gene = Hugo_Symbol, type = Variant_Type,
                start = Start_Position, end = End_Position,
                aliquot_barcode = Tumor_Sample_Barcode, 
                t_alt_count ,t_ref_count) %>%
  plyr::join(.,analysis_dataset) %>%
  filter(case_barcode %in% c(analysis_dataset %>% 
                               group_by(case_barcode) %>% 
                               filter(n() > 2) %>% 
                               select(case_barcode))$case_barcode) %>%
  group_by(case_barcode, gene, start, end) %>%
  filter(n() > 2) %>%
  mutate(vaf = t_alt_count / (t_alt_count + t_ref_count),
         vaf_adj = vaf*purity) %>%
  mutate(timepoint = factor(timepoint,levels = c("1","2","3")),
         gene = factor(gene,levels = c("MT-CO1","MT-ND3","MT-ND1"))) %>%
  ggplot(., aes(x = timepoint, y = vaf_adj) ) +
  geom_line( aes( group = interaction(case_barcode,start) ), linetype = 1, orientation = "x") +
  geom_point( key_glyph = "point") +
  facet_wrap(~ idh_subtype + gene, 
             labeller = labeller(idh_subtype = idh_subtype_labs,
                                 gene = multi_conserved_labs)) +
  labs(y = "Purity adj. VAF", x = NULL, title = NULL) +
  scale_x_discrete(labels = multi_time_labs) +
  scale_y_continuous(limits = c(0,1), breaks = seq(1, from = 0, by = 0.1), expand = c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1))

#--------------------------------------------
# D- Survival and turnover
#--------------------------------------------

#IDH mutant plots
survplot_IDHmut_turnover <- ggsurvplot(survfit_IDHmut_turnover, 
                                       data = coxdf_IDHmut_turnover, 
                                       conf.int = FALSE, 
                                       legend.labs=c("Low mtDNA turnover", "High mtDNA turnover"),
                                       xlab = "Overall Survival (months)",
                                       xlim = c(0,150),
                                      break.x.by = 25,
                                       ggtheme = theme_classic())

hazard_ratio_IDHmut_turnover <- ggforest(fit_cox_IDHmut_turnover, 
                                         data = survdata_IDHmut_turnover,
                                         noDigits = 3, 
                                         fontsize = 0.5,
                                         main = NULL)

#IDH wild-type plots
survplot_IDHwt_turnover <- ggsurvplot(survfit_IDHwt_turnover, 
                                      data = coxdf_IDHwt_turnover,
                                      conf.int = FALSE,
                                      xlab = "Overall Survival (months)",
                                      xlim = c(0,60),
                                      break.x.by = 10,
                                      ggtheme = theme_classic(),
                                      legend = "none")

hazard_ratio_IDHwt_turnover <- ggforest(fit_cox_IDHwt_turnover, 
                                        data = survdata_IDHwt_turnover,
                                        noDigits = 3, 
                                        fontsize = 0.5,
                                        main = NULL)

# Get shared legend
surv_turnover_legend <- get_legend(survplot_IDHmut_turnover$plot)

survplot_IDHmut_turnover$plot$theme$legend.position <- "none"

# Plot legend and survival plots separately
plot_grid( surv_turnover_legend,
           plot_grid(survplot_IDHmut_turnover$plot,survplot_IDHwt_turnover$plot,
                     hazard_ratio_IDHmut_turnover,hazard_ratio_IDHwt_turnover,
                     nrow = 2,ncol = 2, rel_heights = c(15,10)),
           nrow = 2,ncol = 1, rel_heights = c(1,10))