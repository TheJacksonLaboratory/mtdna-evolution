#--------------
# Packages
#--------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(require(circlize))
suppressPackageStartupMessages(library(maftools))

#--------------
# Load data
#--------------

analysis_dataset <- readRDS("~/Desktop/projects/mtDNA/rdata/analysis_dataset.RDS")
mito_genome <- readRDS("~/Desktop/projects/mtDNA/rdata/mito_genome.RDS")
all_variants_maf <- readRDS("~/Desktop/projects/mtDNA/rdata/all_variants_maf.RDS")
oncoplot_column_orders <- readRDS("~/Desktop/projects/mtDNA/rdata/oncoplot_column_orders.RDS")

#--------------
# Labels and colors
#--------------

#Primary and 1st Recurrecnce
timepoint_labs <- c("primary" = "Primary", "recurrent" = "Recurrent")
timepoint_colors <- c("primary" = "#A6611A", "recurrent" = "#018571")

#IDH subtypes- no codel
idh_subtype_labs <- c("IDHmut" = "IDH mutant", "IDHwt" = "IDH wild-type")
idh_colors <- c("IDHmut" = "#FAD000", "IDHwt" = "#158FAD")

#IDH subtypes- with codel
idh_codel_colors <- c("IDHmut-codel" = "#FAD000", "IDHmut-noncodel" = "#FAD000", "IDHwt" = "#158FAD")

#--------------
# A - Workflow
#--------------

# Made with Draw IO online

#--------------------------------------
# B - Num variants timepoint comparison
#--------------------------------------

analysis_dataset %>% 
  filter(timepoint %in% c(1,2)) %>%
  mutate(timepoint = ifelse(timepoint == 1, "primary", "recurrent"))%>%
  ggplot(.,aes(x=timepoint, y=mtDNA_SNV)) +
  geom_violin(aes(fill = timepoint)) +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  facet_wrap( ~ idh_subtype, labeller = labeller(idh_subtype = idh_subtype_labs)) +
  scale_x_discrete(labels = timepoint_labs)+
  scale_y_continuous(limits = c(0,300),breaks = sequence(20, from=0,by=2), expand = c(0,0)) +
  labs(title= NULL, x = NULL, y = "mtDNA variants", fill = "Timepoint") +
  scale_fill_manual(values=timepoint_colors) +
  stat_compare_means(method = "wilcox", paired = F, label.y = 19) +
  coord_cartesian(ylim = c(0,20)) +
  theme(legend.position  = "bottom") 

#----------------------------------
# C - Circular plot of mutations
#----------------------------------

#Create plot data

plotdata_circular <- all_variants_maf@data %>%
  dplyr::select(gene = Hugo_Symbol, type = Variant_Type,start = Start_Position) %>%
  group_by(gene, type, start) %>% summarise(total = n(), sector = "coding") 

#Create plot using circlize

# Outer track (bar plot with # of mutations per point) 
circos.par("points.overflow.warning" = FALSE,gap.after = c(8,0,0))
circos.initialize(c("coding","noncoding_start","noncoding_end"), 
                  xlim = data.frame(a = c(0,576,16024), 
                                    b =c(576,16024,16569), 
                                    row.names = c("noncoding_start","coding","noncoding_end")))
circos.track(ylim = c(0, 7), panel.fun = function(x, y) {
  circos.barplot(plotdata_circular$total, plotdata_circular$start,
                 sector.index = "coding", track.index = 1, bar_width = 0.8, 
                 col = 3)
  
  y = seq(1, 7, by = 1)
  circos.segments(577, y, 16024, y, sector.index = "coding",
                  track.index = 1, col = "#808080")
})
circos.yaxis(side = "right", sector.index = "coding",track.index = 1,
             labels.cex = 0.6, labels = c(1,2,3,4,5,6,7))
highlight.sector(c("noncoding_start","noncoding_end"), track.index = 1, col = "#808080" )


circos.clear()

# Inner plot with labels/colors for mitochondrial genome
par(new = TRUE) 
circos.par(cell.padding = c(0.002, 0, 0.002, 0),
           "points.overflow.warning" = FALSE,
           "canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-2, 2))
circos.initialize(mito_genome$gene, xlim = mito_genome[,2:3])
circos.track(mito_genome$gene,ylim = c(0, 1), track.height = 0.1)
#TYPES
highlight.sector(mito_genome$gene[which(mito_genome$type == "control")], track.index = 1, col = "#808080" )
highlight.sector(mito_genome$gene[which(mito_genome$type == "transferRNA")], track.index = 1, col = "#FFA500" )
highlight.sector(mito_genome$gene[which(mito_genome$type == "ribosomalRNA")], track.index = 1, col = "#4B0082")
highlight.sector(mito_genome$gene[which(mito_genome$type == "proteinCoding")], track.index = 1, col = "#008000")
#GENES
highlight.sector("D-LOOP", text = "D-LOOP",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-ATP6", text = "MT-ATP6",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-ATP8", text = "MT-ATP8",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-CO1", text = "MT-CO1",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-CO2", text = "MT-CO2",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-CO3", text = "MT-CO3",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-CYB", text = "MT-CYB",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-ND1", text = "MT-ND1",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-ND2", text = "MT-ND2",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-ND3", text = "MT-ND3",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-ND4", text = "MT-ND4",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-ND4L", text = "MT-ND4L",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-ND5", text = "MT-ND5",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
highlight.sector("MT-ND6", text = "MT-ND6",facing = "clockwise", niceFacing = TRUE, text.vjust = "8mm", cex = 0.7)
circos.clear()

#----------------------------------------------------
# D - Primary Oncoplot (excluding silent variants)
#----------------------------------------------------

passed_variants_primary_frozen <- subsetMaf(all_variants_maf,
                                            tsb = analysis_dataset$aliquot_barcode[which(analysis_dataset$timepoint == 1)],
                                            query = "Variant_Classification != 'Silent'")

columnOrder_primary_frozen <- oncoplot_column_orders[[1]]

##### D- PRIMARY ########
oncoplot(passed_variants_primary_frozen,
         sampleOrder = columnOrder_primary_frozen,
         genes= mt_genes,
         keepGeneOrder = TRUE,
         GeneOrderSort = FALSE,
         gene_mar = 8,
         fontSize = 0.6,
         clinicalFeatures = c('idh_codel_subtype'),
         annotationColor  = list("idh_codel_subtype"=idh_codel_colors),
         removeNonMutated = FALSE,
         legendFontSize = 1.5,
         legend_height = 6,
         fill = TRUE,
         cohortSize = 39,
         showTitle = FALSE)

#----------------------------------------------------
# E - Recurrent Oncoplot (excluding silent variants)
#----------------------------------------------------
passed_variants_recurrent_frozen <- subsetMaf(all_variants_maf,
                                            tsb = analysis_dataset$aliquot_barcode[which(analysis_dataset$timepoint != 1)],
                                            query = "Variant_Classification != 'Silent'")

columnOrder_recurrent_frozen <- oncoplot_column_orders[[2]]

oncoplot(passed_variants_recurrent_frozen,
         sampleOrder = columnOrder_recurrent_frozen,
         genes= mt_genes,
         keepGeneOrder = TRUE,
         GeneOrderSort = FALSE,
         gene_mar = 8,
         fontSize = 0.6,
         clinicalFeatures = c('idh_codel_subtype'),
         annotationColor  = list("idh_codel_subtype"=idh_codel_colors),
         logColBar = TRUE,
         removeNonMutated = FALSE,
         legendFontSize = 1.5,
         legend_height = 6,
         fill = TRUE,
         cohortSize = 44,
         showTitle = FALSE)
