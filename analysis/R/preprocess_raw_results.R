#--------------
# Packages
#--------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(maftools))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(odbc))

#-------------------------------------
# Connect to DB and source functions
#--------------------------------------

setwd("mtdna")

con <- DBI::dbConnect(odbc::odbc(), "VerhaakDB")

source("R/analysis/functions.R")

#--------------------------
# Get clinical data
#---------------------------

sql_clinical <- "SELECT a.aliquot_barcode, a.sample_barcode, s.*
FROM biospecimen.aliquots as a
INNER JOIN clinical.surgeries as s
ON a.sample_barcode = s.sample_barcode"

clinical_data <- dbGetQuery(con, sql_clinical)

clinical_data <- clinical_data %>% 
    select(Tumor_Sample_Barcode = aliquot_barcode, case_barcode:surgical_interval_mo, grade:codel_status, treatment_tmz, treatment_radiotherapy)

#--------------------------
# Read in raw data
#---------------------------

maf_file <- "data/consensus.vep.maf"
maf <- read.maf(maf_file)

#MAF FILE WITH ALL VARIANTS - read all variants regardless if silent or not
variant_classification <- append(as.character(unique(maf@data$Variant_Classification)),as.character(unique(maf@maf.silent$Variant_Classification)))

maf <- read.maf(maf_file,
                vc_nonSyn = variant_classification,
                clinicalData = clinical_data)

#mtDNA COVERAGE FILE
mtDNA_coverage<- read_csv("data/coverage.csv")

#nucDNA average coverage
gDNA_coverage <- readr::read_csv("data/gDNA_coverage.csv")
gDNA_coverage$aliquot_barcode <- as.character(lapply(gDNA_coverage$file,
                                                     function(x){unlist(strsplit(x,"\\."))[1]}))

#---------------------------------------
# Get titan and sequenza calls from DB
#----------------------------------------

titan <- dbGetQuery(con, "select * from variants.titan_params where pair_barcode like '%-WGS'")

sequenza <- dbGetQuery(con, "select * from variants.seqz_params where pair_barcode like '%-WGS'")


######  RESTART HERE- adding in conserved variants that got filtered out #####


#--------------------------
# Subset MAF File
#---------------------------

passed <- subsetMaf(maf,query = "FILTER == 'PASS'")

conserved_filtered_variants <- maf@data %>% 
  select(aliquot_barcode = Tumor_Sample_Barcode, ref = Reference_Allele, alt = Tumor_Seq_Allele2,
         start = Start_Position, end = End_Position, FILTER) %>% rowwise() %>%
  mutate(aliquot_barcode = as.character(aliquot_barcode),
         case_barcode = case_barcode(aliquot_barcode)) %>%
  group_by(case_barcode,start, ref,alt) %>% 
  mutate(shared_case_variant = ifelse(length(unique(aliquot_barcode)) > 1, TRUE, FALSE )) %>%
  filter(shared_case_variant == TRUE & "PASS" %in% FILTER) %>% 
  ungroup() %>%
  distinct() %>%
  arrange(case_barcode,start,aliquot_barcode) 

#--------------------------
# Intermediate tables
#---------------------------
passed_n_mutations <- passed@variants.per.sample

all_mutations_called <- maf@variants.per.sample

nonsilent_mutations <- passed@data[-which(passed@data$Variant_Classification == "Silent"),]

nonsilent_mutations <- nonsilent_mutations %>%
  select(Tumor_Sample_Barcode,Hugo_Symbol) %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(genes = paste0(unique(Hugo_Symbol), collapse = ";"))

#------------------------------------------------------------------------------
# Make total variant counts table (n mutations, coverage, copy number)
#-------------------------------------------------------------------------------
total_variant_counts <- data.frame(all_WGS$aliquot_barcode[which(all_WGS$aliquot_barcode %in% all_silver$aliquot_barcode)])
names(total_variant_counts) <- c("aliquot_barcode")

#Columns
total_variant_counts$case_barcode <- as.character(lapply(total_variant_counts$aliquot_barcode,case_barcode))

total_variant_counts$variants_passed <- 0
total_variant_counts$variants_called <- 0
total_variant_counts$mutated_genes <- NA
total_variant_counts$mtDNA_coverage <- NA

total_variant_counts$idh_codel_subtype <- NA

total_variant_counts$primary_recurrent[which(total_variant_counts$aliquot_barcode %in% silver$tumor_barcode_a)] <- "primary"
total_variant_counts$primary_recurrent[which(total_variant_counts$aliquot_barcode %in% silver$tumor_barcode_b)] <- "recurrent"

#Fill in all the blanks
for (i in 1:length(total_variant_counts$aliquot_barcode)){

  a <- which( passed_n_mutations$Tumor_Sample_Barcode == total_variant_counts$aliquot_barcode[i] )
  if (!(length(a) == 0)){ total_variant_counts$variants_passed[i] <- passed_n_mutations$Variants[a] }

  b <- which( all_mutations_called$Tumor_Sample_Barcode == total_variant_counts$aliquot_barcode[i] )
  if (!(length(b) == 0)){ total_variant_counts$variants_called[i] <- all_mutations_called$Variants[b] }

  c <- which( coverage$aliquot_barcode == total_variant_counts$aliquot_barcode[i] )
  if (!(length(c) == 0)){ total_variant_counts$coverage[i] <- coverage$coverage_mt[c] }

  d <- which( nonsilent_mutations$Tumor_Sample_Barcode == total_variant_counts$aliquot_barcode[i] )
  if (!(length(d) == 0)){ total_variant_counts$mutated_genes[i] <- nonsilent_mutations$genes[d] }

  e <- which( clinical_subset$Tumor_Sample_Barcode == total_variant_counts$aliquot_barcode[i] )
  if (!(length(e) == 0)){ total_variant_counts$idh_codel_subtype[i] <- clinical_subset$idh_codel_subtype[e] }
  
  f <- which( coverage$aliquot_barcode == total_variant_counts$aliquot_barcode[i] )
  if (!(length(f) == 0)){ total_variant_counts$mtDNA_coverage[i] <- coverage$coverage_mt[f] }
  

}
#------------------------------------------------------------------------------
# Make mutation counts table (SNV wide view for paired samples)
#-------------------------------------------------------------------------------

mutation_counts_wide <- total_variant_counts %>%
  dplyr::select(case_barcode, aliquot_barcode, timepoint = primary_recurrent, variants_called, variants_passed, coverage) %>%
  pivot_wider(names_from = timepoint, values_from = c(aliquot_barcode,variants_called, variants_passed, coverage)) %>%
  dplyr::select(case_barcode, aliquot_barcode_primary, aliquot_barcode_recurrent,
                variants_passed_primary, variants_passed_recurrent,
                variants_called_primary, variants_called_recurrent,
                coverage_primary, coverage_recurrent)

#------------------------------------------------------------------------------
# Calculate mtDNA copy number based on PCAWG formula (for sequenza and titan)
#-------------------------------------------------------------------------------

#Format Sequenza table
sequenza$sample_barcode <- as.character(lapply(sequenza$pair_barcode,sample_barcode))

#Format coverage/ copy number table
coverage <- coverage[which(coverage$aliquot_barcode %in% all_silver$aliquot_barcode),]
coverage$sample_barcode <- as.character(lapply(coverage$aliquot_barcode,sample_barcode))
coverage <- coverage[,c(1,3,2)]

coverage$coverage_gDNA <- coverage$seqz_purity <- coverage$seqz_ploidy <- coverage$seqz_mtDNA_CN <- coverage$titan_purity <- coverage$titan_ploidy <- coverage$titan_mtDNA_CN <- NA

for (i in 1:length(coverage$aliquot_barcode)){
  a <- which(gDNA_coverage$aliquot_barcode == coverage$aliquot_barcode[i])
  if (! (length(a) == 0)){ coverage$coverage_gDNA[i] <- gDNA_coverage$mean_coverage[a] }

  b <- which(titan$aliquot_barcode == coverage$aliquot_barcode[i])
  if (! (length(b) == 0)){ 
    coverage$titan_purity[i] <- titan$purity[b] 
    coverage$titan_ploidy[i] <- titan$ploidy[b] 
    }

  c <- which(sequenza$sample_barcode == coverage$sample_barcode[i])
  if (! (length(c) == 0)){ 
    coverage$seqz_purity[i] <- sequenza$cellularity[c]
    coverage$seqz_ploidy[i] <- sequenza$ploidy[c]
  }

  coverage$seqz_mtDNA_CN[i] <- calc_mtDNA_CN(coverage$coverage_gDNA[i],
                                             coverage$coverage_mt[i],
                                             coverage$seqz_ploidy[i],
                                             coverage$seqz_purity[i])
  coverage$titan_mtDNA_CN[i] <- calc_mtDNA_CN(coverage$coverage_gDNA[i],
                                             coverage$coverage_mt[i],
                                             coverage$titan_ploidy[i],
                                             coverage$titan_purity[i])
}

coverage$case_barcode <- as.character(lapply(coverage$aliquot_barcode,case_barcode))
coverage$primary_recurrent[which(coverage$aliquot_barcode %in% silver$tumor_barcode_a)] <- "primary"
coverage$primary_recurrent[which(coverage$aliquot_barcode %in% silver$tumor_barcode_b)] <- "recurrent"

#------------------------------------------------------------------------------
# Make paired copy number  table (CN wide view for paired samples)
#-------------------------------------------------------------------------------

coverage_wide <- coverage %>%
  dplyr::select(case_barcode, aliquot_barcode, timepoint = primary_recurrent, seqz_mtDNA_CN, titan_mtDNA_CN, coverage_mt) %>%
  pivot_wider(names_from = timepoint, values_from = c(aliquot_barcode,seqz_mtDNA_CN, titan_mtDNA_CN, coverage_mt)) %>%
  dplyr::select(case_barcode, aliquot_barcode_primary, aliquot_barcode_recurrent,
                coverage_mt_primary, coverage_mt_recurrent,
                titan_mtDNA_CN_primary, titan_mtDNA_CN_recurrent,
                seqz_mtDNA_CN_primary, seqz_mtDNA_CN_recurrent)

#--------------------------
# Save work
#---------------------------

write_tsv(mutation_counts_wide, "~/Desktop/projects/mtDNA/data/paired_view_mtDNA_SNV.tsv")
write_tsv(total_variant_counts, "~/Desktop/projects/mtDNA/data/mtDNA_SNV.tsv")

write_tsv(coverage_wide, "~/Desktop/projects/mtDNA/data/paired_view_mtDNA_CN.tsv")
write_tsv(coverage, "~/Desktop/projects/mtDNA/data/mtDNA_CN.tsv")

saveRDS(maf,"~/Desktop/projects/mtDNA/rdata/mtDNA_maf.RDS")
saveRDS(passed,"~/Desktop/projects/mtDNA/rdata/mtDNA_maf_with_only_passed_variants.RDS")
saveRDS(primary,"~/Desktop/projects/mtDNA/rdata/mtDNA_primary_maf.RDS")
saveRDS(recurrent,"~/Desktop/projects/mtDNA/rdata/mtDNA_recurrent_maf.RDS")