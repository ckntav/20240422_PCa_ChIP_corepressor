setwd("/Users/chris/Desktop/20240422_PCa_ChIP_corepressor")

library(tidyverse)
library(knitr)
library(kableExtra)
library(GenomicRanges)
source("scripts/ckn_utils/ckn_utils_overlaps.R")
source("scripts/ckn_utils/ckn_utils_load_vcap_peaks.R")

#
LNCaP_NCOR1 <- rtracklayer::import("output/chip-pipeline_PCA_corepressor-GRCh38/peak_call/LNCaP_DMSO_NCOR1/NCOR1/LNCaP_DMSO_NCOR1.NCOR1_peaks.narrowPeak.stdchr.bed")
LNCaP_SMRT <- rtracklayer::import("output/chip-pipeline_PCA_corepressor-GRCh38/peak_call/LNCaP_DMSO_SMRT/SMRT/LNCaP_DMSO_SMRT.SMRT_peaks.narrowPeak.stdchr.bed")

#
LNCaP_peaks <- GRangesList(
  "LNCaP_NCOR1" = LNCaP_NCOR1,
  "LNCaP_SMRT" = LNCaP_SMRT
)

sapply(LNCaP_peaks, length) %>% as.data.frame() %>% set_names("number_of_peaks") %>% kbl %>% kable_classic_2(full_width = F, html_font = "Arial", font_size = 20)
plotVenn3(LNCaP_peaks, labels = FALSE, legend = TRUE, quantities_val = list(type = c("counts")))
combMat_LNCaP <- generate_comb_mat(LNCaP_peaks)
displayUpSet(combMat_LNCaP)

# AR
AR_EtOH_peaks <- load_vcap_peaks(antibodies = "AR", treatments = "EtOH")
all_peaks <- c(LNCaP_peaks, GRangesList("VCaP_AR_EtOH" = AR_EtOH_peaks))
sapply(all_peaks, length) %>% as.data.frame() %>% set_names("number_of_peaks") %>% kbl %>% kable_classic_2(full_width = F, html_font = "Arial", font_size = 20)

plotVenn3(all_peaks, labels = FALSE, legend = TRUE, quantities_val = list(type = c("counts")))
combMat_all <- generate_comb_mat(all_peaks)
displayUpSet(combMat_all)


sapply(AR_peaks, length) %>% as.data.frame() %>% set_names("number_of_peaks") %>% kbl %>% kable_classic_2(full_width = F, html_font = "Arial", font_size = 20)
plotVenn3(AR_peaks, labels = FALSE, legend = TRUE, quantities_val = list(type = c("counts")))
combMat_AR <- generate_comb_mat(AR_peaks)
displayUpSet(combMat_AR, customSetOrder = paste(sep = "_", "VCap", "AR", condition_list, "rep1"))

# ER
ER_peaks <- load_vcap_peaks(antibodies = "ER")
sapply(ER_peaks, length) %>% as.data.frame() %>% set_names("number_of_peaks") %>% kbl %>% kable_classic_2(full_width = F, html_font = "Arial", font_size = 20)
plotVenn3(ER_peaks, labels = FALSE, legend = TRUE, quantities_val = list(type = c("counts")))
combMat_ER <- generate_comb_mat(ER_peaks)
displayUpSet(combMat_ER, customSetOrder = paste(sep = "_", "VCap", "ER", condition_list, "rep1"))

# MED1
MED1_peaks <- load_vcap_peaks(antibodies = "MED1")
sapply(MED1_peaks, length) %>% as.data.frame() %>% set_names("number_of_peaks") %>% kbl %>% kable_classic_2(full_width = F, html_font = "Arial", font_size = 20)
plotVenn3(MED1_peaks, labels = FALSE, legend = TRUE, quantities_val = list(type = c("counts")))
combMat_MED1 <- generate_comb_mat(MED1_peaks)
displayUpSet(combMat_MED1, customSetOrder = paste(sep = "_", "VCap", "MED1", condition_list, "rep1"))

##### per treatment
# EtOH
EtOH_peaks <- load_vcap_peaks(treatments = "EtOH")
sapply(EtOH_peaks, length) %>% as.data.frame() %>% set_names("number_of_peaks") %>% kbl %>% kable_classic_2(full_width = F, html_font = "Arial", font_size = 20)
plotVenn3(EtOH_peaks, labels = FALSE, legend = TRUE, quantities_val = list(type = c("counts")))
combMat_EtOH <- generate_comb_mat(EtOH_peaks)
displayUpSet(combMat_EtOH, customSetOrder = paste(sep = "_", "VCap", antibody_list, "EtOH", "rep1"))

# R1881
R1881_peaks <- load_vcap_peaks(treatments = "R1881")
sapply(R1881_peaks, length) %>% as.data.frame() %>% set_names("number_of_peaks") %>% kbl %>% kable_classic_2(full_width = F, html_font = "Arial", font_size = 20)
plotVenn3(R1881_peaks, labels = FALSE, legend = TRUE, quantities_val = list(type = c("counts")))
combMat_R1881 <- generate_comb_mat(R1881_peaks)
displayUpSet(combMat_R1881, customSetOrder = paste(sep = "_", "VCap", antibody_list, "R1881", "rep1"))

# E2
E2_peaks <- load_vcap_peaks(treatments = "E2")
sapply(E2_peaks, length) %>% as.data.frame() %>% set_names("number_of_peaks") %>% kbl %>% kable_classic_2(full_width = F, html_font = "Arial", font_size = 20)
plotVenn3(E2_peaks, labels = FALSE, legend = TRUE, quantities_val = list(type = c("counts")))
combMat_E2 <- generate_comb_mat(E2_peaks)
displayUpSet(combMat_E2, customSetOrder = paste(sep = "_", "VCap", antibody_list, "E2", "rep1"))

# R1881_E2
R1881_E2_peaks <- load_vcap_peaks(treatments = "R1881_E2")
sapply(R1881_E2_peaks, length) %>% as.data.frame() %>% set_names("number_of_peaks") %>% kbl %>% kable_classic_2(full_width = F, html_font = "Arial", font_size = 20)
plotVenn3(R1881_E2_peaks, labels = FALSE, legend = TRUE, quantities_val = list(type = c("counts")))
combMat_R1881_E2 <- generate_comb_mat(R1881_E2_peaks)
displayUpSet(combMat_R1881_E2, customSetOrder = paste(sep = "_", "VCap", antibody_list, "R1881_E2", "rep1"))
