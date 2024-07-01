library(GenomicRanges)

#
load_vcap_peaks <- function(antibodies = c("AR", "ER", "MED1"), treatments = c("EtOH", "R1881", "E2", "R1881_E2")) {
  project_dir <- "/Users/chris/Desktop/20240206_VCap_project"
  peaks_dir <- "output/chip-pipeline_AR_ER_MED1_rep1-GRCh38_PE/peak_call_withWCE"
  
  all_peaks <- GRangesList()
  names_all_peaks <- c()
  
  for (antibody in antibodies) {
    for (treatment in treatments) {
      sample_name <- paste(sep = "_", "VCap", antibody, treatment, "rep1")
      message("# ", sample_name)
      
      peaks_path <- file.path(project_dir, peaks_dir, sample_name, antibody, paste0(sample_name, ".", antibody, "_peaks.narrowPeak.stdchr.bed"))
      peaks <- rtracklayer::import(peaks_path)
      # message(" > ", peaks_path)
      message("   > Number of regions : ", length(peaks))
      
      # gather everything
      all_peaks <- append(all_peaks, GRangesList(peaks))
      names_all_peaks <- c(names_all_peaks, sample_name)
    }
  }
  
  names(all_peaks) <- names_all_peaks
  message("#####################################")
  message("Available set of regions: ")
  print(names(all_peaks))
  return(all_peaks)
}