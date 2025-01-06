library(AnnotationDbi)

#####
get_ensembl <- function(version = "v104") {
  gr_project_dir <- "/Users/chris/Desktop/20230130_gr_project"
  if (version == "v104") {
    message("Loading Ensembl104...")
    gencode <- loadDb(file = file.path(gr_project_dir, "input/ensembl/txdb.ensembl104.sqlite"))
  }
  return(gencode)
}

#####
keepStdChr <- function(gr) {
  message("With all chromosomes, including contigs : ", length(gr), " regions")
  # stdChr <- paste0("chr", c(seq(1:22), "X", "Y"))
  stdChr <- paste0(c(seq(1:22), "X", "Y", "MT"))
  gr_StdChr <- keepSeqlevels(gr, stdChr[stdChr %in% seqlevels(gr)], pruning.mode = "coarse")
  message("Keeping standard chromosomes : ", length(gr_StdChr), " regions")
  message("\t--> ", length(gr) - length(gr_StdChr), " regions removed")
  return(gr_StdChr)
}
