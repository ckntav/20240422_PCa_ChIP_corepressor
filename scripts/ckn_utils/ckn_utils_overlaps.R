library(ComplexHeatmap)
library(eulerr)

##### Venn
plotVenn <- function(gr_list, labels = TRUE, legend = FALSE) {
  overlaps <- GenomicOperations::GenomicOverlaps(gr_list)
  mat <- overlaps@matrix > 0
  fit <- euler(mat, shape = "ellipse")
  p <- plot(fit, labels = labels, legend = legend,
            quantities = list(type = c("counts", "percent")),
            fills = list(fill = c("#2B70AB", "#CD3301", "#FFB027", "#449B2B", "#9370DB")),
            edges = list(alpha = 0))
  return(p)
}

##### Venn change colors
plotVenn2 <- function(gr_list, labels = TRUE) {
  overlaps <- GenomicOperations::GenomicOverlaps(gr_list)
  mat <- overlaps@matrix > 0
  fit <- euler(mat, shape = "ellipse")
  p <- plot(fit, labels = labels, legend = list(side = "bottom"),
            quantities = list(type = c("counts", "percent")),
            fills = list(fill = c("#2B70AB", "#FFB027", "#3EA742", "#CD3301", "#9370DB")),
            edges = list(alpha = 0))
  return(p)
}

#####
plotVenn3 <- function(gr_list, labels = TRUE, legend = FALSE,
                      colors = c("#2B70AB", "#FFB027", "#3EA742", "#CD3301", "#9370DB"),
                      quantities_val = list(type = c("counts", "percent"))) {
  overlaps <- GenomicOperations::GenomicOverlaps(gr_list)
  mat <- overlaps@matrix > 0
  fit <- euler(mat, shape = "ellipse")
  p <- plot(fit, labels = labels, legend = legend,
            quantities = quantities_val,
            fills = list(fill = colors),
            edges = list(alpha = 0))
  return(p)
}

##### UpSet plot
generate_comb_mat <- function(peak_sets) {
  inter_peaks <- GenomicOperations::GenomicOverlaps(peak_sets)
  matrix_peaks <- inter_peaks@matrix
  sum(matrix_peaks > 1)
  matrix_peaks[matrix_peaks > 1] <- 1
  sum(matrix_peaks > 1)
  colnames(matrix_peaks)
  
  comb_mat_peaks <- make_comb_mat(matrix_peaks)
  return(comb_mat_peaks)
}

get_all_regions <- function(peak_sets) {
  inter_peaks <- GenomicOperations::GenomicOverlaps(peak_sets)
  all_regions <- inter_peaks@regions
  return(all_regions)
}

#
displayUpSet <- function(combMat, threshold = 1, customSetOrder = NULL) {
  combMat <- combMat[comb_size(combMat) >= threshold]
  annot_top <- HeatmapAnnotation("Intersection\nsize" = anno_barplot(comb_size(combMat), 
                                                                     border = FALSE,
                                                                     gp = gpar(fill = "black"),
                                                                     height = unit(3, "cm")), 
                                 "Size" = anno_text(comb_size(combMat),
                                                    rot = 0,
                                                    just = "center",
                                                    location = 0.25),
                                 annotation_name_side = "left", annotation_name_rot = 0)
  annot_right <- rowAnnotation("Set size" = anno_barplot(set_size(combMat), 
                                                         border = FALSE, 
                                                         gp = gpar(fill = "black"), 
                                                         width = unit(2, "cm")),
                               "Size" = anno_text(set_size(combMat)))
  
  if (!is.null(customSetOrder)) {
    UpSet(combMat, top_annotation = annot_top, right_annotation = annot_right,
          set_order = customSetOrder)
  } else {
    UpSet(combMat, top_annotation = annot_top, right_annotation = annot_right)
  }
}