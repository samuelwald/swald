#' Functions to define outliers for QC
#' 
#' @param seurat_object Your seurat object

#' @param metric Mtric where you want to define ouliers eg mitochondrial genes

#' @param nmads Median Absolute Deviation

#' @export 


is_outlier <- function(seurat_object, metric, nmads = 5){
    M <- seurat_object[[]][,metric]
    outlier <- (M < median(M) - nmads * mad(M)) | (
        median(M) + nmads * mad(M) < M
    )
    return(outlier)
}