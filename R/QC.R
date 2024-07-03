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

#' Basic function to convert human to mouse gene names

#' @export 
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}