#' Extract marker genes as dataframe or gene list

#' @param object Seurat object
#' @param genes_only Logical if True only gene names are returned
#' @param n Number of genes
#' @param marker_genes Results of FindAllMarkers
#' @param per_cluster Logical, if TRUE marker_list contains df for each cluster seperately
#' @export 

marker_list <- function(object, marker_genes, genes_only = FALSE, n = "inf", per_cluster = TRUE){
    
    # Extract marker list
    marker_list <- list()
    for (i in levels(object)){
        marker_list[[i]] <- subset(marker_genes, cluster == i)
        if (n != "inf"){
            marker_list[[i]] <- marker_list[[i]][1:n,]
        }
        else {
            marker_list[[i]] <- marker_list[[i]]
        }
    }
    
    # Return only genes
    if (genes_only == TRUE){
        top_n_genes <- list()
        for (i in names(marker_list)){
            if (n != "inf"){
                top_n_genes[[i]] <- marker_list[[i]]$gene[1:n]
            }
            else {
                top_n_genes[[i]] <- marker_list[[i]]$gene
            }
        }
        if (per_cluster == TRUE){
            return(top_n_genes)
        }
        else {
            top_n_genes <- unique(unname(unlist(top_n_genes)))
        }
    }
    else {
        if (per_cluster == TRUE){
            return(marker_list)
        }
        else {
            marker_list <- do.call(rbind, marker_list)
            return(marker_list)
        }
    }
}