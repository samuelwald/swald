#' Extract marker genes as dataframe or gene list

#' @param object Seurat object
#' @param genes_only Logical if True only gene names are returned
#' @param n Number of genes
#' @export 

marker_list <- function(object, marker_genes, genes_only = FALSE, n = "inf"){
    
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
        return(top_n_genes)
    }
    else {
        return(marker_list)
    }
}