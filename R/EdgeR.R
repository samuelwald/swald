#' Functions for analaysing single cell data with EdgeR
#' 
#' @param names Names of your clusters
#' 
#' @export 
build_contrasts_cluster <- function(names){
    rownames <- names
    colnames <- list()
    for (i in names){
            colnames[[i]] <- paste0(gsub("_tp.*", "", i), "_vs_rest")
            colnames <- unique(unname(unlist(colnames)))
    }
    
    df <- data.frame(matrix(ncol = length(colnames), nrow = length(rownames)))
    rownames(df) <- rownames
    colnames(df) <- colnames
    for (j in rownames(df)){
        for (i in colnames(df)){
            if (startsWith(j, gsub("_vs_rest.*", "", i))){
                df[j, i] <- 1/4
            }
            else{
                df[j, i] <- 1/-(length(rownames(df))-4)
            }
        }
    }
    return(df)
    
}

#' Build contrasts for timepoints

#' @param n_timepoints Number of timepoints

#' @export 
build_contrasts_timepoints <- function(names, n_timepoints){
    rownames <- names
    colnames <- list()
    for (i in names){
        colnames[[i]] <- i
        colnames <- unique(unname(unlist(colnames)))
    }
    
    df <- data.frame(matrix(ncol = length(colnames), nrow = length(rownames)))
    rownames(df) <- rownames
    colnames(df) <- colnames
    for (i in colnames(df)){
        for (j in rownames(df)){
            if (i == j){
                df[j,i] <- 1
            }
            else if (i != j & startsWith(i, gsub("_tp.*", "", j))){
                df[j,i] <- -1/(n_timepoints -1)
            }
            else {
                df[j,i] <- 0
            }

        }
    }
    return(df)
    
}


#' @export 

build_contrasts_cluster_timepoints <- function(names){
    rownames <- names
    colnames <- list()
    for (i in names){
        colnames[[i]] <- i
        colnames <- unique(unname(unlist(colnames)))
    }
    
    df <- data.frame(matrix(ncol = length(colnames), nrow = length(rownames)))
    rownames(df) <- rownames
    colnames(df) <- colnames
    for (i in colnames(df)){
        for (j in rownames(df)){
            if (i == j){
                df[j,i] <- 1
            }
            else {
                df[j,i] <- -1/(nrow(df)-1)
            }

        }
    }
    return(df)
    
}

#' Extract results from DEG
#' @param qlf Result from run_qlf
#' @param p_value P-value, default is 0.05
#' @param n Number of genes
#' @param log_FC Fold change, default is 0.58
#' @param direction Choose between up or downregulated genes
#' @param genes_only Returns only the gene names
#' @export 

DEG_edgeR <- function(
    qlf, 
    p_value = 0.05, 
    n = Inf, 
    log_FC = 0.58, 
    direction = "Up", 
    genes_only = FALSE
){
    
    # Up or Downregulated genes
    sig_genes <- list()
    
    # Upregulated
    if (isTRUE(direction == "Up")){
        for (i in names(qlf)){
            sig_genes[[i]] <- topTags(qlf[[i]],sort.by = "PValue", p.value = p_value, n = Inf)$table
            sig_genes[[i]] <-  sig_genes[[i]][sig_genes[[i]]$logFC > log_FC,]
        }
    }
    
    # Downregulated
    else if (isTRUE(direction == "Down")){
        for (i in names(qlf)){
            sig_genes[[i]] <- topTags(qlf[[i]],sort.by = "PValue", p.value = p_value, n = Inf)$table
            sig_genes[[i]] <-  sig_genes[[i]][sig_genes[[i]]$logFC < -log_FC,]
        }
    }
    
    # Extract only certain number of genes
    if (!isTRUE(n == Inf)){
        for (i in names(qlf)){
            sig_genes[[i]] <-  sig_genes[[i]][1:n,]
        }    
    }
                
                
    # Extract genes only
    if (isTRUE(genes_only == TRUE)){
        for (i in names(qlf)){
            sig_genes[[i]] <- rownames(sig_genes[[i]])
        }
        return(sig_genes)
    }
    else {
        return(sig_genes)
    }
    
    
}