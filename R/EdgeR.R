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
    genes_only = FALSE,
    all_genes = FALSE
){
    genes <- list()
    if (isTRUE(all_genes == T)){
        for (i in names(qlf)){
            genes[[i]] <- topTags(qlf[[i]], n = Inf)$table
        }
        return(genes)
    }
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


#' Run qlf for multiple tests

#' @export 



run_qlf <- function(fit, contrast){
    qlf <- list()
    for (i in colnames(contrast)){
        qlf[[i]] <- glmQLFTest(fit, contrast = contrast[,i])
    }
    return(qlf)
}


#' Statistic for up and downregulated genes
#' @param qlf Result from run_qlf
#' @param p_value P-value, default is 0.05
#' @param log_FC Fold change, default is 0.58
 
#' @export 

stats_edgeR <- function(
    qlf,
    p_value = 0.05,
    log_FC = 0.58
    ){
    
    # Build dataframe
    colnames <- names(qlf)
    rownames <- c("Up", "Down")
    df <- data.frame(matrix(ncol = length(colnames), nrow = length(rownames)))
    rownames(df) <- rownames
    colnames(df) <- colnames
    
    # Upregulated
    sig_genes_up <- list()
    for (i in names(qlf)){
        sig_genes_up[[i]] <- topTags(qlf[[i]],sort.by = "PValue", p.value = p_value, n = Inf)$table
        sig_genes_up[[i]] <-  sig_genes_up[[i]][sig_genes_up[[i]]$logFC > log_FC,]
        df[1,i] <- nrow(sig_genes_up[[i]])
    }

    # Downregulated
    sig_genes_down <- list()
    for (i in names(qlf)){
        sig_genes_down[[i]] <- topTags(qlf[[i]],sort.by = "PValue", p.value = p_value, n = Inf)$table
        sig_genes_down[[i]] <-  sig_genes_down[[i]][sig_genes_down[[i]]$logFC < -log_FC,]
        df[2,i] <- nrow(sig_genes_down[[i]])
    }
    return(df)
    
}


#' Run FGSEA with a list of DEG results
#' 
#' 
#' @export 

run_fgsea <- function(
    species = "mouse", 
    category, 
    subcategory, 
    DEG_results, 
    min_size = 15, 
    max_size = 500
){
    # Load and prepare geneset
    msigdbr_df <- msigdbr(species = species, category = category, subcategory = subcategory)
    msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
    
    # Build ranked gene list
    ranked_genes <- list()
    for (j in names(DEG_results)){
        for (i in 1:nrow((DEG_results[[j]]))){
            DEG_results[[j]]$rank[i] <- -log10(DEG_results[[j]]$PValue[i]) * sign(DEG_results[[j]]$logFC[i])
        }
        ranked_genes[[j]] <- DEG_results[[j]]$rank
        names(ranked_genes[[j]]) <- rownames(DEG_results[[j]])
        ranked_genes[[j]] <- sort(ranked_genes[[j]], decreasing = T)
    }
    
    # Run FGSEA
    res <- list()
    for (i in names(ranked_genes)){
        res[[i]] <- suppressWarnings(fgsea(pathways = msigdbr_list, stats = ranked_genes[[i]], minSize = 10))
    }
    return(res)
    
}


#' Rank genes by fc and p-value

#' @export 

ranked_genes <- function(DEG_results){
    ranked_genes <- list()
    for (j in names(DEG_results)){
        for (i in 1:nrow((DEG_results[[j]]))){
            DEG_results[[j]]$rank[i] <- -log10(DEG_results[[j]]$PValue[i]) * sign(DEG_results[[j]]$logFC[i])
        }
        ranked_genes[[j]] <- DEG_results[[j]]$rank
        names(ranked_genes[[j]]) <- rownames(DEG_results[[j]])
        ranked_genes[[j]] <- sort(ranked_genes[[j]], decreasing = T)
    }
    return(ranked_genes)
}

#' Rank genes by fc and p-value from seurat FindAllMarkers

#' @export 

ranked_genes_seurat <- function(DEG_results){
    ranked_genes <- list()
    for (j in names(DEG_results)){
        for (i in 1:nrow((DEG_results[[j]]))){
            DEG_results[[j]]$rank[i] <- -log10(DEG_results[[j]]$p_val_adj[i]) * sign(DEG_results[[j]]$avg_log2FC[i])
        }
        ranked_genes[[j]] <- DEG_results[[j]]$rank
        names(ranked_genes[[j]]) <- rownames(DEG_results[[j]])
        ranked_genes[[j]] <- sort(ranked_genes[[j]], decreasing = T)
    }
    return(ranked_genes)
}


#' Plot GSEA table
#' @param fgsea_res List with FGSEA results

#' @export 

plot_gsea_table <- function(
    fgsea_res, 
    msigdbr_list, 
    ranked_genes
){
    plot <- list()
    topPathwaysUp <- list()
    topPathwaysDown <- list()
    topPathways <- list()
    for (i in names(fgsea_res)){
        topPathwaysUp[[i]] <- fgsea_res[[i]][ES > 0][head(order(pval), n=10), pathway]
        topPathwaysDown[[i]] <- fgsea_res[[i]][ES < 0][head(order(pval), n=10), pathway]
        topPathways[[i]] <- c(topPathwaysUp[[i]], rev(topPathwaysDown[[i]]))
        plot[[i]] <- plotGseaTable(msigdbr_list[topPathways[[i]]], ranked_genes[[i]], fgsea_res[[i]], gseaParam=0.5)
    }
    return(plot)
}