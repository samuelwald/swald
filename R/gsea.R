#' Run FGSEA with a list of DEG results from EdgeR

#' @export 

run_fgsea_edgeR <- function(
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
        res[[i]] <- suppressWarnings(fgsea(pathways = msigdbr_list, stats = ranked_genes[[i]], minSize = min_size))
    }
    return(res)
    
}


#' Run FGSEA with a list of DEG results from Seurat

#' @export 


run_fgsea_seurat <- function(
    species = "mouse", 
    category, 
    subcategory = NULL, 
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
            DEG_results[[j]]$rank[i] <- -log10(DEG_results[[j]]$p_val_adj[i]) * sign(DEG_results[[j]]$avg_log2FC[i])
        }
        ranked_genes[[j]] <- DEG_results[[j]]$rank
        names(ranked_genes[[j]]) <- DEG_results[[j]]$gene
        ranked_genes[[j]] <- sort(ranked_genes[[j]], decreasing = T)
    }
    
    # Run FGSEA
    res <- list()
    for (i in names(ranked_genes)){
        res[[i]] <- suppressWarnings(fgsea(pathways = msigdbr_list, stats = ranked_genes[[i]], minSize = min_size))
    }
    return(res)
    
}


#' Rank genes (from edgeR) by fc and p-value

#' @export 

ranked_genes_edgeR <- function(DEG_results){
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



#' Rank genes (from Seurat) by fc and p-value

#' @export 
#' 
#' 
ranked_genes_seurat <- function(DEG_results){
    ranked_genes <- list()
    for (j in names(DEG_results)){
        for (i in 1:nrow((DEG_results[[j]]))){
            DEG_results[[j]]$rank[i] <- -log10(DEG_results[[j]]$p_val_adj[i]) * sign(DEG_results[[j]]$avg_log2FC[i])
        }
        ranked_genes[[j]] <- DEG_results[[j]]$rank
        names(ranked_genes[[j]]) <- DEG_results[[j]]$gene
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