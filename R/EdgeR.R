#' Functions for analaysing single cell data with EdgeR
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
                df[j,i] <- 1/4
            }
            else{
                df[j,i] <- 1/-(length(rownames(df))-4)
            }
        }
    }
    return(df)
    
}