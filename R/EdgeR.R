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