#' Save your file without overwriting it
#' @param object_name Name of your object you want to save
#' @param path_to_object Path where you want to store your object, needs to contain the name for your object
#' @param overwrite TRUE for overwriting, default is FALSE hello


#' @export

save_object <- function(object_name, path_to_object, overwrite = FALSE){
    if (file.exists(path_to_object) == FALSE){
        saveRDS(object_name, path_to_object)
    }
    else if (file.exists(path_to_object) == TRUE & overwrite == TRUE){
        saveRDS(object_name, path_to_object)
    }
    else {
        print("Already saved")
    }
}
