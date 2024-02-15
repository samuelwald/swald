#' Create folder structure

#' @export 

create_project <- function(){
    
    # List of folders
    folders <- list(notebooks = "notebooks/",
                n_preprocessing = "notebooks/preprocessing",
                n_quality_control = "notebooks/quality_control",
                n_integration = "notebooks/integration",
                n_annotation = "notebooks/annotation",
                n_analysis = "notebooks/analysis",
                objects = "objects/",
                o_preprocessing = "objects/preprocessing",
                o_quality_control = "objects/quality_control",
                o_integration = "objects/integration",
                o_annotation = "objects/annotation",
                o_analysis = "objects/analysis",
                files = "files/",
                f_preprocessing = "files/preprocessing",
                f_quality_control = "files/quality_control",
                f_integration = "files/integration",
                f_annotation = "files/annotation",
                f_analysis = "files/analysis",
                f_uploads = "files/uploads",
                figures = "figures",
                fi_preprocessing = "figures/preprocessing",
                fi_quality_control = "figures/quality_control",
                fi_integration = "figures/integration",
                fi_annotation = "figures/annotation",
                fi_analysis = "figures/analysis",
                raw_data = "raw_data",
                fastq = "raw_data/FASTQ")
    
    # Create folders
    for (i in names(folders)){
        dir.create(paste(getwd(), folders[[i]], sep = "/"))
    }
}