#
#
#
#' get_trafo_from_WSP
#'
#' This function loads a FlowJo WSP file, extracts the parameters for the 
#' transformation used for every channel and writes them to a table in 
#' csv format. 
#' @param file_name A WSP file (potentially including the path to it)
#' @param sample_name The sample name used to separate the different samples, 
#' defaults to 'samples'
#' @param out_prefix File name prefix for the data table with channel names 
#' as rows and parameters as columns. One file per sample. 
#'
#' @return Write a data table to file.
#' @export
#'
#' @examples
get_trafo_from_WSP <- function(file_name = "test_workspace.wsp",
                            sample_name = 'All Samples',
                            out_prefix = 'test_workspace_params'){
  library(flowWorkspace)
  library(CytoML)
    #open flowjo workspace file without loading data
  ws <- open_flowjo_xml(file = file_name)
  gs <- flowjo_to_gatingset(ws, name = sample_name, execute = FALSE)
  
  #get channel names by using the compensations function, then extract 
  #transformation parameters per sample
  comp_table <- gs_get_compensations(gs) 
  #get sample names from comp_table list
  sample_ids <- names(comp_table)
  chnnls <- list()
  attrib_params_all <- list()
  for (sample_id in sample_ids){
    chnnls[[sample_id]] <- names(comp_table[[sample_id]]@parameters)
    #get transformation parameters for each channel
    attribs_params <- list()
    for (chnnl in chnnls[[sample_id]]){
      trafo = gh_get_transformations(gs[[sample_id]], channel = chnnl)
      attribs_params[[chnnl]] <- attributes(trafo)#$parameters
    }
    attrib_params_all[[sample_id]] <- attribs_params
  }
  
  #write values of the attrib_params_all list into a Python readable format
  #it should be a table per sample with the channels as row names and the 
  #type and parameters as column names
  for (sample_id in sample_ids){
    tmp_table <- list()
    for (chnnl in chnnls[[sample_id]]){
      attrib_table <- attrib_params_all[[sample_id]][[chnnl]]
      tmp_table[[chnnl]] <- data.frame(t(rapply(sapply(attrib_table,c),c)))
      
    }
    
    #concatenate all entries in tmp_table as one data frame
    param_table <- data.frame(t(sapply(tmp_table, rbind)))
    param_table <- apply(param_table, 2, as.character)
    rownames(param_table) <- names(tmp_table)
    out_file <- paste0(out_prefix, '_', sample_id ,'.csv')
    
    write.csv(x=param_table, file=out_file)
  }
  
}






