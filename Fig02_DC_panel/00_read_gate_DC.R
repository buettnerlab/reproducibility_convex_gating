#project: convex gating
#author: Maren Buettner
#date: 18/01/2022
#
#description: Read gating information from WSP FlowJo files 
#data: DC panel

source('./../code/extract_trafo_from_wsp.R')

## set paths
#######
project_path <- './'
#check if project_path is accessible
list.files(project_path)
#set data_path
data_path <- paste0(project_path, 'data/Hofer_IMA/FCS_FJ_files/') 
wsp_file <- paste0(data_path, "V1163 1228 Mono DC2 slan.wsp")
out_prefix <- paste0(data_path,'workspace_hofer_params')
pop <- "root"

## read files
########
get_trafo_from_WSP(file_name = wsp_file, out_prefix =  out_prefix)
out_prefix <- paste0(data_path,'workspace_hofer_gate')
#read gating info
library(flowWorkspace)
library(CytoML)
library(magrittr)
library(Biobase)
#open flowjo workspace file with loading data
ws <- open_flowjo_xml(file = wsp_file)
gs <- flowjo_to_gatingset(ws, name = 'All Samples')
#get populations
gate_list <- gs_get_pop_paths(gs, path = "auto")

#get data
fs <- gs_pop_get_data(obj = gs, 
                      y = pop, 
                      inverse.transform = FALSE) %>% cytoset_to_flowSet()

filenames <- rownames(fs@phenoData)

NumBC <- length(fs)

FFdata <- NULL
#get all measured features (both -A and -H)
OrigNames <-fs[[1]]@parameters$name
#iterate over all samples
for (FFs in filenames){
  #get time of measurement (as proxy for event ID)
  FFa <- exprs(fs[[FFs]]$Time)
  
  # Add the gating info to feature matrix
  for (gate in gate_list) {
    
    FFa <- cbind(FFa, gh_pop_get_indices(gs[[FFs]], y = gate))
    colnames(FFa)[dim(FFa)[2]] <- gate
    
  }
  
  #remove file extension from sample name
  file_id <- unlist(strsplit(FFs, '.fcs', fixed=TRUE))[1]
  
  out_file <- paste0(out_prefix, '_', file_id ,'.csv')
  #write gating info to file
  write.csv(x=FFa, file=out_file)
}

