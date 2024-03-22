#' @rdname readMerscopeSXE
#' 
#' @title Load data from a Vizgen MERSCOPE experiment
#' 
#' @description
#' Creates a \code{\link{SpatialExperiment}} from the downloaded MERSCOPE   
#' directory for Vizgen MERSCOPE spatial gene expression data.
#'
#' @param dirname a directory path to MERSCOPE download that contains files of interest.
#' @param return_type option of \code{"SPE"} or \code{"SCE"}, stands for 
#' \code{SpatialExperiment} or \code{SingleCellExperiment} object. Default value \code{"SPE"}
#' @param countmatfpattern a filename pattern for the count matrix. Default value is 
#' \code{"cell_by_gene.csv"}, and there is no need to change.
#' @param metadatafpattern a filename pattern for the metadata .csv file that 
#' contains spatial coords. Default value is \code{"metadata_file.csv"}, and 
#' there is no need to change.
#' @param coord_names a vector of two strings specify the spatial coord names. 
#' Default value is \code{c("center_x", "center_y")}, and there is no need to change.
#' 
#' @details
#' The constructor assumes the downloaded MERSCOPE count matrix and metadata in 
#' the same folder with the following structure:
#' MERSCOPE_folder/ \cr
#' · | — cell_by_gene.csv \cr
#' · | — cell_metadata.csv \cr
#'
#' @return  a \code{\link{SpatialExperiment}} or a \code{\link{SingleCellExperiment}} object 
#' @export
#' 
#' @author Yixing Estella Dong
#'
#' @examples
#' \dontrun{
#' # A relatively small data download can be from:
#' https://console.cloud.google.com/storage/browser/vz-ffpe-showcase/
#' HumanOvarianCancerPatient2Slice2?pageState=(%22StorageObjectListTable%22:
#' (%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false
#' 
#' A mock counts and mock metadata with spatial location generated for a 9 genes by 
#' 8 cells object is in /extdata: 
#' 
#' merpath <- system.file(
#'   file.path("extdata", "MERSCOPE_small"),
#'   package = "SpatialExperimentIO")
#'   
#' list.files(merpath)
#' 
#' # One of the following depending on your output (`SPE` or `SCE`) requirement.
#' mer_spe <- readMerscopeSXE(dirname = merpath)
#' mer_sce <- readMerscopeSXE(dirname = merpath, return_type = "SCE")
#' 
#' }
#' 
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment rowData counts colData
#' @importFrom methods as
#' @importFrom utils read.csv
readMerscopeSXE <- function(dirname = dirname, 
                            return_type = "SPE",
                            countmatfpattern = "cell_by_gene.csv", 
                            metadatafpattern = "cell_metadata.csv", 
                            coord_names = c("center_x", "center_y")){
  
  if(!return_type %in% c("SPE", "SCE")){
    stop("'return_type' must be one of c('SPE', 'SCE')")
  }
  
  ## Metadata sanity check 
  if(!any(file.exists(file.path(dirname, list.files(dirname, metadatafpattern))))){
    stop("MERSCOPE metadata file does not exist in the directory. Expect 'cell_metadata.csv' in `dirname`")
  }
  
  metadata_file <- file.path(dirname, list.files(dirname, metadatafpattern))
  if(length(metadata_file) > 1){
    stop("More than one metadata file possible with the provided pattern `metadatafpattern`")
  }
  
  ## Count matrix sanity check
  if(!any(file.exists(file.path(dirname, list.files(dirname, countmatfpattern))))){
    stop("MERSCOPE count matrix does not exist in the directory. Expect 'cell_by_gene.csv' in `dirname`")
  }
  
  countmat_file <- file.path(dirname, list.files(dirname, countmatfpattern))
  if(length(countmat_file) > 1){
    stop("More than one count matrix file possible with the provided pattern `countmatfpattern`")
  }
  
  # Read in 
  countmat <- read.csv(countmat_file)
  countmat <- countmat[order(countmat$cell), ]
  
  metadata <- read.csv(metadata_file)
  names(metadata)[names(metadata) == "X"] <- "cell"
  metadata <- metadata[order(metadata$cell), ]
  
  # Count matrix 
  rownames(countmat) <- countmat$cell + 1 # fix the indexing to start at 1
  countmat <- t(subset(countmat, select = -cell))
  
  # rowData (does not exist)
  
  # metadata
  rownames(metadata) <- metadata$cell + 1
  metadata <- subset(metadata, select = -cell)
  
  if(!all(coord_names %in% colnames(metadata))){
    stop("`coord_names` not in columns of `metadatafpattern`. For MERSCOPE, expect c('center_x', 'center_y') in the columns of the metadata 'cell_metadata.csv'. " )
  }
  
  if(return_type == "SPE"){
    # construct 'SpatialExperiment'
    sxe <- SpatialExperiment::SpatialExperiment(
      assays = list(counts = as(countmat, "dgCMatrix")),
      # rowData = rowData,
      colData = metadata,
      spatialCoordsNames = coord_names)
    }else if(return_type == "SCE"){
      # construct 'SingleCellExperiment'
      sxe <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = as(countmat, "dgCMatrix")),
        colData = metadata
      )
    }
  
  # if(any(class(counts(sxe)) != "dgCMatrix")){counts(sxe) <- as(counts(sxe), "dgCMatrix")}
  
  return(sxe)
  
}
