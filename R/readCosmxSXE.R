#' @rdname readCosmxSXE
#' 
#' @title Load data from a Nanostring CosMx experiment
#' 
#' @description
#' Creates a \code{\link{SpatialExperiment}} from the downloaded unzipped CosMx  
#' directory for Nanostring CosMx spatial gene expression data.
#'
#' @param dirname a directory path to CosMx download that contains files of interest.
#' @param return_type option of \code{"SPE"} or \code{"SCE"}, stands for 
#' \code{SpatialExperiment} or \code{SingleCellExperiment} object. Default value \code{"SPE"}
#' @param countmatfpattern a filename pattern for the count matrix. Default value is 
#' \code{"exprMat_file.csv"}, and there is no need to change.
#' @param metadatafpattern a filename pattern for the metadata .csv file that 
#' contains spatial coords. Default value is \code{"metadata_file.csv"}, and 
#' there is no need to change.
#' @param coord_names a vector of two strings specify the spatial coord names. 
#' Default value is \code{c("CenterX_global_px", "CenterY_global_px")}, and 
#' there is no need to change.
#' 
#' @details
#' The constructor assumes the downloaded unzipped CosMx folder has the following
#' structure, with two mandatory files:
#' CosMx_unzipped/optional_default_folder/ \cr
#' · | — *_exprMat_file.csv \cr
#' · | — *_metadata_file.csv \cr
#' 
#'
#' @return  a \code{\link{SpatialExperiment}} or a \code{\link{SingleCellExperiment}} object 
#' @export
#' 
#' @author Yixing Estella Dong
#'
#' @examples
#' \dontrun{
#' # A relatively small data download can be from:
#' https://nanostring.com/resources/smi-ffpe-dataset-lung9-rep1-data/
#' 
#' A mock counts and mock metadata with spatial location generated for a 8 genes by 
#' 9 cells object is in /extdata: 
#' 
#' cospath <- system.file(
#'   file.path("extdata", "CosMx_small"),
#'   package = "SpatialExperimentIO")
#'   
#' list.files(cospath)
#' 
#' # One of the following depending on your output (`SPE` or `SCE`) requirement.
#' cos_spe <- readCosmxSXE(dirname = cospath)
#' cos_sce <- readCosmxSXE(dirname = cospath, return_type = "SCE")
#' 
#' }
#' 
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment rowData counts colData
#' @importFrom methods as
#' @importFrom utils read.csv
readCosmxSXE <- function(dirname = dirname, 
                         return_type = "SPE",
                         countmatfpattern = "exprMat_file.csv", 
                         metadatafpattern = "metadata_file.csv", 
                         coord_names = c("CenterX_global_px", "CenterY_global_px")){
  
  if(!return_type %in% c("SPE", "SCE")){
    stop("'return_type' must be one of c('SPE', 'SCE')")
  }
  
  ## Metadata sanity check 
  if(!any(file.exists(file.path(dirname, list.files(dirname, metadatafpattern))))){
    stop("CosMx metadata file does not exist in the directory. Expect 'metadata_file.csv' in `dirname`")
  }
  
  metadata_file <- file.path(dirname, list.files(dirname, metadatafpattern))
  if(length(metadata_file) > 1){
    stop("More than one metadata file possible with the provided pattern `metadatafpattern`")
  }
  
  ## Count matrix sanity check
  if(!any(file.exists(file.path(dirname, list.files(dirname, countmatfpattern))))){
    stop("CosMx count matrix does not exist in the directory. Expect 'exprMat_file.csv' in `dirname`")
  }
  
  countmat_file <- file.path(dirname, list.files(dirname, countmatfpattern))
  if(length(countmat_file) > 1){
    stop("More than one count matrix file possible with the provided pattern `countmatfpattern`")
  }

  # Read in 
  countmat <- read.csv(countmat_file)
  metadata <- read.csv(metadata_file)
  
  # Count matrix   
  countmat <- merge(countmat, metadata[, c("fov", "cell_ID")])
  countmat <- subset(countmat, select = -c(fov, cell_ID))
  countmat <- t(countmat)
  
  # rowData (does not exist)
  
  # metadata
  metadata <- merge(metadata, countmat[, c("fov", "cell_ID")])
  
  if(!all(coord_names %in% colnames(metadata))){
    stop("`coord_names` not in columns of `metadatafpattern`. For CosMx, expect c('CenterX_global_px', 'CenterY_global_px') in the columns of the metadata 'metadata_file.csv'. " )
  }
  
  colnames(countmat) <- rownames(metadata) <- 1:ncol(countmat)
  
  if(return_type == "SPE"){
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
