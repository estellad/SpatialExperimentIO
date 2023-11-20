#' @rdname readCosmxSPE
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
#' @author Estella Yixing Dong
#'
#' @examples
#' \dontrun{
#' # Data download is from: 
#' # https://nanostring.com/resources/smi-ffpe-dataset-lung9-rep1-data/
#' 
#' # Here as an example, we have downsized the count matrix and meta data file 
#' # to only 453 randomly selected cells and 20 genes, and necessary columns 
#' # of the metadata. 
#' 
#' cospath <- system.file(
#'   file.path("extdata", "NanostringCosMx"),
#'   package = "SpatialExperiment")
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
#' @importFrom SingleCellExperiment SingleCellExperiment
readCosmxSXE <- function(dirname = dirname, 
                         return_type = "SPE",
                         countmatfpattern = "exprMat_file.csv", 
                         metadatafpattern = "metadata_file.csv", 
                         coord_names = c("CenterX_global_px",
                                         "CenterY_global_px")){
  
  if(!return_type %in% c("SPE", "SCE")){
    stop("'return_type' must be one of c('SPE', 'SCE')")
  }
  
  countmat_file <- file.path(dirname, list.files(dirname, countmatfpattern))
  metadata_file <- file.path(dirname, list.files(dirname, metadatafpattern))
  
  # Read in 
  countmat <- read.csv(countmat_file)
  metadata <- read.csv(metadata_file)
  
  # Count matrix   
  counts <- merge(countmat, metadata[, c("fov", "cell_ID")])
  counts <- subset(counts, select = -c(fov, cell_ID))
  counts <- t(counts)
  
  # rowData (does not exist)
  
  # colData
  colData <- merge(metadata, countmat[, c("fov", "cell_ID")])
  
  
  colnames(counts) <- rownames(colData) <- 1:ncol(counts)
  
  if(return_type == "SPE"){
  sxe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    # rowData = rowData,
    colData = colData,
    spatialCoordsNames = coord_names)
  }else if(return_type == "SCE"){
    # construct 'SingleCellExperiment'
    sxe <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = counts),
      colData = colData
    )
  }
  
  return(sxe)
}