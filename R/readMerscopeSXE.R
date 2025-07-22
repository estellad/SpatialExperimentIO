#' @rdname readMerscopeSXE
#' 
#' @title Load data from a Vizgen MERSCOPE experiment
#' 
#' @description
#' Creates a \code{\link{SpatialExperiment}} from the downloaded MERSCOPE   
#' directory for Vizgen MERSCOPE spatial gene expression data.
#'
#' @param dirName a directory path to MERSCOPE download that contains files of interest.
#' @param returnType option of \code{"SPE"} or \code{"SCE"}, stands for 
#' \code{SpatialExperiment} or \code{SingleCellExperiment} object. Default value \code{"SPE"}
#' @param countMatPattern a filename pattern for the count matrix. Default value is 
#' \code{"cell_by_gene.csv"}, and there is no need to change.
#' @param metaDataPattern a filename pattern for the metadata .csv file that 
#' contains spatial coords. Default value is \code{"metadata_file.csv"}, and 
#' there is no need to change.
#' @param coordNames a vector of two strings specify the spatial coord names. 
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
#' # A relatively small data download can be from:
#' # https://console.cloud.google.com/storage/browser/vz-ffpe-showcase/
#' # HumanOvarianCancerPatient2Slice2?pageState=(%22StorageObjectListTable%22:
#' # (%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false
#' 
#' 
#' # A mock counts and mock metadata with spatial location generated for a 9 genes by 
#' # 8 cells object is in /extdata: 
#' 
#' merpath <- system.file(
#'   file.path("extdata", "MERSCOPE_small"),
#'   package = "SpatialExperimentIO")
#'   
#' list.files(merpath)
#' 
#' # One of the following depending on your output (`SPE` or `SCE`) requirement.
#' mer_spe <- readMerscopeSXE(dirName = merpath)
#' mer_sce <- readMerscopeSXE(dirName = merpath, returnType = "SCE")
#' 
#' 
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment rowData counts colData
#' @importFrom methods as
#' @importFrom utils read.csv
readMerscopeSXE <- function(dirName = dirName, 
                            returnType = c("SPE", "SCE"),
                            countMatPattern = "cell_by_gene.csv", 
                            metaDataPattern = "cell_metadata.csv", 
                            coordNames = c("center_x", "center_y")){
  
  returnType <- match.arg(returnType)
  tech <- "Merscope"
  
  # Sanity checks
  countmat_file <- .sanityCheck(tech, filetype = "count matrix", expectfilename = "`cell_by_gene.csv`", dirName = dirName, filepatternvar = countMatPattern)
  metadata_file <- .sanityCheck(tech, filetype = "metadata", expectfilename = "`cell_metadata.csv`", dirName = dirName, filepatternvar = metaDataPattern)
  
  # Read in 
  countmat <- read.csv(countmat_file)
  countmat <- countmat[order(countmat$cell), ]
  
  metadata <- data.table::fread(metadata_file)
  names(metadata)[names(metadata) %in% c("EntityID", "V1", "X")] <- "cell"
  metadata$cell <- as.character(metadata$cell)
  metadata <- metadata[order(metadata$cell), ]
  
  # Count matrix 
  rownames(countmat) <- countmat$cell + 1 # fix the indexing to start at 1
  countmat <- t(subset(countmat, select = -cell))
  
  # rowData (does not exist)
  
  # metadata
  rownames(metadata) <- metadata$cell + 1
  metadata <- subset(metadata, select = -cell)
  
  if(!all(coordNames %in% colnames(metadata))){
    stop("`coordNames` not in columns of `metaDataPattern`. For MERSCOPE, expect c('center_x', 'center_y') in the columns of the metadata 'cell_metadata.csv'. " )
  }
  
  if(returnType == "SPE"){
    # construct 'SpatialExperiment'
    sxe <- SpatialExperiment::SpatialExperiment(
      assays = list(counts = as(countmat, "dgCMatrix")),
      # rowData = rowData,
      colData = metadata,
      spatialCoordsNames = coordNames)
    }else if(returnType == "SCE"){
      # construct 'SingleCellExperiment'
      sxe <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = as(countmat, "dgCMatrix")),
        colData = metadata
      )
    }
  
  return(sxe)
  
}
