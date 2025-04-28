#' @rdname readStarmapplusSXE
#' 
#' @title Load data from a STARmap PLUS experiment
#' 
#' @description
#' Creates a \code{\link{SpatialExperiment}} from the downloaded STARmap PLUS count
#' matrix.csv and metadata.csv
#'
#' @param dirName a directory path to STARmap PLUS download that contains files of interest.
#' @param returnType option of \code{"SPE"} or \code{"SCE"}, stands for 
#' \code{SpatialExperiment} or \code{SingleCellExperiment} object. Default value \code{"SPE"}
#' @param countMatPattern a filename pattern for the count matrix. Default value is 
#' \code{"raw_expression_pd.csv"}, and there is no need to change.
#' @param metaDataPattern a filename pattern for the metadata .csv file that 
#' contains spatial coords. Default value is \code{"spatial.csv"}, and 
#' there is no need to change.
#' @param coordNames a vector of three strings specify the spatial coord names. 
#' Default value is \code{c("X", "Y", "Z")}, and there is no need to change.
#' 
#' @details
#' The constructor assumes the downloaded unzipped STARmap PLUS folder has the following
#' structure, with two mandatory files:
#' STARmap_PLUS_download/ \cr
#' · | — *raw_expression_pd.csv \cr
#' · | — *spatial.csv \cr
#' 
#'
#' @return  a \code{\link{SpatialExperiment}} or a \code{\link{SingleCellExperiment}} object 
#' @export
#' 
#' @author Yixing Estella Dong
#'
#' @examples
#' # A relatively small data download can be from:
#' # https://zenodo.org/records/8327576
#' 
#' 
#' # A mock counts and mock metadata with spatial location generated for a 8 genes by 
#' # 9 cells object is in /extdata: 
#' 
#' starpath <- system.file(
#'   file.path("extdata", "STARmapPLUS_small"),
#'   package = "SpatialExperimentIO")
#' 
#' list.files(starpath)
#' 
#' # One of the following depending on your output (`SPE` or `SCE`) requirement.
#' star_spe <- readStarmapplusSXE(dirName = starpath)
#' star_sce <- readStarmapplusSXE(dirName = starpath, returnType = "SCE")
#' 
#' 
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment rowData counts colData
#' @importFrom methods as
#' @importFrom utils read.csv
readStarmapplusSXE <- function(dirName = dirName, 
                               returnType = c("SPE", "SCE"),
                               countMatPattern = "raw_expression_pd.csv", 
                               metaDataPattern = "spatial.csv", 
                               coordNames = c("X", "Y", "Z")){
  
  returnType <- match.arg(returnType)
  
  tech <- "Starmapplus"
  
  # Sanity checks
  countmat_file <- .sanityCheck(tech, filetype = "count matrix", expectfilename = "`raw_expression_pd.csv`", dirName = dirName, filepatternvar = countMatPattern)
  metadata_file <- .sanityCheck(tech, filetype = "metadata", expectfilename = "`spatial.csv`", dirName = dirName, filepatternvar = metaDataPattern)
  
  # Read in 
  countmat <- read.csv(countmat_file)
  metadata <- read.csv(metadata_file)
  
  # Count matrix
  rownames(countmat) <- countmat$GENE
  countmat <- subset(countmat, select = -GENE)  
  countmat <- as.matrix(countmat)

  # rowData (does not exist)
  
  # metadata
  rownames(metadata) <- metadata$NAME
  metadata <- metadata[rownames(metadata) != "TYPE", ]
  metadata$X <- as.numeric(metadata$X)
  metadata$Y <- as.numeric(metadata$Y)
  metadata$Z <- as.numeric(metadata$Z)
  
  if(!all(coordNames %in% colnames(metadata))){
    stop("`coordNames` not in columns of `metaDataPattern`. For STARmap PLUS, expect c('X', 'Y', 'Z') in the columns of the metadata 'spatial.csv'. " )
  }
  
  if(returnType == "SPE"){
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
