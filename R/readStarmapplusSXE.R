#' @rdname readStarmapplusSXE
#' 
#' @title Load data from a STARmap PLUS experiment
#' 
#' @description
#' Creates a \code{\link{SpatialExperiment}} from the downloaded STARmap PLUS count
#' matrix.csv and metadata.csv
#'
#' @param dirname a directory path to STARmap PLUS download that contains files of interest.
#' @param return_type option of \code{"SPE"} or \code{"SCE"}, stands for 
#' \code{SpatialExperiment} or \code{SingleCellExperiment} object. Default value \code{"SPE"}
#' @param countmatfpattern a filename pattern for the count matrix. Default value is 
#' \code{"raw_expression_pd.csv"}, and there is no need to change.
#' @param metadatafpattern a filename pattern for the metadata .csv file that 
#' contains spatial coords. Default value is \code{"spatial.csv"}, and 
#' there is no need to change.
#' @param coord_names a vector of three strings specify the spatial coord names. 
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
#' \dontrun{
#' A relatively small data download can be from:
#' https://zenodo.org/records/8327576
#' 
#' A mock counts and mock metadata with spatial location generated for a 8 genes by 
#' 9 cells object is in /extdata: 
#' 
#' starpath <- system.file(
#'   file.path("extdata", "STARmapPLUS_small"),
#'   package = "SpatialExperimentIO")
#' 
#' list.files(starpath)
#' 
#' # One of the following depending on your output (`SPE` or `SCE`) requirement.
#' star_spe <- readStarmapplusSXE(dirname = starpath)
#' star_sce <- readStarmapplusSXE(dirname = starpath, return_type = "SCE")
#' 
#' }
#' 
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment rowData counts colData
#' @importFrom methods as
#' @importFrom utils read.csv
readStarmapplusSXE <- function(dirname = dirname, 
                               return_type = "SPE",
                               countmatfpattern = "raw_expression_pd.csv", 
                               metadatafpattern = "spatial.csv", 
                               coord_names = c("X", "Y", "Z")){
  
  if(!return_type %in% c("SPE", "SCE")){
    stop("'return_type' must be one of c('SPE', 'SCE')")
  }
  
  ## Metadata sanity check 
  if(!any(file.exists(file.path(dirname, list.files(dirname, metadatafpattern))))){
    stop("STARmap PLUS metadata file does not exist in the directory. Expect 'spatial.csv' in `dirname`")
  }
  
  metadata_file <- file.path(dirname, list.files(dirname, metadatafpattern))
  if(length(metadata_file) > 1){
    stop("More than one metadata file possible with the provided pattern `metadatafpattern`")
  }
  
  ## Count matrix sanity check
  if(!any(file.exists(file.path(dirname, list.files(dirname, countmatfpattern))))){
    stop("STARmap PLUS count matrix does not exist in the directory. Expect 'raw_expression_pd.csv' in `dirname`")
  }
  
  countmat_file <- file.path(dirname, list.files(dirname, countmatfpattern))
  if(length(countmat_file) > 1){
    stop("More than one count matrix file possible with the provided pattern `countmatfpattern`")
  }
  
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
  
  if(!all(coord_names %in% colnames(metadata))){
    stop("`coord_names` not in columns of `metadatafpattern`. For STARmap PLUS, expect c('X', 'Y', 'Z') in the columns of the metadata 'spatial.csv'. " )
  }
  
  if(return_type == "SPE"){
    sxe <- SpatialExperiment::SpatialExperiment(
      assays = list(counts = countmat),
      # rowData = rowData,
      colData = metadata,
      spatialCoordsNames = coord_names)
  }else if(return_type == "SCE"){
    # construct 'SingleCellExperiment'
    sxe <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = countmat),
      colData = metadata
    )
  }
  
  if(any(class(counts(sxe)) != "dgCMatrix")){counts(sxe) <- as(counts(sxe), "dgCMatrix")}
  
  return(sxe)
}
