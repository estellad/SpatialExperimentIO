#' @rdname readSeqfishSXE
#' 
#' @title Load data from a Spatial Genomics seqFISH experiment
#' 
#' @description
#' Creates a \code{\link{SpatialExperiment}} from the downloaded seqFISH   
#' directory for Spatial Genomics seqFISH spatial gene expression data.
#'
#' @param dirName a directory path to seqFISH download that contains files of interest.
#' @param returnType option of \code{"SPE"} or \code{"SCE"}, stands for 
#' \code{SpatialExperiment} or \code{SingleCellExperiment} object. Default value \code{"SPE"}
#' @param countMatPattern a filename pattern for the count matrix. Default value is 
#' \code{"CellxGene.csv"}, and there is no need to change.
#' @param metaDataPattern a filename pattern for the metadata .csv file that 
#' contains spatial coords. Default value is \code{"CellCoordinates.csv"}, and 
#' there is no need to change.
#' @param coordNames a vector of two strings specify the spatial coord names. 
#' Default value is \code{c("center_x", "center_y")}, and there is no need to change.
#' 
#' @details
#' The constructor assumes the downloaded seqFISH count matrix and metadata in 
#' the same folder with the following structure:
#' seqFISH_folder/ \cr
#' · | — *_CellxGene.csv \cr
#' · | — *_CellCoordinates.csv \cr
#' 
#'
#' @return  a \code{\link{SpatialExperiment}} or a \code{\link{SingleCellExperiment}} object 
#' @export
#' 
#' @author Yixing Estella Dong
#'
#' @examples
#' # A relatively small data download can be from:
#' # https://spatialgenomics.com/data/#kidney-data
#' 
#' 
#' # A mock counts and mock metadata with spatial location generated for a 9 genes by 
#' # 13 cells object is in /extdata: 
#' 
#' seqfpath <- system.file(
#'   file.path("extdata", "seqFISH_small"),
#'   package = "SpatialExperimentIO")
#'   
#' list.files(seqfpath)
#' 
#' # One of the following depending on your output (`SPE` or `SCE`) requirement.
#' seqf_spe <- readSeqfishSXE(dirName = seqfpath)
#' seqf_sce <- readSeqfishSXE(dirName = seqfpath, returnType = "SCE")
#' 
#' 
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment rowData counts colData
#' @importFrom methods as
#' @importFrom utils read.csv
readSeqfishSXE <- function(dirName = dirName, 
                           returnType = c("SPE", "SCE"),
                           countMatPattern = "CellxGene.csv", 
                           metaDataPattern = "CellCoordinates.csv", 
                           coordNames = c("center_x", "center_y")){
  
  returnType <- match.arg(returnType)
  
  tech <- "Seqfish"
  
  # Sanity checks
  countmat_file <- .sanityCheck(tech, filetype = "count matrix", expectfilename = "`CellxGene.csv`", dirName = dirName, filepatternvar = countMatPattern)
  metadata_file <- .sanityCheck(tech, filetype = "metadata", expectfilename = "`CellCoordinates.csv`", dirName = dirName, filepatternvar = metaDataPattern)
  
  # Read in 
  countmat <- read.csv(countmat_file)
  names(countmat)[names(countmat) == "X"] <- "cell"
  countmat <- countmat[order(countmat$cell), ]
  
  metadata <- read.csv(metadata_file)
  names(metadata)[names(metadata) == "label"] <- "cell"
  metadata <- metadata[order(metadata$cell), ]
  
  # Count matrix 
  countmat <- t(subset(countmat, select = -cell))
  
  # rowData (does not exist)
  
  # metadata
  metadata <- subset(metadata, select = -cell)
  
  if(!all(coordNames %in% colnames(metadata))){
    stop("`coordNames` not in columns of `metaDataPattern`. For seqFISH, expect c('center_x', 'center_y') in the columns of the metadata 'cell_metadata.csv'. " )
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
