#' @rdname readXeniumSXE
#' 
#' @title Load data from a 10x Geonomics Xenium experiment
#' 
#' @description
#' Creates a \code{\link{SpatialExperiment}} from the downloaded unzipped Xenium 
#' Output Bundle directory for 10x Genomics Xenium spatial gene expression data.
#'
#' @param dirname a directory path to Xenium Output Bundle download that contains 
#' files of interest.
#' @param return_type option of \code{"SPE"} or \code{"SCE"}, stands for 
#' \code{SpatialExperiment} or \code{SingleCellExperiment} object. Default value \code{"SPE"}
#' @param countfname a folder directory or the h5 file for the count matrix. 
#' Default value is \code{"cell_feature_matrix.h5"}, alternative value is 
#' \code{"cell_feature_matrix"} that takes a bit longer. The count matrix is 
#' read in and stored in a \code{SingleCellExperiment} object, using 
#' \code{DropletUtils::read10xCounts()}
#' @param coordfpattern a filename pattern of the zipped .csv file that 
#' contains spatial coords. Default value is \code{"cells.csv.gz"}, and there is no 
#' need to change.
#' @param coord_names a vector of two strings specify the spatial coord names. 
#' Default value is \code{c("x_centroid", "y_centroid")}, and there is no need to change.
#' 
#'
#' @details
#' The constructor assumes the downloaded unzipped Xenium Output Bundle has the 
#' following structure, with mandatory file of cells.csv.gz and either folder 
#' /cell_feature_matrix or .h5 file cell_feature_matrix.h5:
#' Xenium_unzipped \cr
#' · | — cell_feature_matrix.h5 \cr
#' · | — cell_feature_matrix \cr
#' · · | - barcodes.tsv.gz \cr
#' · · | - features.tsv.gz \cr
#' · · | - matrix.mtx.gz \cr
#' · | — cells.csv.gz \cr
#'
#' @return  a \code{\link{SpatialExperiment}} or a \code{\link{SingleCellExperiment}} object 
#' @export
#'
#' @author Estella Yixing Dong
#'
#' @examples
#' \dontrun{
#' # Data download is from: 
#' # https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_
#' # Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
#' 
#' xepath <- system.file(
#'   file.path("extdata", "10xXenium"),
#'   package = "SpatialExperimentIO")
#'   
#' list.files(xepath)
#' 
#' # One of the following depending on your input (.h5 or folder) and output 
#' # (`SPE` or `SCE`) requirement.
#' xe_spe <- readXeniumSXE(dirname = xepath)
#' xe_spe <- readXeniumSXE(dirname = xepath, countfname = "cell_feature_matrix")
#' xe_sce <- readXeniumSXE(dirname = xepath, return_type = "SCE")
#' 
#' # Subset to no control genes, and the same needed for `xe_sce` if read in as 
#' # `SCE`.
#' xe_spe <- xe_spe[rowData(xe_spe)$Type == "Gene Expression"]
#'
#' }
#' @importFrom DropletUtils read10xCounts
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
readXeniumSXE <- function(dirname, 
                          return_type = "SPE",
                          countfname = "cell_feature_matrix.h5",
                          coordfpattern = "cells.csv.gz", 
                          coord_names = c("x_centroid", "y_centroid")){
  
  if(!return_type %in% c("SPE", "SCE")){
    stop("'return_type' must be one of c('SPE', 'SCE')")
  }
  
  countfpath <- file.path(dirname, countfname)
  coord_file <- file.path(dirname, list.files(dirname, coordfpattern))
  
  # Count matrix + rowData
  sce <- DropletUtils::read10xCounts(countfpath, col.names = TRUE)
  
  # Spatial and colData
  colData <- read.csv(gzfile(coord_file), header = TRUE)
  
  if(return_type == "SPE"){
    # construct 'SpatialExperiment'
    sxe <- SpatialExperiment::SpatialExperiment(
      assays = assays(sce),
      rowData = rowData(sce),
      colData = colData,
      spatialCoordsNames = coord_names
    )
  }else if(return_type == "SCE"){
    # construct 'SingleCellExperiment'
    rownames(colData) <- colnames(sce)
    colData(sce) <- as(colData, "DFrame")
    sxe <- sce
  }
  
  if(class(counts(sxe)) != "dgCMatrix"){counts(sxe) <- as(counts(sxe), "dgCMatrix")}
  
  return(sxe)
}