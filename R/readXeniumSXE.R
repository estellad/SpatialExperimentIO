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
#' @param countmatfpattern a folder directory or the h5 file pattern for the count matrix. 
#' Default value is \code{"cell_feature_matrix.h5"}, alternative value is 
#' \code{"cell_feature_matrix"} that takes a bit longer. The count matrix is 
#' read in and stored in a \code{SingleCellExperiment} object, using 
#' \code{DropletUtils::read10xCounts()}
#' @param metadatafpattern a filename pattern of the zipped .csv file that 
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
#' @author Yixing Estella Dong
#'
#' @examples
#' \dontrun{
#' A relatively small data set is the Xenium mouse brain data that can be downloaded from 10X website.
#' 
#' A mock .h5 and mock metadata with spatial location generated for a 4 genes by 
#' 6 cells object is in /extdata: 
#' 
#' xepath <- system.file(
#'   file.path("extdata", "Xenium_small"),
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
#' @importFrom SingleCellExperiment SingleCellExperiment rowData counts colData
#' @importFrom methods as
#' @importFrom utils read.csv
readXeniumSXE <- function(dirname, 
                          return_type = "SPE",
                          countmatfpattern = "cell_feature_matrix.h5",
                          metadatafpattern = "cells.csv.gz", 
                          coord_names = c("x_centroid", "y_centroid")){
  
  if(!return_type %in% c("SPE", "SCE")){
    stop("'return_type' must be one of c('SPE', 'SCE')")
  }
  
  ## Metadata sanity check 
  if(!any(file.exists(file.path(dirname, list.files(dirname, metadatafpattern))))){
    stop("Xenium metadata file does not exist in the directory. Expect 'cells.csv.gz' in `dirname`")
  }
  
  metadata_file <- file.path(dirname, list.files(dirname, metadatafpattern))
  if(length(metadata_file) > 1){
    stop("More than one metadata file possible with the provided pattern `metadatafpattern`")
  }
  
  ## Count matrix sanity check
  if(!any(file.exists(file.path(dirname, list.files(dirname, countmatfpattern))))){
    stop("Xenium count matrix .h5 file or directory does not exist in the directory. Expect 'cell_feature_matrix.h5' or folder `/cell_feature_matrix` in `dirname`")
  }
  
  countmat_file <- file.path(dirname, list.files(dirname, countmatfpattern))
  
  # .h5 file 
  if(grepl(".h5", countmatfpattern)){
    if(length(countmat_file) > 1){
      stop("More than one count matrix .h5 file possible with the provided pattern `countmatfpattern`")
    }
  }
  
  # folder 
  if(!grepl(".h5", countmatfpattern)){
    folders <- list.files(dirname, countmatfpattern)[!grepl(".h5", list.files(dirname, countmatfpattern))]
    countmat_file <- file.path(dirname, folders)
    
    if(length(dir.exists(countmat_file)) > 1){
      stop("More than one count matrix folder possible with the provided pattern `countmatfpattern`")
    }
    
    if(!all(c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz") %in% list.files(countmat_file))){
      stop("For Xenium with count matrix directory input, expect '/cell_feature_matrix' folder contains files 'barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz'")
    }
    
  }

  
  # Count matrix + rowData
  sce <- DropletUtils::read10xCounts(countmat_file, col.names = TRUE)
  
  # Spatial and metadata
  metadata <- read.csv(gzfile(metadata_file), header = TRUE)
  
  if(!all(coord_names %in% colnames(metadata))){
    stop("`coord_names` not in columns of `metadatafpattern`. For Xenium, expect c('x_centroid', 'y_centroid') in the columns of the metadata 'cells.csv.gz'. " )
  }
  
  if(return_type == "SPE"){
    # construct 'SpatialExperiment'
    sxe <- SpatialExperiment::SpatialExperiment(
      assays = assays(sce),
      rowData = rowData(sce),
      colData = metadata,
      spatialCoordsNames = coord_names
    )
  }else if(return_type == "SCE"){
    # construct 'SingleCellExperiment'
    rownames(metadata) <- colnames(sce)
    colData(sce) <- as(metadata, "DFrame")
    sxe <- sce
  }
  
  if(any(class(counts(sxe)) != "dgCMatrix")){counts(sxe) <- as(counts(sxe), "dgCMatrix")}
  
  return(sxe)
}
