#' @rdname readXeniumSXE
#' 
#' @title Load data from a 10x Geonomics Xenium experiment
#' 
#' @description
#' Creates a \code{\link{SpatialExperiment}} from the downloaded unzipped Xenium 
#' Output Bundle directory for 10x Genomics Xenium spatial gene expression data.
#'
#' @param dirName a directory path to Xenium Output Bundle download that contains 
#' files of interest.
#' @param returnType option of \code{"SPE"} or \code{"SCE"}, stands for 
#' \code{SpatialExperiment} or \code{SingleCellExperiment} object. Default value \code{"SPE"}
#' @param countMatPattern a folder directory or the h5 file pattern for the count matrix. 
#' Default value is \code{"cell_feature_matrix.h5"}, alternative value is 
#' \code{"cell_feature_matrix"} that takes a bit longer. The count matrix is 
#' read in and stored in a \code{SingleCellExperiment} object, using 
#' \code{DropletUtils::read10xCounts()}
#' @param metaDataPattern a filename pattern of the zipped .csv file that 
#' contains spatial coords. Default value is \code{"cells.csv.gz"}, and there is no 
#' need to change.
#' @param coordNames a vector of two strings specify the spatial coord names. 
#' Default value is \code{c("x_centroid", "y_centroid")}, and there is no need to change.
#' 
#' @param addExperimentXenium to add experiment.xenium parameters to \code{metadata(sxe)} or not. 
#' Default value is TRUE. 
#' @param altExps gene names contains these strings will be moved to \code{altExps(sxe)} as separate sxe-s. 
#' Default is \code{c("NegControlProbe", "UnassignedCodeword", "NegControlCodeword", "antisense", "BLANK")}. 
#' @param addParquetPaths to add parquet paths to \code{metadata(sxe)} or not. If TRUE, 
#' transcripts, cell_boundaries, and nucleus_boundaries .parquet paths will be added to `metadata()`. 
#' If, for instance, no cell_boundaries file is available, and transcript and nucleus_boundaries
#' files are available, please set this argument to TRUE and adjust `addCellBound = FALSE` 
#' in the \code{...} argument. Default is TRUE.
#' @param ... extra parameters to pass to \code{\link{addParquetPathsXenium}()}, including 
#' `addTx`, `txMetaNames`, `txPattern`, `addCellBound`, `cellBoundMetaNames`, 
#' `cellBoundPattern`, `addNucBound`, `NucBoundMetaNames`, `NucBoundPattern`.
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
#' · | — cells.parquet \cr
#'
#' Optional files to add to the metadata() as a list of paths (will be converted to parquet):
#' · | — transcripts.parquet \cr
#' · | — cell_boundaries.parquet \cr
#' · | — nucleus_boundaries.parquet \cr
#' · | — experiment.xenium \cr
#' See addParquetPathsXenium()
#'
#' @return  a \code{\link{SpatialExperiment}} or a \code{\link{SingleCellExperiment}} object 
#' @export
#'
#' @author Yixing Estella Dong
#'
#' @examples
#' # A relatively small data set is the Xenium mouse brain data that can be 
#' # downloaded from 10X website.
#' 
#' # A mock .h5 and mock metadata with spatial location generated for a 4 genes by 
#' # 6 cells object is in /extdata: 
#' 
#' xepath <- system.file(
#'   file.path("extdata", "Xenium_small"),
#'   package = "SpatialExperimentIO")
#'   
#' list.files(xepath)
#' 
#' # One of the following depending on your input (.h5 or folder) and output 
#' # (`SPE` or `SCE`) requirement.
#' xe_spe <- readXeniumSXE(dirName = xepath)
#' \dontrun{
#' xe_spe <- readXeniumSXE(dirName = xepath, countMatPattern = "cell_feature_matrix")
#' }
#' xe_sce <- readXeniumSXE(dirName = xepath, returnType = "SCE")
#' 
#' xe_spe <- readXeniumSXE(dirName = xepath, addParquetPaths = TRUE)
#' xe_spe <- readXeniumSXE(dirName = xepath, addParquetPaths = TRUE, addNucBound = FALSE)
#'
#' @importFrom DropletUtils read10xCounts
#' @importFrom SpatialExperiment SpatialExperiment 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SingleCellExperiment SingleCellExperiment rowData counts colData altExps altExps<-
#' @importFrom methods as
#' @importFrom arrow read_parquet read_json_arrow
#' 
readXeniumSXE <- function(dirName, 
                          returnType = c("SPE", "SCE"),
                          countMatPattern = c("cell_feature_matrix.h5", "cell_feature_matrix"),
                          metaDataPattern = c("cells.parquet", "cells.csv.gz"), 
                          coordNames = c("x_centroid", "y_centroid"), 
                          addExperimentXenium = TRUE,
                          altExps = c("NegControlProbe", "UnassignedCodeword", 
                                      "NegControlCodeword", "antisense", "BLANK"),
                          addParquetPaths = TRUE,
                          ...){

  returnType <- match.arg(returnType)
  countMatPattern <- match.arg(countMatPattern)
  metaDataPattern <- match.arg(metaDataPattern)
  tech <- "Xenium"
  
  # Sanity checks
  countmat_file <- .sanityCheck(tech, filetype = "count matrix", expectfilename = "`cell_feature_matrix.h5`", 
                                dirName = dirName, filepatternvar = countMatPattern)
  metadata_file <- .sanityCheck(tech, filetype = "metadata", expectfilename = "`cells.parquet` or `cells.csv.gz`", 
                                dirName = dirName, filepatternvar = metaDataPattern)
  if(addExperimentXenium) expxe_file <- .sanityCheck(tech, filetype = "experiment.xenium", expectfilename = "`experiment.xenium`", 
                                                     dirName = dirName, filepatternvar = "experiment.xenium")
  
  # folder 
  if(!grepl(".h5", countMatPattern)){
    folders <- list.files(dirName, countMatPattern)[!grepl(".h5", list.files(dirName, countMatPattern))]
    countmat_file <- file.path(dirName, folders)
    
    if(length(dir.exists(countmat_file)) > 1){
      stop("More than one count matrix folder possible with the provided pattern `countMatPattern`")
    }
    
    if(!all(c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz") %in% list.files(countmat_file))){
      stop("For Xenium with count matrix directory input, expect '/cell_feature_matrix' 
           folder contains files 'barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz'")
    }
    
  }

  # Count matrix + rowData
  sce <- DropletUtils::read10xCounts(countmat_file, col.names = TRUE)
  
  # Spatial and metadata
  if(grepl(".parquet", metaDataPattern)){
    metadata <- as.data.frame(read_parquet(metadata_file)) 
  }else if(grepl(".csv.gz", metaDataPattern)){
    metadata <- read.csv(gzfile(metadata_file), header = TRUE)
  }else{
    stop("Expect `.parquet` or `.csv.gz` format for metadata. Please convert and 
         save file needed first with e.g. arrow::write_parquet() or R.utils::gzip(). ")
  }

  if(!all(coordNames %in% colnames(metadata))){
    stop("`coordNames` not in columns of `metaDataPattern`. For Xenium, expect 
         c('x_centroid', 'y_centroid') in the columns of the metadata 'cells.csv.gz'. " )
  }
  
  if(returnType == "SPE"){
    sxe <- SpatialExperiment(
      assays = list(counts = as(counts(sce), "dgCMatrix")),
      rowData = rowData(sce),
      colData = metadata,
      spatialCoordsNames = coordNames
    )
  }else if(returnType == "SCE"){
    sxe <- SingleCellExperiment(
      assays = list(counts = as(counts(sce), "dgCMatrix")),
      rowData = rowData(sce),
      colData = metadata
    )
  }
  rownames(sxe) <- rowData(sxe)$Symbol
  
  if(!is.null(altExps)){
    names(altExps) <- altExps
    idx <- lapply(altExps, grep, rownames(sxe))
    idx <- idx[vapply(idx, length, numeric(1)) > 0]
    if(length(idx)){                              
      alt <- lapply(idx, \(.) sxe[., ]) 
      altExps(sxe)[names(alt)] <- alt
      sxe <- sxe[-unlist(idx), ]
    }
  }
  
  # Experiment Xenium
  if(addExperimentXenium) metadata(sxe)$experiment.xenium <- read_json_arrow(expxe_file)
  
  # Add Parquet Paths
  if(addParquetPaths) sxe <- addParquetPathsXenium(sxe, dirName, ...)
                                                        
  return(sxe)
}
