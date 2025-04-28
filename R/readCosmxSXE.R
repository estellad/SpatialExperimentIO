#' @rdname readCosmxSXE
#' 
#' @title Load data from a Nanostring CosMx experiment
#' 
#' @description
#' Creates a \code{\link{SpatialExperiment}} from the downloaded unzipped CosMx  
#' directory for Nanostring CosMx spatial gene expression data.
#'
#' @param dirName a directory path to CosMx download that contains files of interest.
#' @param returnType option of \code{"SPE"} or \code{"SCE"}, stands for 
#' \code{SpatialExperiment} or \code{SingleCellExperiment} object. Default value \code{"SPE"}
#' @param countMatPattern a filename pattern for the count matrix. Default value is 
#' \code{"exprMat_file.csv"}, and there is no need to change.
#' @param metaDataPattern a filename pattern for the metadata .csv file that 
#' contains spatial coords. Default value is \code{"metadata_file.csv"}, and 
#' there is no need to change.
#' @param coordNames a vector of two strings specify the spatial coord names. 
#' Default value is \code{c("CenterX_global_px", "CenterY_global_px")}, and 
#' there is no need to change.
#' @param addFovPos to read in fov_position_list.csv and add the data frame to 
#' \code{metadata(sxe)$fov_positions} or not. Default is TRUE.
#' @param fovPosPattern .csv pattern of fov_position_list.csv files in the raw download. 
#' Default value is \code{"fov_positions_file.csv"}.
#' @param altExps gene names contains these strings will be moved to \code{altExps(sxe)} 
#' as separate sxe-s. Default is \code{c("NegPrb", "Negative", "SystemControl", "FalseCode")}. 
#' @param addParquetPaths to add parquet paths to \code{metadata(sxe)} or not. If TRUE, 
#' transcripts and polygon .csv files will be converted to .parquet, and the paths 
#' will be added. If, for instance, no polygon file is available, and only transcript 
#' file is available, please set this argument to TRUE and adjust `addPolygon = FALSE` 
#' in the \code{...} argument. Default is TRUE.
#' @param ... extra parameters to pass to \code{\link{addParquetPathsCosmx}()}, including 
#' `addTx`, `txMetaNames`, `txPattern`, `addPolygon`, `polygonMetaNames`,
#' `polygonPattern`.
#' 
#' @details
#' The constructor assumes the downloaded unzipped CosMx folder has the following
#' structure, with two mandatory files:
#' CosMx_unzipped/optional_default_folder/ \cr
#' · | — *_exprMat_file.csv \cr
#' · | — *_metadata_file.csv \cr
#' 
#' Optional files to add to the metadata() as a list of paths (will be converted to parquet):
#' · | — *_fov_positions_file.csv \cr
#' · | — *_tx_file.csv \cr
#' · | — *_polygons.csv \cr
#' If no optional files, need to set `addFovPos = FALSE` and `addParquetPaths = FALSE`.
#' If only one of `*_tx_file.csv` or `*_polygons.csv` exists, set `addParquetPaths = TRUE` but 
#' set the not available `addTx` or `addPolygon` to `FALSE`. 
#' See \code{\link{addParquetPathsCosmx}()}
#'
#' @return  a \code{\link{SpatialExperiment}} or a \code{\link{SingleCellExperiment}} object 
#' @export
#' 
#' @author Yixing Estella Dong
#'
#' @examples
#' # A relatively small data download can be from:
#' # https://nanostring.com/resources/smi-ffpe-dataset-lung9-rep1-data/
#' 
#' 
#' # A mock counts and mock metadata with spatial location generated for a 8 genes by 
#' # 9 cells object is in /extdata: 
#' 
#' cospath <- system.file(
#'   file.path("extdata", "CosMx_small"),
#'   package = "SpatialExperimentIO")
#'   
#' list.files(cospath)
#' 
#' # One of the following depending on your output (`SPE` or `SCE`) requirement.
#' cos_spe <- readCosmxSXE(dirName = cospath, addPolygon = FALSE)
#' cos_sce <- readCosmxSXE(dirName = cospath, returnType = "SCE", addPolygon = FALSE)
#' cos_spe <- readCosmxSXE(dirName = cospath, addParquetPaths = FALSE)
#' 
#' 
#' @importFrom SpatialExperiment SpatialExperiment 
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SingleCellExperiment SingleCellExperiment rowData counts colData altExps altExps<-
#' @importFrom methods as
#' @importFrom data.table fread
#' @importFrom arrow write_parquet
#' 
readCosmxSXE <- function(dirName = dirName, 
                         returnType = c("SPE", "SCE"),
                         countMatPattern = "exprMat_file.csv", 
                         metaDataPattern = "metadata_file.csv", 
                         coordNames = c("CenterX_global_px", "CenterY_global_px"),
                         addFovPos = TRUE,
                         fovPosPattern = "fov_positions_file.csv",
                         altExps = c("NegPrb", "Negative", "SystemControl", "FalseCode"),
                         addParquetPaths = TRUE,
                         ...){
  
  returnType <- match.arg(returnType)
  tech <- "CosMx"
  
  # Sanity checks
  countmat_file <- .sanityCheck(tech, filetype = "count matrix", expectfilename = "`exprMat_file.csv`", 
                                dirName = dirName, filepatternvar = countMatPattern)
  metadata_file <- .sanityCheck(tech, filetype = "metadata", expectfilename = "`metadata_file.csv`", 
                                dirName = dirName, filepatternvar = metaDataPattern)
  if(addFovPos) fov_pos_file <- .sanityCheck(tech, filetype = "fov position", expectfilename = "`fov_positions_file.csv`", 
                                              dirName = dirName, filepatternvar = fovPosPattern)
  
  # Read in 
  countmat <- as.data.frame(fread(countmat_file))
  metadata <- as.data.frame(fread(metadata_file))
  
  overlap_cols <- intersect(colnames(countmat), colnames(metadata))
  
  # Count matrix   
  counts_ <- merge(countmat, metadata[, overlap_cols])
  counts_ <- counts_[, !(colnames(counts_) %in% overlap_cols)]
  counts_ <- t(counts_)
  
  # rowData (does not exist)
  
  # metadata
  metadata <- merge(metadata, countmat[, overlap_cols])

  if(!all(coordNames %in% colnames(metadata))){
    stop("`coordNames` not in columns of `metaDataPattern`. For CosMx, expect 
         c('CenterX_global_px', 'CenterY_global_px') in the columns of the metadata 'metadata_file.csv'. " )
  }
  
  colnames(counts_) <- rownames(metadata) <- seq_len(ncol(counts_))
  
  if(returnType == "SPE"){
    sxe <- SpatialExperiment(
      assays = list(counts = as(counts_, "dgCMatrix")),
      # rowData = rowData,
      colData = metadata,
      spatialCoordsNames = coordNames)
    #altExps=alt)
  }else if(returnType == "SCE"){
    sxe <- SingleCellExperiment(
      assays = list(counts = as(counts_, "dgCMatrix")),
      colData = metadata
      #altExps=alt
    )
  }
  
  if(!is.null(altExps)){
    names(altExps) <- altExps
    idx <- lapply(altExps, grep, rownames(sxe))     # get features matching input pattern(s)
    idx <- idx[vapply(idx, length, numeric(1)) > 0] # drops elements without match
    if(length(idx)){                                # skip if there aren't any matches
      alt <- lapply(idx, \(.) sxe[., ])
      altExps(sxe)[names(alt)] <- alt
      sxe <- sxe[-unlist(idx), ]
    }
  }
  
  # Add Parquet Paths
  if(addParquetPaths) sxe <- addParquetPathsCosmx(sxe, dirName, ...)
  
  # Read in and add FovPos to metadata()
  if(addFovPos){
    fovpos <- as.data.frame(fread(fov_pos_file))
    metadata(sxe)$fov_positions <- fovpos
  }
  
  return(sxe)
}

