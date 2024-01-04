#' @rdname readStarmapplusSXE
#' 
#' @title Load data from a STARmap PLUS experiment
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
#' # https://zenodo.org/records/8327576
#' 
#' TODO: not this ~ ---------------------------------
#' starpath <- system.file(
#'   file.path("extdata", "STARmap_PLUS"),
#'   package = "SpatialExperimentIO")
#' --------------------------------------------------
#'   
#' TODO: But something like this! -------------------
#' starpath <- "path/to/a/defined/folder/"
#' 
#' eh <- ExperimentHub()
#' query(eh, c("scSpatialExperimentData", "STARmap_PLUS"))
#' 
#' count.dat <- eh[["EH7543"]]
#' write.csv(count.dat, file.path(starpath, "raw_expression_pd.csv")
#' col.dat <- eh[["EH7544"]]
#' write.csv(col.dat, file.path(starpath, "spatial.csv")
#' --------------------------------------------------
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
#' @importFrom SingleCellExperiment SingleCellExperiment
readStarmapplusSXE <- function(dirname = dirname, 
                               return_type = "SPE",
                               countmatfpattern = "raw_expression_pd.csv", 
                               metadatafpattern = "spatial.csv", 
                               coord_names = c("X", "Y", "Z")){
  
  if(!return_type %in% c("SPE", "SCE")){
    stop("'return_type' must be one of c('SPE', 'SCE')")
  }
  
  countmat_file <- file.path(dirname, list.files(dirname, countmatfpattern))
  metadata_file <- file.path(dirname, list.files(dirname, metadatafpattern))
  
  # Read in 
  countmat <- read.csv(countmat_file)
  metadata <- read.csv(metadata_file)
  
  # Count matrix
  rownames(stmpp_counts) <- stmpp_counts$GENE
  stmpp_counts <- subset(stmpp_counts, select = -GENE)  

  # rowData (does not exist)
  
  # colData
  rownames(stmpp_coord) <- stmpp_coord$NAME
  stmpp_coord <- stmpp_coord[rownames(stmpp_coord) != "TYPE", ]
  stmpp_coord$X <- as.numeric(stmpp_coord$X)
  stmpp_coord$Y <- as.numeric(stmpp_coord$Y)
  stmpp_coord$Z <- as.numeric(stmpp_coord$Z)
  stmpp_coord$Main_molecular_cell_type <- as.factor(stmpp_coord$Main_molecular_cell_type)
  stmpp_coord$Sub_molecular_cell_type <- as.factor(stmpp_coord$Sub_molecular_cell_type)
  stmpp_coord$Main_molecular_tissue_region <- as.factor(stmpp_coord$Main_molecular_tissue_region)
  stmpp_coord$Sub_molecular_tissue_region <- as.factor(stmpp_coord$Sub_molecular_tissue_region)
  stmpp_coord$Molecular_spatial_cell_type <- as.factor(stmpp_coord$Molecular_spatial_cell_type)
  
  if(return_type == "SPE"){
    sxe <- SpatialExperiment::SpatialExperiment(
      assays = list(counts = stmpp_counts),
      # rowData = rowData,
      colData = stmpp_coord,
      spatialCoordsNames = coord_names)
  }else if(return_type == "SCE"){
    # construct 'SingleCellExperiment'
    sxe <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = stmpp_counts),
      colData = stmpp_coord
    )
  }
  
  if(class(counts(sxe)) != "dgCMatrix"){counts(sxe) <- as(counts(sxe), "dgCMatrix")}
  
  return(sxe)
}