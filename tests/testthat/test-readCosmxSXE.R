dir <- system.file(
  file.path("extdata", "CosMx_small"), 
  package="SpatialExperimentIO")

test_that("example data folders uniquely contains needed files", {
  expect_true("lung_p9s1_exprMat_file.csv" %in% list.files(dir))
  expect_true("lung_p9s1_metadata_file.csv" %in% list.files(dir)) 
  
  expect_length(list.files(dir, "exprMat_file.csv"), 1)
  expect_length(list.files(dir, "metadata_file.csv"), 1)
})

test_that("needed files contains spatial columns of interest", {
  metadata <- read.csv(file.path(dir, "lung_p9s1_metadata_file.csv"))
  
  expect_true(all(c("CenterX_global_px", "CenterY_global_px") %in% colnames(metadata))) 
  expect_true(is.numeric(metadata$CenterX_global_px))
  expect_true(is.numeric(metadata$CenterY_global_px))
})

test_that("data are read correctly to SpatialExperiment class", {
  x <- readCosmxSXE(dirname = dir, 
                    return_type = "SPE",
                    countmatfpattern = "exprMat_file.csv", 
                    metadatafpattern = "metadata_file.csv", 
                    coord_names = c("CenterX_global_px",
                                    "CenterY_global_px"))
  
  expect_s4_class(x, "SpatialExperiment")
  expect_true(all(colnames(SpatialExperiment::spatialCoords(x)) == c("CenterX_global_px", "CenterY_global_px")))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(8, 9)))
  expect_s4_class(SingleCellExperiment::counts(x), "dgCMatrix")
})

test_that("data are read correctly to SingleCellExperiment class", {
  x <- readCosmxSXE(dirname = dir, 
                    return_type = "SCE",
                    countmatfpattern = "exprMat_file.csv",
                    metadatafpattern = "metadata_file.csv", 
                    coord_names = c("CenterX_global_px", 
                                    "CenterY_global_px"))
  
  expect_s4_class(x, "SingleCellExperiment")
  expect_true(all(c("CenterX_global_px", "CenterY_global_px") %in% colnames(SingleCellExperiment::colData(x))))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(8, 9)))
  expect_s4_class(SingleCellExperiment::counts(x), "dgCMatrix")
})
