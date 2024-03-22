dir <- system.file(
  file.path("extdata", "MERSCOPE_small"), 
  package="SpatialExperimentIO")

test_that("example data folders uniquely contains needed files", {
  expect_true("ovarian_p2s2_cell_by_gene.csv" %in% list.files(dir))
  expect_true("ovarian_p2s2_cell_metadata.csv" %in% list.files(dir)) 
  
  expect_length(list.files(dir, "cell_by_gene.csv"), 1)
  expect_length(list.files(dir, "cell_metadata.csv"), 1)
})

test_that("needed files contains spatial columns of interest", {
  colData <- read.csv(file.path(dir, "ovarian_p2s2_cell_metadata.csv"))
  
  expect_true(all(c("center_x", "center_y") %in% colnames(colData))) 
  expect_true(is.numeric(colData$center_x))
  expect_true(is.numeric(colData$center_y))
})

test_that("data are read correctly to SpatialExperiment class", {
  x <- readMerscopeSXE(dirname = dir, 
                       return_type = "SPE",
                       countmatfpattern = "cell_by_gene.csv", 
                       metadatafpattern = "cell_metadata.csv", 
                       coord_names = c("center_x", "center_y"))
  
  expect_s4_class(x, "SpatialExperiment")
  expect_true(all(colnames(SpatialExperiment::spatialCoords(x)) == c("center_x", "center_y")))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(9, 8)))
  expect_s4_class(SpatialExperiment::counts(x), "dgCMatrix")
})

test_that("data are read correctly to SingleCellExperiment class", {
  x <- readMerscopeSXE(dirname = dir, 
                       return_type = "SCE",
                       countmatfpattern = "cell_by_gene.csv", 
                       metadatafpattern = "cell_metadata.csv", 
                       coord_names = c("center_x", "center_y"))
  
  expect_s4_class(x, "SingleCellExperiment")
  expect_true(all(c("center_x", "center_y") %in% colnames(SpatialExperiment::colData(x))))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(9, 8)))
  expect_s4_class(SpatialExperiment::counts(x), "dgCMatrix")
})
