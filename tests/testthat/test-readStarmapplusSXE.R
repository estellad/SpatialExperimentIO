dir <- system.file(
  file.path("extdata", "STARmapPLUS_small"), 
  package="SpatialExperimentIO")

test_that("example data folders uniquely contains needed files", {
  expect_true("mockraw_expression_pd.csv" %in% list.files(dir))
  expect_true("mock_spatial.csv" %in% list.files(dir)) 
  
  expect_length(list.files(dir, "spatial.csv"), 1)
  expect_length(list.files(dir, "raw_expression_pd.csv"), 1)
})

test_that("needed files contains spatial columns of interest", {
  metadata <- read.csv(file.path(dir, "mock_spatial.csv"))
  
  expect_true(all(c("X", "Y", "Z") %in% colnames(metadata))) 
})

test_that("data are read correctly to SpatialExperiment class", {
  x <- readStarmapplusSXE(dirname = dir, 
                          return_type = "SPE",
                          countmatfpattern = "raw_expression_pd.csv", 
                          metadatafpattern = "spatial.csv", 
                          coord_names = c("X", "Y", "Z"))
  
  expect_s4_class(x, "SpatialExperiment")
  expect_true(all(colnames(SpatialExperiment::spatialCoords(x)) == c("X", "Y", "Z")))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(8, 9)))
  expect_s4_class(SingleCellExperiment::counts(x), "dgCMatrix")
})

test_that("data are read correctly to SingleCellExperiment class", {
  x <- readStarmapplusSXE(dirname = dir, 
                          return_type = "SCE",
                          countmatfpattern = "raw_expression_pd.csv", 
                          metadatafpattern = "spatial.csv", 
                          coord_names = c("X", "Y", "Z"))
  
  expect_s4_class(x, "SingleCellExperiment")
  expect_true(all(c("X", "Y", "Z") %in% colnames(SingleCellExperiment::colData(x))))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(8, 9)))
  expect_s4_class(SingleCellExperiment::counts(x), "dgCMatrix")
})
