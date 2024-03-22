dir <- system.file(
  file.path("extdata", "Xenium_small"), 
  package="SpatialExperimentIO")

test_that("example data folders uniquely contains needed files", {
  expect_true("cell_feature_matrix.h5" %in% list.files(dir))
  expect_true("cells.csv.gz" %in% list.files(dir)) 
  
  expect_length(list.files(dir, "cell_feature_matrix.h5"), 1)
  expect_length(list.files(dir, "cells.csv.gz"), 1)
})

test_that("needed files contains spatial columns of interest", {
  metadata <- read.csv(gzfile(file.path(dir, "cells.csv.gz")), header = TRUE)
  
  expect_true(all(c("x_centroid", "y_centroid") %in% colnames(metadata))) 
  expect_true(is.numeric(metadata$x_centroid))
  expect_true(is.numeric(metadata$y_centroid))
})

test_that("data are read correctly to SpatialExperiment class", {
  x <- readXeniumSXE(dirname = dir, 
                     return_type = "SPE",
                     countmatfpattern = "cell_feature_matrix.h5",
                     metadatafpattern = "cells.csv.gz", 
                     coord_names = c("x_centroid", "y_centroid"))
  
  expect_s4_class(x, "SpatialExperiment")
  expect_true(all(colnames(SpatialExperiment::spatialCoords(x)) == c("x_centroid", "y_centroid")))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(4, 6)))
  expect_s4_class(SingleCellExperiment::counts(x), "dgCMatrix")
})

test_that("data are read correctly to SingleCellExperiment class", {
  x <- readXeniumSXE(dirname = dir, 
                     return_type = "SCE",
                     countmatfpattern = "cell_feature_matrix.h5",
                     metadatafpattern = "cells.csv.gz", 
                     coord_names = c("x_centroid", "y_centroid"))
  
  expect_s4_class(x, "SingleCellExperiment")
  expect_true(all(c("x_centroid", "y_centroid") %in% colnames(SingleCellExperiment::colData(x))))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(4, 6)))
  expect_s4_class(SingleCellExperiment::counts(x), "dgCMatrix")
})
