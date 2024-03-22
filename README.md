# Introduction

The `SpatialExperimentIO` package provides a set of functions to import Xenium (10x Genomics), CosMx (Nanostring), MERSCOPE (Vizgen), STARmapPLUS (Wang et al., 2023, Broad Institute) data into a `SpatialExperiment` or `SingleCellExperiment`class object.

# Installation

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SpatialExperimentIO")
```

# Development version

You can also install the development version of *SpatialExperimentIO* from [GitHub](https://github.com/estellad/SpatialExperimentIO) with:

``` r
# install.packages("devtools")
devtools::install_github("estellad/SpatialExperimentIO")
```

# Load package

``` r
library(SpatialExperimentIO)
```

# Xenium ouptut folder structure

A standard Xenium output folder should contain these files for the function `readXeniumSXE()`. Either `cell_feature_matrix.h5` or `/cell_feature_matrix`, and cells.csv.gz are required.

```         
    Xenium_unzipped
        └── outs 
            ├── cells.csv.gz 
            ├── cell_feature_matrix.h5 
            └── cell_feature_matrix 
                ├── barcodes.tsv 
                ├── features.tsv 
                └── matrix.mtx
    
```

# CosMx output folder structure

A standard CosMx output folder should contain these files for the function `readCosMxSXE()`. Both `exprMat_file.csv` and `metadata_file.csv` is required.

```         
    CosMx 
        ├── metadata_file.csv 
        └── exprMat_file.csv
```

# MERSCOPE output folder structure

A standard MERSCOPE output folder should contain these files for the function `readMerscopeSXE()`. Both `cell_by_gene.csv` and `cell_metadata.csv` is required.

```         
    MERSCOPE 
        ├── cell_metadata.csv 
        └── cell_by_gene.csv
```

# STARmap PLUS output folder structure

A standard STARmap PLUS output folder should contain these files for the function `readStarmapplusSXE()`. Both `raw_expression_pd.csv` and `spatial.csv` is required.

```         
    STARmap_PLUS 
        ├── spatial.csv 
        └── raw_expression_pd.csv
```

# Usage

Taking Xenium as an example, providing a path to the folder that stores all the required files (i.e. `/outs` ) would return a `SpatialExperiment` object.

``` r
spe <- readXeniumSXE(dir)
spe
# class: SpatialExperiment 
# dim: 4 6 
# metadata(0):
# assays(1): counts
# rownames(4): AATK ABL1 ACKR3 ACKR4
# rowData names(3): ID Symbol Type
# colnames(6): 1 2 ... 5 6
# colData names(9): X cell_id ... nucleus_area sample_id
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(0):
```
