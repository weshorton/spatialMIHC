# spatialMIHC

Package to create a Seurat spatial object from mIHC outputs and run various analyses, generate plots, etc.

# Overview

## Required Inputs

In order to build the Seurat object, the following are required:

1. Classified_CSV - there is one `.csv` file for each ROI in a given study. These are the output of <insert processing tool> and contain a cell classification for each cell based on the hierarchical gating structure
1. Metadata - there should be a metadata file with one entry for each ROI. The Sample_ID column should correspond to the file name in Classified_CSV
1. Color file - a file that maps cell classifications to specific colors is also required.

## Object Creation

Once the required inputs are gathered, the seurat object can be made.

1. Remove Duplicates - `removeAnnotationDuplicates()`
  - Due to the nature of the gating strategy, some cells can be classified as two different populations and show up twice
  - As of now, we are only taking the first (most general) classification for each cell

1. Reformat For Seurat - `reformatForSeurat()`
  - Place the X/Y coordinates from Classified_CSV into appropriate slot
  - Place marker intensity from Classified_CSV into appropriate slot
  - Add cell-level metadata (i.e. class) to object
  - `idCol_v`: default of `ObjectNumber` should work
  - `coordCols_v`: default of `c("OBJECTID", "Location_Center_X", "Location_Center_Y")` should work
  - `exprCols_v`: anything that's not `idCol_v`, `coordCols_v`, `_func` columns, etc. Just the marker intensity for each population class
  - `metaCols_v`: cell-level attributes, such as `class` and any `_func` markers.
  
1. Build Object - `buildSeurat()`
  - Create a Seurat S4 object using the provided input data from reformat
  - The slideName_v argument is the most granular ID (i.e. ROI's ID). This should correspond to the file name of classified CSV.

## Analysis

This object can then be used to run various spatial analyses in R.  
Other cell-level annotations from different tools (like Katie's) can be loaded and added to the metadata as well.  