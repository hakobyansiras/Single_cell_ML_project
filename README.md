# Single_cell_ML_project
This repository stores scripts for single cell analysis project for Intro to ML class of Computational Bilogy Master's degree program at Russian-Armenian University.

## Data:
Zenodo URL: https://zenodo.org/record/4569496#.Y8emVHb7SF4
Samples: Single cell RNA-Seq adult skin samples from healthy, eczema, and psoriasis groups

## Workflow:

### 1. Seurat Analysis
Complete Seurat workflow for preprocessing and clustering the data. <br>
Script:: seurat_analysis.R (part 1) <br>
Input:: Data in AnnData (h5ad) format <br>

#### Steps
* Normalization and Scaling
* Principal Component Analysis (PCA)
* Clustering 
* Visualizations (UMAP)


### 2. Metacell Generation
Downsampling using subdividing and  refinement k-means clustering. <br>
Script:: seurat_analysis.R (part 2) <br>
Input:: Data in AnnData (h5ad) format <br>

#### Steps
* unsupervised cell clustering as provided by Seurat
* each cluster divided into patient and cell type specific sub-clusters (cluster X cell_type X patient_id)
* each subgroup is subclustered by a refinement k-means clustering (using centroids as average of each subcluster)


### 3. SOM Analysis
Self-organizing maps for feature reduction (compression into another featureset) and downstream analysis <br>
Script:: som_analysis.R  <br>
Input:: metecell matrix (5970 x 20k)  <br>

#### Steps
* low-dimensional representation of a higher dimensional data 
* condition specific gene-expression portraits
* definition of metagenes (groups of co-differentially expressed genes)
* visualizations (heatmaps, correlation networks, SOM maps)

### 4. Classification
ML Classification algorithms for condition and cell-type classification  <br>
Script(s):: decision_tree.py, random_forest.py, svm.py, adaboost.py, xgboostclassifier.py  <br>
Input:: 3 different datasets  <br>
  * Metacell, metagene data - reduced sample and feature size (5970x3600)
  * Metacell, variable features - reduced sample size, filtered features (5970x2000)
  * Random cells, variable features - filtered samples and features (5970x2000)


#### Steps
* training models
* hyperparameter tuning (Grid Search)
* measuring performance on test set (accuracy, balanced accuarcy f1 score, per class f1 score)







