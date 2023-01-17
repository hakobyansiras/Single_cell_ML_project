# loading required packages
library(caTools)
library(e1071)
library(MLmetrics)

svm_function <- function(data_matrix, scale_data = FALSE) {
  
  print("Processing data")
  ### splitting dataset
  set.seed(123)
  split = sample.split(data_matrix$labels, SplitRatio = 0.75)
  training_set = subset(data_matrix, split == TRUE)
  test_set = subset(data_matrix, split == FALSE)
  
  ### scaling data
  if(scale_data) {
    training_set[-ncol(data_matrix)] = scale(training_set[-ncol(data_matrix)])
    test_set[-ncol(data_matrix)] = scale(test_set[-ncol(data_matrix)])
  }
  
  print("Running SVM")
  ### running svm on train dataset
  classifierR = svm(formula = labels ~ .,
                    data = training_set,
                    type = 'C-classification', # this is because we want to make a regression classification
                    kernel = 'radial',
                    degree = 3
  )
  
  print("Running prediction")
  ### running prediction 
  y_predR = predict(classifierR, newdata = test_set[-ncol(data_matrix)])
  
  cmR <- table(test_set[, ncol(data_matrix)], y_predR)
  
  f1_score <- MLmetrics::F1_Score(test_set[, ncol(data_matrix)], y_predR)
 
  return(list(cmR = cmR, f1_score = f1_score, classifierR = classifierR))
   
}

### loading SOM metagene data
load("skin_metacell_som/skin_metacell_res/skin_metacells.RData")

svm_metadata <- env$metadata

### converting labels
sapply(metacell_list[[1]]$metacellLabels, function(x) {unlist(strsplit(x, split = " "))[1]})

metacell_labels <- sapply(metacell_list[[1]]$metacellLabels, function(x) {unlist(strsplit(x, split = " "))[1]})

metacell_lables_to_cell_id <- sapply(sapply(colnames(svm_metadata), function(x) {paste(unlist(strsplit(x, split = "_"))[1:3], collapse = "_")}), function(x) {
  
  names(metacell_labels)[which(metacell_labels %in% x)]
  
})

## creating disaese state labels
disease_state <- unname(setNames(object = c(1, 1, 0), nm = c("Eczema", "Psoriasis", "Healthy"))[as.character(skin_metadata[sapply(metacell_lables_to_cell_id, function(x) {x[1]}), "Status"])])

#### svm on metagene metacell data ####
svm_metadata_labeled <- rbind(svm_metadata, disease_state)

svm_metadata_labeled <- t(svm_metadata_labeled)
svm_metadata_labeled <- as.data.frame(svm_metadata_labeled)
colnames(svm_metadata_labeled)[ncol(svm_metadata_labeled)] <- "labels"

metagene_metacell_svm <- svm_function(data_matrix = svm_metadata_labeled, scale_data = TRUE)


#### svm on metagene metacell data for each cell type ####
cell_type_in_metacells <- as.character(skin_metadata[sapply(metacell_lables_to_cell_id, function(x) {x[1]}), "full_clustering"])
## filter by these cell types
cell_type_svm_f1_score <- lapply(c("F1", "F2", "F3", "Proliferating_KC",  "Differentiated_KC", "Pericyte"), function(x) {
  
  svm_function(data_matrix = svm_metadata_labeled[grep(x, cell_type_in_metacells),], scale_data = TRUE)$f1_score
  
})

names(cell_type_svm_f1_score) <- c("F1", "F2", "F3", "Proliferating_KC",  "Differentiated_KC", "Pericyte")

#### svm on metacell variable feature data ####
metacell_data_var_features <- metacell_list[[2]][which(rownames(metacell_list[[2]]) %in% VariableFeatures(skin_data_seurat)),]
metacell_data_var_features <- rbind(metacell_data_var_features, disease_state)
metacell_data_var_features <- t(metacell_data_var_features)
metacell_data_var_features <- as.data.frame(metacell_data_var_features)
colnames(metacell_data_var_features)[ncol(metacell_data_var_features)] <- "labels"

### Exporting data for external analysis
write.table(metacell_data_var_features, file = "metacell_data_var_features.tsv", sep = "\t", quote = F, row.names = F) 
### running svm
metacell_data_var_features_svm <- svm_function(data_matrix = metacell_data_var_features, scale_data = FALSE)


#### svm on variable features without compression ####
variable_features <- skin_data_seurat@assays$RNA@scale.data

variable_features <- rbind(variable_features, unname(setNames(object = c(1, 1, 0), nm = c("Eczema", "Psoriasis", "Healthy"))[as.character(skin_metadata[colnames(variable_features), "Status"])]))

variable_features_subset <- t(variable_features)
variable_features_subset <- as.data.frame(variable_features_subset)

colnames(variable_features_subset)[ncol(variable_features_subset)] <- "labels"

random_cell_matrixes <- lapply(1:10, function(y) {variable_features_subset[sapply(metacell_lables_to_cell_id, function(x) {sample(x, size = 1)}),]})

random_var_features_svm <- lapply(random_cell_matrixes, function(x) {
  
  svm_function(data_matrix = x, scale_data = FALSE)
  
})


#### svm on metacell PSF results ####

metacell.data <- metacell_list[[2]]

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
symbol_to_entrez <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                           mart = ensembl)
symbol_to_entrez <- symbol_to_entrez[which(!is.na(symbol_to_entrez$entrezgene_id)), ]


symbol_to_entrez_filtered <- symbol_to_entrez[which(symbol_to_entrez$hgnc_symbol %in% rownames(metacell.data)),]

symbol_to_entrez_filtered <- symbol_to_entrez_filtered[which(!duplicated(symbol_to_entrez_filtered$entrezgene_id)),]

symbol_to_entrez_filtered <- setNames(object = symbol_to_entrez_filtered$entrezgene_id, nm = symbol_to_entrez_filtered$hgnc_symbol)

metacell.data <- metacell.data[which(rownames(metacell.data) %in% names(symbol_to_entrez_filtered)),]

rownames(metacell.data) <- symbol_to_entrez_filtered[rownames(metacell.data)]

metacell.data <- metacell.data + 1

metacell.data <- metacell.data/rowMeans(metacell.data)

load(system.file("extdata", "edited_pathways_new.RData", package="psf"))

metacell_psf <- psf::run_psf(entrez.fc = metacell.data, kegg.collection = edited_pathways_new, calculate.significance = F, ncores = 30)

save(metacell_psf, file = "metacell_psf.RData")

psf_matrix <- simplify2array(Reduce(rbind, lapply(metacell_psf[names(unlist(per_pathway_svm)[which(unlist(per_pathway_svm) > 0.2)])], function(x) {
  
  x$psf_activities[x$sink.nodes,]
  
})))



disease_state <- unname(setNames(object = c(1, 1, 0), nm = c("Eczema", "Psoriasis", "Healthy"))[as.character(skin_metadata[sapply(metacell_lables_to_cell_id, function(x) {x[1]}), "Status"])])

metacell.data_filtered <- metacell.data[which(rownames(metacell.data) %in% unique(unlist(lapply(edited_pathways_new, function(x) {unique(unname(unlist(graph::nodeData(x$graph, attr = "genes"))))})))),]

psf_matrix_with_labels <- rbind(metacell.data_filtered, disease_state)

psf_matrix_with_labels <- t(psf_matrix_with_labels)
psf_matrix_with_labels <- as.data.frame(psf_matrix_with_labels)


### splitting dataset
library(caTools)
set.seed(123)
split = sample.split(psf_matrix_with_labels$disease_state, SplitRatio = 0.75)
training_set = subset(psf_matrix_with_labels, split == TRUE)
test_set = subset(psf_matrix_with_labels, split == FALSE)

### need to discuss if we need scaling
training_set[-ncol(psf_matrix_with_labels)] = scale(training_set[-ncol(psf_matrix_with_labels)])
test_set[-ncol(psf_matrix_with_labels)] = scale(test_set[-ncol(psf_matrix_with_labels)])

### running svm on train dataset
library(e1071)
classifierR = svm(formula = disease_state ~ .,
                  data = training_set,
                  type = 'C-classification', # this is because we want to make a regression classification
                  kernel = 'radial'
)


y_predR = predict(classifierR, newdata = test_set[-ncol(psf_matrix_with_labels)])

cmR = table(test_set[, ncol(psf_matrix_with_labels)], y_predR)
cmR

library("MLmetrics")
MLmetrics::F1_Score(test_set[, ncol(psf_matrix_with_labels)], y_predR)

### running svm for each pathway
per_pathway_svm <- lapply(metacell_psf, function(x) {
  
  print(x$attrs$title)
  
  psf_matrix <- x$psf_activities[x$sink.nodes,]
  
  psf_matrix_with_labels <- rbind(psf_matrix, disease_state)
  
  psf_matrix_with_labels <- t(psf_matrix_with_labels)
  psf_matrix_with_labels <- as.data.frame(psf_matrix_with_labels)
  
  ### splitting dataset
  set.seed(123)
  split = sample.split(psf_matrix_with_labels$disease_state, SplitRatio = 0.75)
  training_set = subset(psf_matrix_with_labels, split == TRUE)
  test_set = subset(psf_matrix_with_labels, split == FALSE)
  
  ### need to discuss if we need scaling
  training_set[-ncol(psf_matrix_with_labels)] = scale(training_set[-ncol(psf_matrix_with_labels)])
  test_set[-ncol(psf_matrix_with_labels)] = scale(test_set[-ncol(psf_matrix_with_labels)])
  
  ### running svm on train dataset
  classifierR = svm(formula = disease_state ~ .,
                    data = training_set,
                    type = 'C-classification', # this is because we want to make a regression classification
                    kernel = 'radial',
                    degree = 3
  )
  
  
  y_predR = predict(classifierR, newdata = test_set[-ncol(psf_matrix_with_labels)])
  
  MLmetrics::F1_Score(test_set[, ncol(psf_matrix_with_labels)], y_predR)
  
})

names(per_pathway_svm) <- names(metacell_psf)
