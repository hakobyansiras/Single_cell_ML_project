library(caTools)
library(e1071)
library(MLmetrics)

svm_function <- function(data_matrix, labels, scale_data = FALSE) {
  
  data_matrix <- rbind(data_matrix, labels)
  
  data_matrix <- t(data_matrix)
  data_matrix <- as.data.frame(data_matrix)
  
  ### splitting dataset
  set.seed(123)
  split = sample.split(data_matrix$labels, SplitRatio = 0.75)
  training_set = subset(data_matrix, split == TRUE)
  test_set = subset(data_matrix, split == FALSE)
  
  ### need to discuss if we need scaling
  if(scale_data) {
    training_set[-ncol(data_matrix)] = scale(training_set[-ncol(data_matrix)])
    test_set[-ncol(data_matrix)] = scale(test_set[-ncol(data_matrix)])
  }
  
  ### running svm on train dataset
  classifierR = svm(formula = labels ~ .,
                    data = training_set,
                    type = 'C-classification', # this is because we want to make a regression classification
                    kernel = 'radial',
                    degree = 3
  )
  
  
  y_predR = predict(classifierR, newdata = test_set[-ncol(data_matrix)])
  
  cmR <- table(test_set[, ncol(data_matrix)], y_predR)
  
  f1_score <- MLmetrics::F1_Score(test_set[, ncol(data_matrix)], y_predR)
 
  return(list(cmR = cmR, f1_score = f1_score))
   
}

### loading SOM metagene data
load("skin_metacell_som/skin_metacells - Results/skin_metacells pre.RData")

svm_metadata <- env$metadata

# converting labels
sapply(metacell_list[[1]]$metacellLabels, function(x) {unlist(strsplit(x, split = " "))[1]})

metacell_labels <- sapply(metacell_list[[1]]$metacellLabels, function(x) {unlist(strsplit(x, split = " "))[1]})

metacell_lables_to_cell_id <- sapply(sapply(colnames(svm_metadata), function(x) {paste(unlist(strsplit(x, split = "_"))[1:3], collapse = "_")}), function(x) {
  
  names(metacell_labels)[which(metacell_labels %in% x)]
  
})

## creating disaese state labels
disease_state <- unname(setNames(object = c(1, 2, 0), nm = c("Eczema", "Psoriasis", "Healthy"))[as.character(skin_metadata[sapply(metacell_lables_to_cell_id, function(x) {x[1]}), "Status"])])

## svm on metagene metacell data
metagene_metacell_svm <- svm_function(data_matrix = svm_metadata, labels = disease_state, scale_data = TRUE)


cell_type_in_metacells <- as.character(skin_metadata[sapply(metacell_lables_to_cell_id, function(x) {x[1]}), "full_clustering"])
## filter by these cell types
cell_type_svm_f1_score <- lapply(c("F1", "F2", "F3", "Proliferating_KC",  "Differentiated_KC", "Pericyte"), function(x) {
  
  svm_function(data_matrix = svm_metadata[,grep(x, cell_type_in_metacells)], labels = disease_state[grep(x, cell_type_in_metacells)], scale_data = TRUE)$f1_score
  
})

names(cell_type_svm_f1_score) <- c("F1", "F2", "F3", "Proliferating_KC",  "Differentiated_KC", "Pericyte")


#### svm on variable features without compression
variable_features <- skin_data_seurat@assays$RNA@scale.data

variable_features <- rbind(variable_features, unname(setNames(object = c(1, 2, 0), nm = c("Eczema", "Psoriasis", "Healthy"))[as.character(skin_metadata[colnames(variable_features), "Status"])]))

variable_features_subset <- variable_features[,sapply(metacell_lables_to_cell_id, function(x) {x[1]})]

variable_features_subset <- t(variable_features_subset)
variable_features_subset <- as.data.frame(variable_features_subset)
colnames(variable_features_subset)[ncol(variable_features_subset)] <- "disease_state"

write.table(variable_features_subset, file = "variable_features_subset.tsv", sep = "\t", quote = F, row.names = F)

library(caTools)
set.seed(123)
split = sample.split(variable_features_subset$disease_state, SplitRatio = 0.75)
training_set = subset(variable_features_subset, split == TRUE)
test_set = subset(variable_features_subset, split == FALSE)

### need to discuss if we need scaling
# training_set[-2001] = scale(training_set[-2001])
# test_set[-2001] = scale(test_set[-2001])

### running svm on train dataset
library(e1071)
classifierR = svm(formula = disease_state ~ .,
                  data = training_set,
                  type = 'C-classification', # this is because we want to make a regression classification
                  kernel = 'radial',
                  degree = 3
)


y_predR = predict(classifierR, newdata = test_set[-2001])

cmR = table(test_set[, 2001], y_predR)
cmR

library("MLmetrics")
MLmetrics::F1_Score(test_set[, 2001], y_predR)
