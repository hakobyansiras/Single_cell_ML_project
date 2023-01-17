import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import xgboost as xgb
from sklearn.metrics import accuracy_score, f1_score, balanced_accuracy_score, classification_report
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV


## reading data
variable_feature_data = pd.read_csv('variable_features_subset_cell_type_labeled.tsv', delimiter='\t')
# variable_feature_data["binary_disease_state"] = np.where(variable_feature_data['disease_state']!= 0, 1, 0) 

## encoding character labels into numeric
le = preprocessing.LabelEncoder()
le.fit(variable_feature_data.cell_type)
variable_feature_data['categorical_label'] = le.transform(variable_feature_data.cell_type)

## target names for per class accuracy
target_names =  ['DC1','DC2','Differentiated_KC','Differentiated_KC*','F1','F2','F3'
,'ILC1_3','ILC1_NK','ILC2','Inf_mac','LC_1','LC_2','LC_3','LC_4','LE1'
,'LE2','Macro_1','Macro_2','Mast_cell','Melanocyte','MigDC','Mono_mac'
,'NK','Pericyte_1','Pericyte_2','Plasma','Proliferating_KC','Schwann_1'
,'Schwann_2','Tc','Tc17_Th17','Tc_IL13_IL22','Th','Treg'
,'Undifferentiated_KC','VE1','VE2','VE3','moDC_1','moDC_2','moDC_3']

## split
X_train, X_test, y_train, y_test = train_test_split(variable_feature_data.drop(columns=['cell_type', 'categorical_label']),
                                                    variable_feature_data.categorical_label,
                                                    test_size=0.2,
                                                    random_state=0)
## declare model and train
xgc = xgb.XGBClassifier()
xgc.fit(X_train, y_train)

#predict and evaluate 
preds = xgc.predict(X_test)
print("Test Accuracy:", accuracy_score(y_test, preds))
print("Test F1 score:", f1_score(y_test, preds, average='weighted'))
print("Test Balanced Accuracy:", balanced_accuracy_score(y_test, preds))
print("Test Classification Report:", classification_report(y_test, preds, target_names=target_names, digits=4))