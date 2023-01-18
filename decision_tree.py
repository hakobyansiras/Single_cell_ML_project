import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.tree import DecisionTreeClassifier
from sklearn import preprocessing
from sklearn.metrics import accuracy_score, f1_score, balanced_accuracy_score, classification_report
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

#MULTI CLASS
#Load data
variable_feature_random_cell = pd.read_csv('variable_features_subset.tsv', delimiter='\t')
metadata_metacell = pd.read_csv('svm_metadata_with_labels.tsv', delimiter='\t')
variable_feature_metacell = pd.read_csv('metacell_data_var_features.tsv', delimiter='\t')

def run_decision_tree(data, title=''):
    print('Running Decision tree on ', title)
    #split
    X_train, X_test, y_train, y_test = train_test_split(data.drop(columns='disease_state'),
                                                        data.disease_state,
                                                        test_size=0.2,
                                                        random_state=0)

    #parameter tuning
    parameters = {'criterion':('gini', 'entropy'),
                'max_depth': [3, 5, 10, 25],
                'min_samples_split': [2, 4, 6, 10],
                'min_impurity_decrease': [0, 0.1, 0.15]}

    dt = DecisionTreeClassifier()
    clf = GridSearchCV(dt, parameters)
    clf.fit(X_train, y_train)

    print("Best Estimator:", clf.best_estimator_)
    print("Train Accuracy:", clf.best_score_)

    #predict
    preds = clf.best_estimator_.predict(X_test)
    print("Test Accuracy:", accuracy_score(y_test, preds))
    print("Test F1 score:", f1_score(y_test, preds, average='weighted'))
    print('\n\n')


run_decision_tree(data=variable_feature_random_cell, title='variable_feature - random_cell')
run_decision_tree(data=metadata_metacell, title='metadata - metacell')
run_decision_tree(data=variable_feature_metacell, title='variable_feature - metacell')


BINARY CLASS
variable_feature_random_cell['disease_state'] = np.where(variable_feature_random_cell['disease_state'] == 0, 0, 1)
metadata_metacell['disease_state'] =  np.where(metadata_metacell['disease_state'] == 0, 0, 1)
variable_feature_metacell['disease_state'] =  np.where(variable_feature_metacell['disease_state'] == 0, 0, 1)

def run_decision_tree_binary(data, title=''):
    print('Running Decision tree binary classification on ', title)
    #split
    X_train, X_test, y_train, y_test = train_test_split(data.drop(columns='disease_state'),
                                                        data.disease_state,
                                                        test_size=0.2,
                                                        random_state=0)

    #parameter tuning
    parameters = {'criterion':('gini', 'entropy'),
                'max_depth': [3, 5, 10, 25],
                'min_samples_split': [2, 4, 6, 10],
                'min_impurity_decrease': [0, 0.1, 0.15]}

    dt = DecisionTreeClassifier()
    clf = GridSearchCV(dt, parameters)
    clf.fit(X_train, y_train)

    print("Best Estimator:", clf.best_estimator_)
    print("Train Accuracy:", clf.best_score_)

    #predict
    preds = clf.best_estimator_.predict(X_test)
    print("Test Accuracy:", accuracy_score(y_test, preds))
    print("Test F1 score:", f1_score(y_test, preds, average='weighted'))
    print('\n\n')


run_decision_tree_binary(data=variable_feature_random_cell, title='variable_feature - random_cell')
run_decision_tree_binary(data=metadata_metacell, title='metadata - metacell')
run_decision_tree_binary(data=variable_feature_metacell, title='variable_feature - metacell')


#CELL TYPE CLASSIFICATION
variable_feature_random_cell_type = pd.read_csv('variable_features_subset_cell_type_labeled.tsv', delimiter='\t')
metacell_metagene_cell_type = pd.read_csv('metacell_metagene_cell_type_labeled.tsv', delimiter='\t')


le = preprocessing.LabelEncoder()
le.fit(variable_feature_random_cell_type.cell_type)
variable_feature_random_cell_type['categorical_label'] = le.transform(variable_feature_random_cell_type.cell_type)

le = preprocessing.LabelEncoder()
le.fit(metacell_metagene_cell_type.cell_type)
metacell_metagene_cell_type['categorical_label'] = le.transform(metacell_metagene_cell_type.cell_type)

target_names =  ['DC1','DC2','Differentiated_KC','Differentiated_KC*','F1','F2','F3'
,'ILC1_3','ILC1_NK','ILC2','Inf_mac','LC_1','LC_2','LC_3','LC_4','LE1'
,'LE2','Macro_1','Macro_2','Mast_cell','Melanocyte','MigDC','Mono_mac'
,'NK','Pericyte_1','Pericyte_2','Plasma','Proliferating_KC','Schwann_1'
,'Schwann_2','Tc','Tc17_Th17','Tc_IL13_IL22','Th','Treg'
,'Undifferentiated_KC','VE1','VE2','VE3','moDC_1','moDC_2','moDC_3']

def run_decision_tree_cell_type(data, title=''):
    print('Running Decision tree cell type classification on ', title)
    #split
    X_train, X_test, y_train, y_test = train_test_split(data.drop(columns=['cell_type', 'categorical_label']),
                                                        data.categorical_label,
                                                        test_size=0.2,
                                                        random_state=0)

    #parameter tuning
    parameters = {'criterion':('gini', 'entropy'),
                'max_depth': [3, 5, 10, 25],
                'min_samples_split': [2, 4, 6, 10],
                'min_impurity_decrease': [0, 0.1, 0.15]}

    dt = DecisionTreeClassifier()
    clf = GridSearchCV(dt, parameters)
    clf.fit(X_train, y_train)

    print("Best Estimator:", clf.best_estimator_)
    print("Train Accuracy:", clf.best_score_)

    #predict
    preds = clf.best_estimator_.predict(X_test)
    print("Test Accuracy:", accuracy_score(y_test, preds))
    print("Test F1 score:", f1_score(y_test, preds, average='weighted')) #, labels=np.unique(preds) doesn't use labels that were not predicted
    print("Test Balanced Accuracy:", balanced_accuracy_score(y_test, preds))
    print("Test Classification Report:", classification_report(y_test, preds, target_names=target_names, digits=4))
    print('\n\n')

run_decision_tree_cell_type(data=variable_feature_random_cell_type, title='variable_feature - random_cell')
run_decision_tree_cell_type(data=metacell_metagene_cell_type, title='metagene - metacell')
