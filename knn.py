import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.neighbors import KNeighborsClassifier
from sklearn import preprocessing
from sklearn.metrics import accuracy_score, f1_score, balanced_accuracy_score, classification_report
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

#Load data
variable_feature_random_cell = pd.read_csv('variable_features_subset.tsv', delimiter='\t')
metadata_metacell = pd.read_csv('svm_metadata_with_labels.tsv', delimiter='\t')
variable_feature_metacell = pd.read_csv('metacell_data_var_features.tsv', delimiter='\t')

#Binary class
variable_feature_random_cell['disease_state'] = np.where(variable_feature_random_cell['disease_state'] == 0, 0, 1)
metadata_metacell['disease_state'] =  np.where(metadata_metacell['disease_state'] == 0, 0, 1)
variable_feature_metacell['disease_state'] =  np.where(variable_feature_metacell['disease_state'] == 0, 0, 1)

def run_knn(data, title=''):
    print('Running knn binary classification on ', title)
    #split
    X_train, X_test, y_train, y_test = train_test_split(data.drop(columns='disease_state'),
                                                        data.disease_state,
                                                        test_size=0.2,
                                                        random_state=0)

    #parameter tuning
    param_grid = {'n_neighbors': [3, 5, 7, 9],
              'weights': ['uniform', 'distance'],
              'metric': ['euclidean', 'manhattan']}

    dt = KNeighborsClassifier()
    clf = GridSearchCV(dt, param_grid, cv=5)
    clf.fit(X_train, y_train)

    print("Best Estimator:", clf.best_estimator_)
    print("Train Accuracy:", clf.best_score_)

    #predict
    preds = clf.best_estimator_.predict(X_test)
    print("Test Accuracy:", accuracy_score(y_test, preds))
    print("Test F1 score:", f1_score(y_test, preds, average='weighted'))
    print('\n\n')


run_knn(data=variable_feature_random_cell, title='variable_feature - random_cell')
run_knn(data=metadata_metacell, title='metadata - metacell')
run_knn(data=variable_feature_metacell, title='variable_feature - metacell')
