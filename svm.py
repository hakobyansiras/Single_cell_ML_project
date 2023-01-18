import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV

from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

#Load data
variable_feature_random_cell = pd.read_csv('metacell_metagene_cell_type_labeled.tsv', delimiter='\t')
metadata_metacell = pd.read_csv('variable_features_subset_cell_type_labeled.tsv', delimiter='\t')

 
#BINARY CLASS
metacell_metagene_cell_type_labeled['disease_state'] = np.where(metacell_metagene_cell_type_labeled['disease_state'] == 0, 0, 1)
variable_features_subset_cell_type_labeled['disease_state'] =  np.where(variable_features_subset_cell_type_labeled['disease_state'] == 0, 0, 1)

def run_svm(data, title=''):
    print('Running svm on ', title)
    #split
    X_train, X_test, y_train, y_test = train_test_split(data.drop(columns='disease_state'),
                                                        data.disease_state,
                                                        test_size=0.2,
                                                        random_state=0)

    #parameter tuning
    parameters = {'C': [0.1, 1, 10],
              'kernel': ['linear', 'poly', 'rbf', 'sigmoid'],
              'degree': [2, 3, 4],
              'gamma': ['scale', 'auto']}

    svm = SVC()
    clf = GridSearchCV(svm, parameters, cv=5)
    clf.fit(X_train, y_train)

    print("Best Estimator:", clf.best_estimator_)
    print("Train Accuracy:", clf.best_score_)

    #predict
    preds = clf.best_estimator_.predict(X_test)
    print("Test Accuracy:", accuracy_score(y_test, preds))
    print("Test F1 score:", f1_score(y_test, preds, average='weighted'))
    print('\n\n')


run_svm(data=metacell_metagene_cell_type_labeled, title='metacell_metagene_cell_type_labeled')
run_svm(data=variable_features_subset_cell_type_labeled, title='variable_features_subset_cell_type_labeled')
