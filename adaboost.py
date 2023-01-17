import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import accuracy_score, f1_score, balanced_accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

#Load data
#variable feature-random cell
variable_feature_data = pd.read_csv('svm_metadata_with_labels.tsv', delimiter='\t')
variable_feature_data["binary_disease_state"] = np.where(variable_feature_data['disease_state']!= 0, 1, 0) 

#split
X_train, X_test, y_train, y_test = train_test_split(variable_feature_data.drop(columns=['disease_state', 'binary_disease_state']),
                                                    variable_feature_data.binary_disease_state,
                                                    test_size=0.2,
                                                    random_state=0)

#parameter tuning
parameters = {'n_estimators': [50, 100],
              'learning_rate': [0.1, 1],
              'random_state': [0]}

ab = AdaBoostClassifier()
clf = GridSearchCV(ab, parameters)
clf.fit(X_train, y_train)

print("Best Estimator:", clf.best_estimator_)
print("Train Accuracy:", clf.best_score_)

#predict and evaluate
preds = clf.best_estimator_.predict(X_test)
print("Test Accuracy:", accuracy_score(y_test, preds))
print("Test F1 score:", f1_score(y_test, preds, average='weighted'))
print("Test Balanced Accuracy:", balanced_accuracy_score(y_test, preds))



