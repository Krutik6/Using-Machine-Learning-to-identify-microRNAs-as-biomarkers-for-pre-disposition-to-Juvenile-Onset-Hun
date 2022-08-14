# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import os, sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import structuring

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold 
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.model_selection import cross_val_score

#global parameters
PLOT_DIR = 'C:\\Users\\nkp68\\OneDrive - Newcastle University\\Writing_Papers\\HD\\ML\\Attempt\\Info\\Plot'

def rffunc():
  print("rf module loaded")
  

def paramfind(traindf, testdf):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)

    n_estimators = [10, 50, 100, 500, 1000]
    forest = RandomForestClassifier(random_state = 1)
    hyperF = dict(n_estimators = n_estimators)

    gridF = GridSearchCV(forest, hyperF, cv = 3, verbose = 1, n_jobs = -1)
    bestF = gridF.fit(X_train, y_train)
    print(bestF.best_params_)
    

def predict(traindf, testdf, ne, md=5, ms=2, ml=1):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
    
    k = 5
    kf = KFold(n_splits=k, random_state=None)
    model = RandomForestClassifier(random_state = 1, n_estimators=ne, 
                                   max_depth=md, min_samples_split=ms,
                                   min_samples_leaf=ml)

    result = cross_val_score(model , X_train, y_train, cv = kf)
     
    print("Avg accuracy: {}".format(result.mean()))
    model.fit(X_train, y_train)
    expected = y_test
    predicted = model.predict(X_test)

    print(confusion_matrix(expected, predicted), ": is the confusion matrix")
    print(accuracy_score(expected, predicted), ": is the accuracy score")
    print(precision_score(expected, predicted, pos_label='HD'), ": is the precision score")
    print(recall_score(expected, predicted, pos_label='HD'), ": is the recall score")
    print(f1_score(expected, predicted, pos_label='HD'), ": is the f1 score")
    

def predict2(traindf, testdf, ne, md=5, ms=2, ml=1):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
    
    k = 5
    kf = KFold(n_splits=k, random_state=None)
    model = RandomForestClassifier(random_state = 1, n_estimators=ne, 
                                   max_depth=md, min_samples_split=ms,
                                   min_samples_leaf=ml)

    result = cross_val_score(model , X_train, y_train, cv = kf)
     
    model.fit(X_train, y_train)
    expected = y_test
    predicted = model.predict(X_test)
    
    print(accuracy_score(expected, predicted), ": is the accuracy score")
    return model, X_train, y_train, X_test, y_test, predicted

