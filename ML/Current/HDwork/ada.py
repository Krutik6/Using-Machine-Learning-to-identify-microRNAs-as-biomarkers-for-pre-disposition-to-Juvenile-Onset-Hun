# -*- coding: utf-8 -*-
import os, sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import structuring


from numpy import mean
from numpy import std
from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.ensemble import AdaBoostClassifier
from matplotlib import pyplot
from sklearn.tree import DecisionTreeClassifier
from numpy import arange
from sklearn.model_selection import GridSearchCV
import pandas as pd
from sklearn.model_selection import KFold 
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score
from sklearn import metrics
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score

#global parameters
PLOT_DIR = 'C:\\Users\\nkp68\\OneDrive - Newcastle University\\Writing_Papers\\HD\\ML\\Currect\\Info\\Plot'

def adafunc():
  print("ada module loaded")
  

def gridlr(traindf, testdf, ne=500):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
    
    model = AdaBoostClassifier(n_estimators=ne)
    grid = dict()
    grid['learning_rate'] = [0.0001, 0.001, 0.01, 0.1, 1.0]
    cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
    grid_search = GridSearchCV(estimator=model, param_grid=grid, n_jobs=-1, 
                               cv=cv, scoring='accuracy')
    grid_result = grid_search.fit(X_train, y_train)
    print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))
    means = grid_result.cv_results_['mean_test_score']
    stds = grid_result.cv_results_['std_test_score']
    params = grid_result.cv_results_['params']
    for mean, stdev, param in zip(means, stds, params):
        print("%f (%f) with: %r" % (mean, stdev, param))
        
    
def gridne(dftrain, dftest):
    X_train, y_train, X_test, y_test = structuring.toArrays(dftrain, dftest)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
    
    model = AdaBoostClassifier()
    grid = dict()
    grid['n_estimators'] = [10, 50, 100, 500, 1000]
    cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
    grid_search = GridSearchCV(estimator=model, param_grid=grid, n_jobs=-1, 
                               cv=cv, scoring='accuracy')
    grid_result = grid_search.fit(X_train, y_train)
    print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))
    means = grid_result.cv_results_['mean_test_score']
    stds = grid_result.cv_results_['std_test_score']
    params = grid_result.cv_results_['params']
    for mean, stdev, param in zip(means, stds, params):
        print("%f (%f) with: %r" % (mean, stdev, param))
        
        

def predict(traindf, testdf, lr, ne=500):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    TrainAcc = []
    TestAcc = []
    
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)

    k = 5
    kf = KFold(n_splits=k, random_state=None)
    model = AdaBoostClassifier(n_estimators=ne, learning_rate=lr)

     
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
    
def predict2(traindf, testdf, lr, ne=500):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    TrainAcc = []
    TestAcc = []
    
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)

    k = 5
    kf = KFold(n_splits=k, random_state=None)
    model = AdaBoostClassifier(n_estimators=ne, learning_rate=lr)

    result = cross_val_score(model , X_train, y_train, cv = kf)
     
    model.fit(X_train, y_train)
    expected = y_test
    predicted = model.predict(X_test)
    
    print(accuracy_score(expected, predicted), ": is the accuracy score")
    
    return model, X_train, y_train, X_test, y_test, predicted
