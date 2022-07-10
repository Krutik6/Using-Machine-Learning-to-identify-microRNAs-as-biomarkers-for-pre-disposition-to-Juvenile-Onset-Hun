# -*- coding: utf-8 -*-
import os, sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import structuring

from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold 
import numpy as np

PLOT_DIR = 'C:\\Users\\nkp68\\OneDrive - Newcastle University\\Writing_Papers\\HD\\ML\\Currect\\Info\\Plot'

def gnbfunc():
  print("gnb module loaded")


def predict(traindf, testdf):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)

    param_grid_nb = {
        'var_smoothing': np.logspace(0,-9, num=100)
    }
    nbModel_grid = GridSearchCV(estimator=GaussianNB(), 
                                param_grid=param_grid_nb, verbose=1,
                                cv=5, n_jobs=-1)

    nbModel_grid.fit(X_train, y_train)
    print(nbModel_grid.best_estimator_)
    vs=nbModel_grid.best_estimator_.var_smoothing
    
    k = 5
    kf = KFold(n_splits=k, random_state=None)
    print("With variant best smoother parameter identified")
    model = GaussianNB(priors=None, var_smoothing=vs)
    result = cross_val_score(model , X_train, y_train, cv = kf)
    print("Avg accuracy: {}".format(result.mean()))
    
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    print(confusion_matrix(y_test, y_pred), ": is the confusion matrix")
    print(accuracy_score(y_test, y_pred), ": is the accuracy score")
    print(precision_score(y_test, y_pred, pos_label='HD'), ": is the precision score")
    print(recall_score(y_test, y_pred, pos_label='HD'), ": is the recall score")
    print(f1_score(y_test, y_pred, pos_label='HD'), ": is the f1 score")
    
    print("Without using best variant smoother")
    model = GaussianNB()
    result = cross_val_score(model , X_train, y_train, cv = kf)
    print("Avg accuracy: {}".format(result.mean()))
    
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    print(confusion_matrix(y_test, y_pred), ": is the confusion matrix")
    print(accuracy_score(y_test, y_pred), ": is the accuracy score")
    print(precision_score(y_test, y_pred, pos_label='HD'), ": is the precision score")
    print(recall_score(y_test, y_pred, pos_label='HD'), ": is the recall score")
    print(f1_score(y_test, y_pred, pos_label='HD'), ": is the f1 score")
    
    
def predict2(traindf, testdf):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)

    param_grid_nb = {
        'var_smoothing': np.logspace(0,-9, num=100)
    }
    nbModel_grid = GridSearchCV(estimator=GaussianNB(), 
                                param_grid=param_grid_nb, verbose=1,
                                cv=5, n_jobs=-1)

    nbModel_grid.fit(X_train, y_train)
    vs=nbModel_grid.best_estimator_.var_smoothing
    
    k = 5
    kf = KFold(n_splits=k, random_state=None)
    model = GaussianNB(priors=None, var_smoothing=vs)
    result = cross_val_score(model , X_train, y_train, cv = kf)
    
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    
    return model, X_train, y_train, X_test, y_test, y_pred
    


   
    
    