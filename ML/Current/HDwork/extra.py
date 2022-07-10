# -*- coding: utf-8 -*-
import os, sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import structuring

from numpy import mean
from numpy import std
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.model_selection import KFold 
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score

#global parameters
PLOT_DIR = 'C:\\Users\\nkp68\\OneDrive - Newcastle University\\Writing_Papers\\HD\\ML\\Currect\\Info\\Plot'

def extrafunc():
  print("extra module loaded")

def numtrees(traindf, testdf):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)

    def get_models():
    	models = dict()
    	# define number of trees to consider
    	n_trees = [10, 50, 100, 500, 1000]
    	for n in n_trees:
    		models[str(n)] = ExtraTreesClassifier(n_estimators=n)
    	return models

    def evaluate_model(model, X, y):
    	# define the evaluation procedure
    	cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
    	# evaluate the model and collect the results
    	scores = cross_val_score(model, X, y, scoring='accuracy', cv=cv, n_jobs=-1)
    	return scores

    models = get_models()

    results, names = list(), list()
    for name, model in models.items():
    	# evaluate the model
    	scores = evaluate_model(model, X_train, y_train)
    	# store the results
    	results.append(scores)
    	names.append(name)
    	# summarize the performance along the way
    	print('>%s %.3f (%.3f)' % (name, mean(scores), std(scores)))
        

def numsplits(traindf, testdf, ne=10):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
    # get a list of models to evaluate
    def get_models():
    	models = dict()
    	# explore the number of samples per split from 2 to 14
    	for i in range(2, 15):
    		models[str(i)] = ExtraTreesClassifier(n_estimators=ne, min_samples_split=i)
    	return models
     
    # evaluate a given model using cross-validation
    def evaluate_model(model, X, y):
    	# define the evaluation procedure
    	cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=None)
    	# evaluate the model and collect the results
    	scores = cross_val_score(model, X, y, scoring='accuracy', cv=cv, n_jobs=-1)
    	return scores
     
        
    # get the models to evaluate
    models = get_models()
    # evaluate the models and store results
    results, names = list(), list()
    for name, model in models.items():
    	# evaluate the model
    	scores = evaluate_model(model, X_train, y_train)
    	# store the results
    	results.append(scores)
    	names.append(name)
    	# summarize the performance along the way
    	print('>%s %.3f (%.3f)' % (name, mean(scores), std(scores)))
        

def predict(traindf, testdf, ne=10, ms=2):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)

    k = 5
    kf = KFold(n_splits=k, random_state=None)
    model = ExtraTreesClassifier(n_estimators=ne, min_samples_split=ms)

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
    
def predict2(traindf, testdf, ne=10, ms=2):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)

    k = 5
    kf = KFold(n_splits=k, random_state=None)
    model = ExtraTreesClassifier(n_estimators=ne, min_samples_split=ms)

    result = cross_val_score(model , X_train, y_train, cv = kf)
     
    model.fit(X_train, y_train)
    expected = y_test
    predicted = model.predict(X_test)
    
    print(accuracy_score(expected, predicted), ": is the accuracy score")

    return model, X_train, y_train, X_test, y_test, predicted