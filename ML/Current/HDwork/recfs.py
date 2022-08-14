# -*- coding: utf-8 -*-
import os, sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import structuring

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.datasets import make_classification
from copy import deepcopy
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE
from numpy import mean
from numpy import std
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.pipeline import Pipeline
from matplotlib import pyplot
from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.linear_model import Perceptron

#global parameters
PLOT_DIR = 'C:\\Users\\nkp68\\OneDrive - Newcastle University\\Writing_Papers\\HD\\ML\\Current\\Info\\Plot'

def recfsfunc():
  print("recfs module loaded")


def reccySVC(traindf, testdf, minfeat=1, strats=10):
    #set up data
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
    features = list(traindf)
    X_train = pd.DataFrame(X_train, columns = np.array(features[:-1]))
    X_test = pd.DataFrame(X_test, columns = np.array(features[:-1]))
    #perform svc stratified
    svc = SVC(kernel="linear")
    min_features_to_select = minfeat
    rfecv = RFECV(
        estimator=svc,
        step=1,
        cv=StratifiedKFold(strats),
        scoring="accuracy",
        min_features_to_select=min_features_to_select,)
    rfe = rfecv.fit(X_train, y_train)
    print("Optimal number of features : %d" % rfe.n_features_)
   
    plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Cross validation score (accuracy)")
    plt.plot(
    range(min_features_to_select, len(rfecv.grid_scores_) + min_features_to_select),
        rfe.grid_scores_,)
    plt.savefig(os.path.join(PLOT_DIR, 'RFE_skfolds_SVC.png')) 
    plt.close()
    
    rfesupp = rfe.support_
    
    return rfesupp


def retreiveFeat(dftrain, dftest, supp):
    training = deepcopy(dftrain)
    testing = deepcopy(dftest)
    Xtrain, ytrain, Xtest, ytest = structuring.splitsies(training, testing)
    Xbool = Xtrain.append(pd.DataFrame(supp.reshape(1,-1), columns=list(Xtrain)),
                     ignore_index=False)
    Xbool2 = Xbool.loc[:, Xbool.iloc[-1].ne(0)]
    XrfeTrain = Xbool2[:-1]
    rfelist = list(XrfeTrain.columns.values)
    XrfeTest = Xtest[Xtest.columns.intersection(rfelist)]
    training = XrfeTrain.join(ytrain)
    testing = XrfeTest.join(ytest)
    
    return training, testing


def reccyRF(traindf, testdf, minfeat=10, strats=5):
    #set up data
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    features = list(traindf)
    X_train = pd.DataFrame(X_train, columns = np.array(features[:-1]))
    X_test = pd.DataFrame(X_test, columns = np.array(features[:-1]))
    #perform svc stratified
    rf = RFECV(RandomForestClassifier(), cv=strats, scoring='f1_weighted')
    min_features_to_select = minfeat
    rfecv = RFECV(
        estimator=rf,
        step=1,
        cv=StratifiedKFold(strats),
        scoring="f1_weighted",
        min_features_to_select=min_features_to_select,)
    rfe = rfecv.fit(X_train, y_train)
    print("Optimal number of features : %d" % rfe.n_features_)
   
    plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Cross validation score (accuracy)")
    plt.plot(
    range(min_features_to_select, len(rfecv.grid_scores_) + min_features_to_select),
        rfe.grid_scores_,)
    plt.savefig(os.path.join(PLOT_DIR, 'RFE_skfolds_SVC.png')) 
    plt.close()
    
    rfesupp = rfe.support_
    rferank = rfe.ranking_ 
    rfeabs = np.absolute(rfe.estimator_.coef_)
    
    return rfesupp, rferank, rfeabs 


def boxplotRE(traindf, testdf):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
     
    # get a list of models to evaluate
    def get_models():
    	models = dict()
    	# lr
    	rfe = RFE(estimator=LogisticRegression(), n_features_to_select=5)
    	model = DecisionTreeClassifier()
    	models['lr'] = Pipeline(steps=[('s',rfe),('m',model)])
    	# perceptron
    	rfe = RFE(estimator=Perceptron(), n_features_to_select=5)
    	model = DecisionTreeClassifier()
    	models['per'] = Pipeline(steps=[('s',rfe),('m',model)])
    	# cart
    	rfe = RFE(estimator=DecisionTreeClassifier(), n_features_to_select=5)
    	model = DecisionTreeClassifier()
    	models['cart'] = Pipeline(steps=[('s',rfe),('m',model)])
    	# rf
    	rfe = RFE(estimator=RandomForestClassifier(), n_features_to_select=5)
    	model = DecisionTreeClassifier()
    	models['rf'] = Pipeline(steps=[('s',rfe),('m',model)])
    	# gbm
    	rfe = RFE(estimator=GradientBoostingClassifier(), n_features_to_select=5)
    	model = DecisionTreeClassifier()
    	models['gbm'] = Pipeline(steps=[('s',rfe),('m',model)])
    	return models
     
    # evaluate a give model using cross-validation
    def evaluate_model(model, X, y):
    	cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
    	scores = cross_val_score(model, X, y, scoring='accuracy', cv=cv, n_jobs=-1)
    	return scores
     
    # get the models to evaluate
    models = get_models()
    # evaluate the models and store results
    results, names = list(), list()
    for name, model in models.items():
    	scores = evaluate_model(model, X_train, y_train)
    	results.append(scores)
    	names.append(name)
    	print('>%s %.3f (%.3f)' % (name, mean(scores), std(scores)))
    # plot model performance for comparison
    pyplot.boxplot(results, labels=names, showmeans=True)
    pyplot.show()