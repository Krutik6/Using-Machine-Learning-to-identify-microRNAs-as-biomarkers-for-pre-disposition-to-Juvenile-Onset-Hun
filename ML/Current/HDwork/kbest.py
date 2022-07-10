# -*- coding: utf-8 -*-
import os, sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import structuring
import kbest

import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import KFold 
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score
from sklearn import metrics
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import StratifiedKFold
import numpy as np

#global parameters
PLOT_DIR = 'C:\\Users\\nkp68\\OneDrive - Newcastle University\\Writing_Papers\\HD\\ML\\Current\\Info\\Plot'

def kbestfunc():
  print("kbest module loaded")

def findkbest(dftrain, dftest, feats=15):
    X_train, y_train, X_test, y_test = structuring.toArrays(dftrain, dftest)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
    
    selector = SelectKBest(f_classif, k = feats)
    X_new = selector.fit_transform(X_train, y_train)
    x_names = dftrain.iloc[:,:-1]
    names = x_names.columns.values[selector.get_support()]
    scores = selector.scores_[selector.get_support()]
    names_scores = list(zip(names, scores))
    ns_df = pd.DataFrame(data = names_scores, columns=['Feat_names', 'F_Scores'])
    ns_df_sorted = ns_df.sort_values(['F_Scores', 'Feat_names'], ascending = [False, True])
    print(ns_df_sorted)
    
    return ns_df_sorted


def ktodf(fdf, dftrain, dftest):
    arr = fdf.iloc[:,0]
    arr = arr.append(pd.Series(["Samples"]))
    dftrain_drop = dftrain[dftrain.columns.intersection(arr)]
    dftest_drop = dftest[dftest.columns.intersection(arr)]
    
    return dftrain_drop, dftest_drop


def kopt(traindf, testdf, maxfeat=100):
    
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
    X_train_copy = pd.DataFrame(data=X_train, columns=traindf.columns[:-1])
    rows, columns = X_train.shape
    NS = {}
    ada = {}
    gnb = {}
    et = {}
    rf = {}
    
    for i in range(maxfeat):
        kbest=i+1
        selector = SelectKBest(f_classif, k = kbest)
        X_new = selector.fit_transform(X_train, y_train)
        x_names = traindf.iloc[:,:-1]
        names = x_names.columns.values[selector.get_support()]
        scores = selector.scores_[selector.get_support()]
        names_scores = list(zip(names, scores))
        ns_df = pd.DataFrame(data = names_scores, columns=['Feat_names', 'F_Scores'])
        ns_df_sorted = ns_df.sort_values(['F_Scores', 'Feat_names'], ascending = [False, True])
        NS[i+1] = ns_df_sorted
        
        arr = ns_df_sorted.iloc[:,0]
        traindf_drop = X_train_copy[X_train_copy.columns.intersection(arr)]      
        #test accuracy of k features with adaboost
        kf = KFold(n_splits=5, random_state=None)
        model = AdaBoostClassifier(n_estimators=10)
        result = cross_val_score(model , traindf_drop, y_train, cv = kf)
        meanres = format(result.mean())
        ada[i+1] = meanres
        
        model = GaussianNB()
        result = cross_val_score(model , traindf_drop, y_train, cv = kf)
        meanres = format(result.mean())
        gnb[i+1] = meanres
        
        model = ExtraTreesClassifier(n_estimators=10, min_samples_split=2)
        result = cross_val_score(model , traindf_drop, y_train, cv = kf)
        meanres = format(result.mean())
        et[i+1] = meanres
        
        model = RandomForestClassifier(max_depth=10, n_estimators=10)
        result = cross_val_score(model , traindf_drop, y_train, cv = kf)
        meanres = format(result.mean())
        rf[i+1] = meanres
        
    return ada, gnb, et, rf, NS


def ktop(knames, score_ada, score_gnb, score_et, score_rf, w="Samples"):
    bestrow = max(score_ada, key=score_ada.get)
    bestscores = score_ada[bestrow]
    bestnames = knames[bestrow]
    print("Using ADABoost: For %s Best training found at %s using %s features" % (w, bestscores,
                                                                  bestrow))
    
    bestrow = max(score_gnb, key=score_gnb.get)
    bestscores = score_gnb[bestrow]
    bestnames = knames[bestrow]
    print("Using GaussianNB: For %s Best training found at %s using %s features" % (w, bestscores,
                                                                  bestrow))
    
    bestrow = max(score_et, key=score_et.get)
    bestscores = score_et[bestrow]
    bestnames = knames[bestrow]
    print("Using ExtraTrees: For %s Best training found at %s using %s features" % (w, bestscores,
                                                                  bestrow))
    
    bestrow = max(score_rf, key=score_rf.get)
    bestscores = score_rf[bestrow]
    bestnames = knames[bestrow]
    print("Using RandomForest: For %s Best training found at %s using %s features" % (w, bestscores,
                                                                  bestrow))
    


def kopt2(dftrain, dftest,m=1):
    X_train, y_train, X_test, y_test = structuring.toArrays(df10train, df10test)
    scores_ada = []
    scores_gnb = []
    scores_et = []
    scores_rf = []
    scores_names = []
    cv = StratifiedKFold(n_splits=3, random_state=42, shuffle=True)
    
    for train_idx, test_idx, in cv.split(X_train, y_train):
      X_train, y_train = structuring.smotey(X_train, y_train)
      X_train, X_test = structuring.scale1to0(X_train, X_test)
      X_train_copy = pd.DataFrame(data=X_train, columns=df10train.columns[:-1])
    
      rows, columns = X_train.shape
      ada = {}
      gnb = {}
      et = {}
      rf = {}
      NS = {}
           
      for i in range(m,100):
        kbest=i+1
        selector = SelectKBest(f_classif, k = kbest)
        X_new = selector.fit_transform(X_train, y_train)
        x_names = df10train.iloc[:,:-1]
        names = x_names.columns.values[selector.get_support()]
        scores = selector.scores_[selector.get_support()]
        names_scores = list(zip(names, scores))
        ns_df = pd.DataFrame(data = names_scores, columns=['Feat_names', 'F_Scores'])
        ns_df_sorted = ns_df.sort_values(['F_Scores', 'Feat_names'], ascending = [False, True])
        NS[i+1] = ns_df_sorted
          
        arr = ns_df_sorted.iloc[:,0]
        traindf_drop = X_train_copy[X_train_copy.columns.intersection(arr)]      
          
        model = AdaBoostClassifier(n_estimators=10)
        model.fit(traindf_drop, y_train)
        result = model.score(traindf_drop, y_train)
        ada[i+1] = result
          
        model = GaussianNB()
        model.fit(traindf_drop, y_train)
        result = model.score(traindf_drop, y_train)
        gnb[i+1] = result
          
        model = ExtraTreesClassifier(n_estimators=10, min_samples_split=2)
        model.fit(traindf_drop, y_train)
        result = model.score(traindf_drop, y_train)
        et[i+1] = result
          
        model = RandomForestClassifier(max_depth=10, n_estimators=10)
        model.fit(traindf_drop, y_train)
        result = model.score(traindf_drop, y_train)
        rf[i+1] = result
             
      scores_ada.append(ada)
      scores_gnb.append(gnb)
      scores_et.append(et)
      scores_rf.append(rf)
      scores_names.append(NS)
      namelist = scores_names[:1]
      namedict=namelist[0]

    ada = structuring.getmean(scores_ada)
    gnb = structuring.getmean(scores_gnb)
    et = structuring.getmean(scores_et)
    rf = structuring.getmean(scores_rf)

    return ada, gnb, et, rf, namedict



def findbestvalk(dftrain, dftest, reps=10):
    a = dict()
    a["ada"] = []
    a["gnb"] = []
    a["et"] = []
    a["rf"] = []

    for i in range(reps):
        ada_opt, gnb_opt, et_opt, rf_opt, names_opt = kbest.kopt(dftrain, dftest)
        a["ada"].append(ada_opt)
        a["gnb"].append(gnb_opt)
        a["et"].append(et_opt)
        a["rf"].append(rf_opt)

    ada = a["ada"]
    ada_df = pd.DataFrame(ada)
    ada_df = ada_df.astype('float64') 
    ada_maxiter = ada_df.idxmax(axis=1, skipna=True)
    ada_maxval = ada_df.max(axis=1)
    ada_concat  = pd.concat([ada_maxiter,ada_maxval], axis=1)
    
    gnb = a["gnb"]
    gnb_df = pd.DataFrame(gnb)
    gnb_df = gnb_df.astype('float64') 
    gnb_maxiter = gnb_df.idxmax(axis=1, skipna=True)
    gnb_maxval = gnb_df.max(axis=1)
    gnb_concat  = pd.concat([gnb_maxiter,gnb_maxval], axis=1)
    
    et = a["et"]
    et_df = pd.DataFrame(et)
    et_df = et_df.astype('float64') 
    et_maxiter = et_df.idxmax(axis=1, skipna=True)
    et_maxval = et_df.max(axis=1)
    et_concat  = pd.concat([et_maxiter,et_maxval], axis=1)
    
    rf = a["rf"]
    rf_df = pd.DataFrame(rf)
    rf_df = rf_df.astype('float64') 
    rf_maxiter = rf_df.idxmax(axis=1, skipna=True)
    rf_maxval = rf_df.max(axis=1)
    rf_concat  = pd.concat([gnb_maxiter,rf_maxval], axis=1)
    
    order = ["Highest", "Second Highest", "Third Highest"]
    
    ada_concat.sort_values(by=[1], inplace=True)
    for i in range(3):
       print(f"For ADABoost : {order[i]} accuracy was {ada_concat[1].loc[i]} with {ada_concat[0].loc[i]} features")
    
    gnb_concat.sort_values(by=[1], inplace=True)
    for i in range(3):
       print(f"For GaussianNB : {order[i]} accuracy was {gnb_concat[1].loc[i]} with {gnb_concat[0].loc[i]} features")
       
    et_concat.sort_values(by=[1], inplace=True)
    for i in range(3):
       print(f"For ExtraTrees : {order[i]} accuracy was {et_concat[1].loc[i]} with {et_concat[0].loc[i]} features")   
    
    rf_concat.sort_values(by=[1], inplace=True)
    for i in range(3):
       print(f"For RandomForest : {order[i]} accuracy was {rf_concat[1].loc[i]} with {rf_concat[0].loc[i]} features")
    
    return ada_concat, gnb_concat, et_concat, rf_concat