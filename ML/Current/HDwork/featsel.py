# -*- coding: utf-8 -*-
import os, sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import structuring

from sklearn.feature_selection import mutual_info_classif as MIC
import numpy as np
from sklearn.tree import DecisionTreeClassifier as DTC

#global parameters
PLOT_DIR = 'C:\\Users\\nkp68\\OneDrive - Newcastle University\\Writing_Papers\\HD\\ML\\Currect\\Info\\Plot'

def featselfun():
  print("featsel module loaded")


def MIclass(traindf, testdf, score=0.2):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
    
    mi_score = MIC(X_train, y_train)
    
    mi_score_selected_index = np.where(mi_score >score)[0]
    X_train_02 = X_train[:,mi_score_selected_index]
    X_test_02 = X_test[:,mi_score_selected_index]
    
    model_1 = DTC().fit(X_train,y_train)
    model_2 = DTC().fit(X_train_02,y_train)
    score_1 = model_1.score(X_test,y_test)
    score_2 = model_2.score(X_test_02,y_test)
    print(f"score_1:{score_1}\n score_2:{score_2}")