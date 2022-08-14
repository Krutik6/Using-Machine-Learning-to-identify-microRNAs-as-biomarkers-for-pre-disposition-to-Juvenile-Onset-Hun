# -*- coding: utf-8 -*-
import os, sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import sklearn.metrics as metrics
import numpy as np
import pandas as pd

#global parameters
PLOT_DIR = 'C:\\Users\\nkp68\\OneDrive - Newcastle University\\Writing_Papers\\HD\\ML\\Current\\Info\\Plot'
y_DIR = 'C:\\Users\\nkp68\\OneDrive - Newcastle University\\Writing_Papers\\HD\\ML\\Current\\Info\\besty'

def plotbestfunc():
  print("plotbest module loaded")
  
  
def plotROC(model, X_test, y_test, name):  
    SMALL_SIZE = 12
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 17
    
    plotname="AUC - "+name
    name.replace(" ", "_")
    savename=name+".png"
    plt.figure(figsize=(20,20))
    metrics.plot_roc_curve(model, X_test, y_test, pos_label="HD") 
    plt.title(plotname)
    plt.xlabel('False Positive Rate (Pos = JOHD)')
    plt.ylabel('True Positive Rate (Pos = JOHD)')
    plt.rc('font', size=MEDIUM_SIZE)         
    plt.rc('axes', titlesize=MEDIUM_SIZE)    
    plt.rc('axes', labelsize=MEDIUM_SIZE)   
    plt.rc('xtick', labelsize=SMALL_SIZE)   
    plt.rc('ytick', labelsize=SMALL_SIZE)  
    plt.rc('legend', fontsize=MEDIUM_SIZE)   
    plt.rc('figure', titlesize=BIGGER_SIZE)  
    plt.plot([0, 1], [0, 1],'r--', alpha=0.7)
    plt.savefig(os.path.join(PLOT_DIR, savename), dpi=300) 
    plt.close()


def forR(X_test, y_test, y_pred, name, f):
    y = np.c_[y_test,y_pred]
    y = np.c_[X_test, y]
    y = pd.DataFrame(y)
    y = y.replace('HD' ,'JOHD')
    newnames = f.iloc[:,0]
    ser = pd.Series(['True Labels', 'Predicted Labels'])
    newnames = newnames.append(ser)
    y.columns = [newnames]
    y.to_csv(os.path.join(y_DIR, name), sep=',')