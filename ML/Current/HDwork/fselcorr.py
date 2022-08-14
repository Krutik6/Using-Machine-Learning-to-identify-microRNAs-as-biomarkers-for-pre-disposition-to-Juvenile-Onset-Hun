# -*- coding: utf-8 -*-
import os, sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import structuring

from copy import deepcopy
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
from sklearn.feature_selection import f_regression
import numpy as np
from collinearity import SelectNonCollinear 

#global parameters
PLOT_DIR = 'C:\\Users\\nkp68\\OneDrive - Newcastle University\\Writing_Papers\\HD\\ML\\Currect\\Info\\Plot'

def fselectionfunc():
  print("fselection module loaded")
  
  
def createCorr(dftrain):
    corr10 = deepcopy(dftrain)
    corr10.Samples.replace(('HD', 'WT'), (1, 0), inplace=True)
    corr10 = corr10.corr(method='spearman')
    f, ax = plt.subplots(figsize=(32, 32))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    heatmap = sns.heatmap(corr10, cmap=cmap, center=0.0, vmax=1,
                           linewidths=0, ax=ax)
    ax.invert_yaxis()
    ax.set_xlabel('microRNAs')
    ax.set_ylabel('microRNAs')
    plt.savefig(os.path.join(PLOT_DIR, 'hmap_10.png')) 
    plt.close()
    return corr10
    

def correlation(dataset, threshold):
    col_corr = set()  # Set of all the names of correlated columns
    corr_matrix = dataset.corr()
    for i in range(len(corr_matrix.columns)):
        for j in range(i):
            if abs(corr_matrix.iloc[i, j]) > threshold:  # we are interested in absolute coeff value
                colname = corr_matrix.columns[i]  # getting the name of column
                col_corr.add(colname)
    return col_corr


def dropHighCorr(dftrain, dftest, corrf):
    corrFE10 = dftrain.drop(corrf, axis=1)
    corrFE2 = dftest.drop(corrf, axis=1)
    Samples = dftrain.iloc[:,-1]
    corrFE10 = pd.concat([corrFE10, Samples], axis=1)
    corrFE2 = pd.concat([corrFE2, Samples], axis=1)
    return corrFE10, corrFE2


def colin(dftrain, dftest, threshold=0.4):
    #copy training and split into X and y
    copytrain = deepcopy(dftrain)
    copytrain.Samples.replace(('HD', 'WT'), (1, 0), inplace=True)
    features = list(copytrain)
    X, y, X_test, y_test = structuring.toArrays(copytrain, dftest)
    X = X.astype(np.int64)
    y = y.astype(np.int64)
    df = pd.DataFrame(X,columns=features[:-1])
    #plot initial heatmap
    f, ax = plt.subplots(figsize=(12, 8))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    heatmap = sns.heatmap(df, cmap=cmap, center=0.0, vmax=1,
                           linewidths=0, ax=ax)
    ax.invert_yaxis()
    ax.set_xlabel('microRNAs')
    ax.set_ylabel('microRNAs')
    plt.savefig(os.path.join(PLOT_DIR, 'hmap_colin.png')) 
    plt.close()
    #apply mask to remove collinear data
    selector = SelectNonCollinear(correlation_threshold=threshold,scoring=f_regression)
    selector.fit(X,y)
    mask = selector.get_support()
    df2 = pd.DataFrame(X[:,mask],columns = np.array(features[:-1])[mask])
    #make new plot
    heatmap = sns.heatmap(df2, cmap=cmap, center=0.0, vmax=1,
                           linewidths=0, ax=ax)
    ax.invert_yaxis()
    ax.set_xlabel('microRNAs')
    ax.set_ylabel('microRNAs')
    plt.savefig(os.path.join(PLOT_DIR, 'hmap_colin_maskapplied.png')) 
    plt.close()
    #save new pd objects
    new_train = pd.DataFrame(X[:,mask],columns = np.array(features[:-1])[mask])
    new_train['Samples'] = dftrain.Samples
    new_test = pd.DataFrame(X_test[:,mask],columns = np.array(features[:-1])[mask])
    new_test['Samples'] = dftrain.Samples
    
    return new_train, new_test