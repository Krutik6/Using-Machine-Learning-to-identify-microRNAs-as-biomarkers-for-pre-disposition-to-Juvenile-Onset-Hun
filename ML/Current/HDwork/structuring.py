# -*- coding: utf-8 -*-
import os, sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import structuring
from sklearn.model_selection import train_test_split
import pandas as pd
from sklearn import preprocessing
from sklearn import metrics
import matplotlib.pyplot as plt 
import os
from imblearn.over_sampling import RandomOverSampler
from collections import Counter

#global parameters
PLOT_DIR = 'C:\\Users\\nkp68\\OneDrive - Newcastle University\\Writing_Papers\\HD\\ML\\Current\\Info\\Plot'

def structuringfunc():
  print("structuring module loaded")

def readData(path, valpath):
    df10 = pd.read_csv(path)
    df10 = df10.drop(df10.columns[[0]], axis=1)
    print(df10.shape)
    df2 = pd.read_csv(valpath)
    df2 = df2.drop(df2.columns[[0]], axis=1)
    print(df2.shape)
    return df10, df2


def splitsies(train, test):
    X10 = train.loc[:, train.columns != 'Samples']
    y10 = train.loc[:, train.columns == 'Samples']
    X2 = test.loc[:, test.columns != 'Samples']
    y2 = test.loc[:, test.columns == 'Samples']
    return X10, y10, X2, y2


def scale1to0(x10, x2):
    plt.plot(x10)
    plt.savefig(os.path.join(PLOT_DIR, 'unscaled_x10.png'))
    plt.close()
    scaler = preprocessing.MinMaxScaler().fit(x10)
    x10_scaled = scaler.transform(x10)
    plt.plot(x10_scaled)
    plt.savefig(os.path.join(PLOT_DIR, 'scaled_x10.png'))
    plt.close()
    plt.plot(x2)
    plt.savefig(os.path.join(PLOT_DIR, 'unscaled_x2.png'))
    plt.close()
    scaler = preprocessing.MinMaxScaler().fit(x2)
    x2_scaled = scaler.transform(x2)
    plt.plot(x2_scaled)
    plt.savefig(os.path.join(PLOT_DIR, 'scaled_x2.png'))
    plt.close()
    return x10_scaled, x2_scaled


def toArrays(FE10, FE2):
    array10 = FE10.values
    X_train = array10[:, :-1]
    y_train = array10[:, -1]
    array2 = FE2.values
    X_test = array2[:, :-1]
    y_test = array2[:, -1]
    return X_train, y_train, X_test, y_test


def smotey(X_train, y_train):
  ROS = RandomOverSampler(random_state=42)
  x_train_oversampled, y_train_oversampled = ROS.fit_resample(X_train, y_train)
  return x_train_oversampled, y_train_oversampled
  

def tentwo(df10, df2, size=0.2):
    df10x, df10y, df2x, df2y = structuring.splitsies(df10, df2)
    
    X_train, X_test, y_train, y_test = train_test_split(df10x, df10y, test_size=size, random_state=42)
    df10train = pd.concat([X_train, y_train], axis=1)
    df10test = pd.concat([X_test, y_test], axis=1)

    X_train, X_test, y_train, y_test = train_test_split(df2x, df2y, test_size=size, random_state=42)
    df2train = pd.concat([X_train, y_train], axis=1)
    df2test = pd.concat([X_test, y_test], axis=1)
    
    return df10train, df10test, df2train, df2test

def getmean(scores):
    sums = Counter()
    counters = Counter()
    for itemset in [scores[0], scores[1], scores[2]]:
        sums.update(itemset)
        counters.update(itemset.keys())
        
    ret = {x: float(sums[x])/counters[x] for x in sums.keys()}

    return ret  