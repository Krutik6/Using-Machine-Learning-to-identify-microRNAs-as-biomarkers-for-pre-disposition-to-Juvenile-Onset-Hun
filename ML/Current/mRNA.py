# -*- coding: utf-8 -*-
import os
from HDwork import *
structuring.structuringfunc()
multiclasstest.multiclasstestfunc()
recfs.recfsfunc()
kbest.kbestfunc()
ada.adafunc()
rf.rffunc()
extra.extrafunc()

#set global variables
WORKING_DIRECTORY = os.path.dirname(__file__)
DATA_DIR = os.path.join(WORKING_DIRECTORY, 'Data')
INFO_DIR = os.path.join(WORKING_DIRECTORY, 'Info')
PLOT_DIR = os.path.join(INFO_DIR, "Plot")
y_DIR = os.path.join(INFO_DIR, "besty")
Data = os.path.join(DATA_DIR, "10m_mRNA.csv")
ValData = os.path.join(DATA_DIR, "2m_mRNA.csv")
#read in and sort data
df10, df2 = structuring.readData(Data, ValData)
df10train, df10test, df2train, df2test = structuring.tentwo(df10, df2)
################################################################################
#feature selection based on recursive feature elimination
################################################################################
#try only with 10m
##Feature Selection
# supp = recfs.reccySVC(df10train, df10test, minfeat=5)
# rfe10train,rfe10test = recfs.retreiveFeat(df10train, df10test, supp)
#45 features found

##Train Multiple Classes
# trainclasses10, testclasses10 = multiclasstest.mulitanalysis2(rfe10train, rfe10test, 5)
#multiclasstest.makeconfusions(rfe10train, rfe10test)
#ADAboost, Neural Net, polynomial SVM, Guassian Process best
################################################################################
#try only with 2m
##Feature Selection
# supp = recfs.reccySVC(df2train, df2test, minfeat=5)
# rfe2train,rf2etest = recfs.retreiveFeat(df2train, df2test, supp)
#185 features found

##Train Multiple Classes
# trainclasses2, testclasses2 = multiclasstest.mulitanalysis2(rfe2train,rf2etest, 5)
#multiclasstest.makeconfusions(rfe2train,rf2etest)
#RF is best
################################################################################
#recursive feature elimination 
#use old data to predict young data
##Feature Selection
# supp = recfs.reccySVC(df10, df2, minfeat=5)
# rfeagetrain,rfeagetest = recfs.retreiveFeat(df10, df2, supp)
#24 features found

##Train Multiple Classes
# trainclassesage, testclassesage = multiclasstest.mulitanalysis2(rfeagetrain,rfeagetest, 5)
#multiclasstest.makeconfusions(rfeagetrain,rfeagetest)
#Linear SVM was best
################################################################################
#kbest only using ADA, Extratrees, NGB and RF - because they were best for miRNA
###############################################################################
#try only with 10m data
##Feature Selection
# ada10, gnb10, et10, rf10 = kbest.findbestvalk(df10train, df10test, reps=100)
###############################################################################
#try only with 2m data
##Feature Selection
# ada2, gnb2, et2, rf2 = kbest.findbestvalk(df2train, df2test, reps=100)
################################################################################
#use old data to predict young data - pretty poor
##Feature Selection
# ada, gnb, et, rf = kbest.findbestvalk(df10, df2, reps=100)
###############################################################################
#Hypertuning parameters
#10 - 7, 8, 11, 25, 27, 31
#2 - 17, 25, 26, 36, 52, 79, 87
#age - 4, 31, 47, 70, 76, 89
###############################################################################
#For ADAboost get best params and training and testing results
# f10 = kbest.findkbest(df10train, df10test, feats=31)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# ada.gridne(k10train, k10test) 
#ada.gridlr(k10train, k10test, ne=50) 
# ada.predict(k10train, k10test, lr=0.01, ne=50)
#7 0.96, 0.83, 9,0,2,1, 500, 0.01
#8, 0.91, 0.83, 8,1,1,2, 100, 1
#9, 0.96, 0.83, 9,0,2,1, 10, 0.1
#11, 0.96, 0.83, 9,0,2,1, 100, 0.01 
#25, 0.95, 0.83, 9,0,2,1, 50, 0.01
#27, 0.95, 0.83, 9,0,2,1, 50, 0.01
#31,  0.95, 0.83, 9,0,2,1, 50, 0.01 
#36, 0.93, 0.83, 9,0,2,1, 100, 0.01

# f2 = kbest.findkbest(df2train, df2test, feats=87)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
# ada.gridne(k2train, k2test)
# ada.gridlr(k2train, k2test, ne=500)
# ada.predict(k2train, k2test, lr=1, ne=500)
#3, 0.93, 0.63, 5,2,2,2, 10, 0.1
#17, 0.92, 0.72, 6,1,2,2, 50, 0.01
#25, 0.90, 0.81, 6,1,1,3, 10, 1
#26, 0.90, 0.81, 6,1,1,3, 10, 1
#36, 0.89, 0.81,  6,1,1,3, 100, 1
#52,  0.93, 0.81, 6,1,2,2, 10, 1
#79, 0.93, 0.72, 6,1,2,2, 500,1
#87, 0.93, 0.90, 6,1,0,4, 500, 1

# fage = kbest.findkbest(df10, df2, feats=89)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# ada.gridne(kagetrain, kagetest) 
# ada.gridlr(kagetrain, kagetest, ne=100)
# ada.predict(kagetrain, kagetest, lr=1, ne=100)
#4, 0.92, 0.8, 36,3,8,8, 10, 0.1
#10, 0.87, 0.76, 36,4,10,6, 50, 1
#31, 0.82, 0.81, 38,1,9,7, 10, 0.1
#47, 0.85, 0.74, 37,2,12,4, 50,1
#70, 0.8, 0.72, 39,0,15,1, 50, 1
#76 0.9, 0.74, 39,0,14,2, 100, 1
#89, 0.8, 0.70, 39,0,15,1, 500, 0.01
###############################################################################
#For GaussianNegativeBinomial get best params and training and testing results
# f10 = kbest.findkbest(df10train, df10test, feats=8)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# gnb.predict(k10train, k10test) 
#7, 0.96, 0.91, 8,1,0,3, smoother
#8, 0.96, 1, 9,0,3,0, smoother
#9, 0.96, 1, 9,0,3,0, smoother
#11, 0.96, 0.91, 9,0,1,2 smoother
#25, 0.96, 0.91, 9,0,1,2, smoother
#27, 0.96, 0.91, 9,0,1,2, smoother
#31, 0.96, 1, 1,9,0,0,3, smoother
#36, 0.98, 0.91, 9,0,1,2, no smoother

# f2 = kbest.findkbest(df2train, df2test, feats=79)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
# gnb.predict(k2train, k2test) 
#3, 0.92, 0.81, 5,2,0,4, same
#17, 1, 0.81, 7,0,2,2, no smoother
#25, 0.98, 0.96, 6,1,3,1, no smoother
#26, 0.96, 0.72, 6,1,2,2, smoother
#36, 0.96, 0.72, 6,1,2,2, smoother
#52, 0.96, 0.81, 6,1,1,3, smoother
#79, 0.95, 0.72, 6,1,2,2, smoother
#87, 95, 0.63, 6,1,3,1, same

# fage = kbest.findkbest(df10, df2, feats=76)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# gnb.predict(kagetrain, kagetest) 
#4, 0.95, 0.83, 33,6,3,13, smoother
#10, 0.96, 0.76, 39,0,13,3, no smoother 
#31, 0.96, 0.70, 39,0,16,0, same
#47, 0.97, 0.70, 39,0,16,0, no smoother
#70, 0.96, 0.70, 39,0,16,0, same
#76, 0.96, 0.70, 39,0,16,0, same
#89, 0.98, 0.70, 39,0,16,0, no smoother 
###############################################################################
#For ExtraTrees get best params and training and testing results
# f10 = kbest.findkbest(df10train, df10test, feats=36)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# extra.numtrees(k10train, k10test) 
# extra.numsplits(k10train, k10test, ne=1000) 
# extra.predict(k10train, k10test, ne=1000, ms=6)
#7, 0.98, 0.91, 9,0,1,2, 10,3
#8, 0.96, 0.91, 9,0,1,2, 100,14
#9, 0.96, 0.83, 9,0,2,1, 50,5
#11, 0.95, 0.83, 9,0,2,1, 50,3
#25, 0.96, 0.83, 9,0,2,1, 1000,9
#27, 0.96, 0.83, 9,0,1,2, 10,10
#31, 0.95, 1, 9,0,0,3, 50,4
#36, 0.96, 0.91, 9,0,1,2, 1000,3

# f2 = kbest.findkbest(df2train, df2test, feats=79)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
# extra.numtrees(k2train, k2test) 
# extra.numsplits(k2train, k2test, ne=100)
# extra.predict(k2train, k2test, ne=100, ms=7)
#3, 0.93, 0.63, 5,2,2,2, 100,6
#17, 0.93, 0.72, 6,1,2,2, 100,8
#25, 0.96, 0.72, 6,1,2,2, 100,11
#26, 0.96, 0.72, 6,1,2,2, 100,11
#36, 0.92, 0.81, 7,0,2,2, 100,11
#52, 0.95, 0.72, 7,0,3,1, 10,2
#79, 0.96, 0.63, 6,1,3,1, 100,7
#87, 0.90, 0.72, 6,1,2,2, 10,13

# fage = kbest.findkbest(df10, df2, feats=4)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# extra.numtrees(kagetrain, kagetest) 
# extra.numsplits(kagetrain, kagetest, ne=100) 
# extra.predict(kagetrain, kagetest, ne=500, ms=7)
#4, 0.92, 0.85, 37,2,6,10, 500,7
#10, 0.93, 0.74, 38,1,13,3, 500,5
#31, 0.96, 0.70, 39,0,16,0, 10,6
#47, 0.93, 0.72, 39,0,16,0, 1000,7
#70, 0.93, 0.70, 39,0,16,0, 1000,3
#76, 0.96, 0.70, 39,0,16,0, 100,9
#89, 0.91, 0.70, 39,0,16,0, 100,2
###############################################################################
#For RandomForest get best params and training and testing results
# f10 = kbest.findkbest(df10train, df10test, feats=36)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# rf.paramfind(k10train, k10test)
# rf.predict(k10train, k10test, ne=10, ml=1, ms=2, md=5)
#7, 0.96, 0.91, 9,0,1,2, 10
#8, 0.95, 0.91, 9,0,1,2, 50
#9, 1, 0.83, 9,0,2,1, 100
#11, 0.96, 0.83, 9,0,2,1, 50
#25, 0.91, 0.83, 9,0,2,1, 10
#27, 0.98, 0.83, 9,0,2,1, 10
#31, 0.95, 0.91, 9,0,1,2, 100
#36, 0.98, 0.83, 9,0,2,1, 10

# f2 = kbest.findkbest(df2train, df2test, feats=79)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
# rf.paramfind(k2train, k2test)
# rf.predict(k2train, k2test, ne=50, ml=1, ms=2, md=5)
#3, 0.90, 0.63, 5,2,2,2, 50
#17, 0.92, 0.72, 6,1,2,2,50 
#25, 0.90, 0.81, 6,1,1,3, 500
#26, 0.95, 0.72, 6,1,2,2, 10
#36, 0.92, 0.81, 6,1,1,3, 50
#52, 0.93, 0.81, 6,1,1,3, 50
#79, 0.92, 0.91, 7,0,2,2, 50
#87, 0.92, 0.81, 6,1,1,3, 100

# fage = kbest.findkbest(df10, df2, feats=76)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# rf.paramfind(kagetrain, kagetest)
# rf.predict(kagetrain, kagetest, ne=500, ml=1, ms=2, md=5)
#4, 0.90, 0.80, 36,3,8,8, 50
#10, 0.88, 0.78, 37,2,10,6, 50
#31, 0.90, 0.72, 37,2,13,3, 10
#47, 0.78, 0.74, 39,0,14,2, 10
#70, 0.75, 0.72, 39,0,15,1, 50
#76, 0.87, 0.72, 39,0,15,1, 500
#89, 0.81, 0.72, 39,0,15,1, 50
###############################################################################
#best cases
#Print ROC_CURVES
#10m
# f10 = kbest.findkbest(df10train, df10test, feats=8)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# m, X_train, y_train, X_test, y_test, y_pred = gnb.predict2(k10train, k10test) 
# plotROC(m, X_test, y_test, "mRNA Aged Samples")

#error occurs when plotting mRNA Aged Samples ROC - so did here
# SMALL_SIZE = 12
# MEDIUM_SIZE = 14
# BIGGER_SIZE = 17

# plotname="ROC - "+"mRNA Aged Samples"
# "mRNA Aged Samples".replace(" ", "_")
# savename="mRNA Aged Samples"+".png"
# plt.figure(figsize=(20,20))
# metrics.plot_roc_curve(m, X_test, y_test, pos_label="HD") 
# plt.title(plotname)
# plt.xlabel('False Positive Rate (Pos = JOHD)')
# plt.ylabel('True Positive Rate (Pos = JOHD)')
# plt.rc('font', size=MEDIUM_SIZE)         
# plt.rc('axes', titlesize=MEDIUM_SIZE)    
# plt.rc('axes', labelsize=MEDIUM_SIZE)   
# plt.rc('xtick', labelsize=SMALL_SIZE)   
# plt.rc('ytick', labelsize=SMALL_SIZE)  
# plt.rc('legend', fontsize=MEDIUM_SIZE)   
# plt.rc('figure', titlesize=BIGGER_SIZE)  
# plt.plot([0, 1], [0, 1],'r--', alpha=0.7)
# plt.show()
#plt.savefig(os.path.join(PLOT_DIR, savename), dpi=300) 

#2m
# f2 = kbest.findkbest(df2train, df2test, feats=87)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
# m, X_train, y_train, X_test, y_test, y_pred = ada.predict2(k2train, k2test, lr=1, ne=500)
# plotbest.plotROC(m, X_test, y_test, "mRNA Young Samples")

#age
# fage = kbest.findkbest(df10, df2, feats=4)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# m, X_train, y_train, X_test, y_test, y_pred = extra.predict2(kagetrain, kagetest, ne=500, ms=7)
# plotbest.plotROC(m, X_test, y_test, "mRNA Predisposed Samples")
################################################################################
# MAKE PCAs - save data and make in R
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

f10 = kbest.findkbest(df10train, df10test, feats=8)
k10train, k10test = kbest.ktodf(f10, df10train, df10test)
m, X_train, y_train, X_test, y_test, y_pred = gnb.predict2(k10train, k10test) 
plotbest.forR(X_test, y_test, y_pred, "mRNA_Aged.csv", f10)

f2 = kbest.findkbest(df2train, df2test, feats=87)
k2train, k2test = kbest.ktodf(f2, df2train, df2test)
m, X_train, y_train, X_test, y_test, y_pred = ada.predict2(k2train, k2test, lr=1, ne=500)
plotbest.forR(X_test, y_test, y_pred, "mRNA_Young.csv", f2)

fage = kbest.findkbest(df10, df2, feats=4)
kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
m, X_train, y_train, X_test, y_test, y_pred = extra.predict2(kagetrain, kagetest, ne=500, ms=7)
plotbest.forR(X_test, y_test, y_pred, "mRNA_Predisposed.csv", fage)