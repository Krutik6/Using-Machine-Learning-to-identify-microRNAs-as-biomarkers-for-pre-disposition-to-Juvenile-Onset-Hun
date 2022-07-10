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
plotbest.plotbestfunc()

#set global variables
WORKING_DIRECTORY = os.path.dirname(__file__)
DATA_DIR = os.path.join(WORKING_DIRECTORY, 'Data')
INFO_DIR = os.path.join(WORKING_DIRECTORY, 'Info')
PLOT_DIR = os.path.join(INFO_DIR, "Plot")
Data = os.path.join(DATA_DIR, "10m_miRNA.csv")
ValData = os.path.join(DATA_DIR, "2m_miRNA.csv")
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
#224 features found

##Train Multiple Classes
# trainclasses10, testclasses10 = multiclasstest.mulitanalysis(rfe10train, rfe10test, 5)
#multiclasstest.makeconfusions(rfe10train, rfe10test)
#ExtraTrees, RF and ADA are best
################################################################################
#try only with 2m
##Feature Selection
# supp = recfs.reccySVC(df2train, df2test, minfeat=5)
# rfe2train,rf2etest = recfs.retreiveFeat(df2train, df2test, supp)
#86 features found

##Train Multiple Classes
# trainclasses2, testclasses2 = multiclasstest.mulitanalysis2(rfe2train,rf2etest, 3)
#multiclasstest.makeconfusions(rfe2train,rf2etest)
#ADA, NB and Extratees are best
################################################################################
#recursive feature elimination 
#use old data to predict young data - better than kbest
##Feature Selection
# supp = recfs.reccySVC(df10, df2, minfeat=5)
# rfeagetrain,rfeagetest = recfs.retreiveFeat(df10, df2, supp)
#50 features found

##Train Multiple Classes
# trainclassesage, testclassesage = multiclasstest.mulitanalysis2(rfeagetrain,rfeagetest, 3)
#multiclasstest.makeconfusions(rfeagetrain,rfeagetest)
#ADA, Extratrees and RF are best
################################################################################
#kbest only using ADA, Extratrees, NGB and RF
###############################################################################
#try only with 10m data
##Feature Selection
#ada10, gnb10, et10, rf10 = kbest.findbestvalk(df10train, df10test, reps=100)
###############################################################################
#try only with 2m data
##Feature Selection
#ada2, gnb2, et2, rf2 = kbest.findbestvalk(df2train, df2test, reps=100)
################################################################################
#use old data to predict young data - pretty poor
##Feature Selection
#ada, gnb, et, rf = kbest.findbestvalk(df10, df2, reps=100)
###############################################################################
#Hypertuning parameters
#10 - 20, 36, 44, 47 
#2 - 6, 37, 38, 52, 100  
#age - 15, 18, 41, 54
###############################################################################
#For ADAboost get best params and training and testing results
#f10 = kbest.findkbest(df10train, df10test, feats=36)
#k10train, k10test = kbest.ktodf(f10, df10train, df10test)
#ada.gridne(k10train, k10test) 
#ada.gridlr(k10train, k10test, ne=50) 
#ada.predict(k10train, k10test, lr=1, ne=50)
#8, 0.88, 0.58, 4,5,0,3, 100, 1
#20, 0.92, 0.58, 4,5,0,3, 100, 1
#34, 0.92, 0.58, 4,5,0,3, 50, 1
#36, 0.93, 0.75,  6,3,0,3, 50, 1
#44, 0.90, 0.66, 5,4,0,3, 100, 1
#47, 0.92, 0.75, 6,3,0, 3, 500, 1

#f2 = kbest.findkbest(df2train, df2test, feats=100)
#k2train, k2test = kbest.ktodf(f2, df2train, df2test)
#ada.gridne(k2train, k2test)
#ada.gridlr(k2train, k2test, ne=500)
#ada.predict(k2train, k2test, lr=0.01, ne=500)
#6, 0.87 0.58, 7,2,3,0, 100,0.1
#37, 0.80, 0.5, 6,3,3,0, 1000, 1
#38, 0.79, 0.58, 7,2,3,0, 500, 1
#48, 0.81, 0.58, 7,2,3,0, 1000,1
#52, 0.82, 0.5, 6,3,3,0, 500, 1
#53, 0.76, 0.5, 6,3,3,0, 1000,1
#74, 0.77, 0.5, 9,0,3,0, 1000,1
#100, 0.79, 0.75, 9,0,3,0, 500,0.01

# fage = kbest.findkbest(df10, df2, feats=54)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# ada.gridne(kagetrain, kagetest) 
# ada.gridlr(kagetrain, kagetest, ne=500)
# ada.predict(kagetrain, kagetest, lr=1, ne=500)
#15, 0.8, 0.57, 29,11,13,3, 500,0.1
#17, 0.8, 0.55, 27,13,11,5, 1000,1
#18, 0.83, 0.60, 29,11,11,5, 100,1 
#24, 0.78, 0.55, 27,13,12,4, 500,1
#41, 0.78, 0.6, 30,10,12,4, 500,1
#54, 0.75, 0.62, 32,8,13,3, 500,1
###############################################################################
#For GaussianNegativeBinomial get best params and training and testing results
#f10 = kbest.findkbest(df10train, df10test, feats=36)
#k10train, k10test = kbest.ktodf(f10, df10train, df10test)
#gnb.predict(k10train, k10test)
#8, 0.74, 0.5, 4,5,1,2, no smoother
#20, 0.83, 0.58, 4,5,0,3, no smoother
#34, 0.90, 0.83, 8,1,1,2, no smoother
#36, 0.93, 0.75,  6,3,0,3, 50, 1
#44, 0.92, 0.75, 6,3,0,3, no smoother
#47, 0.93, 0.75, 6,3,0,3, no smoother

#f2 = kbest.findkbest(df2train, df2test, feats=100)
#k2train, k2test = kbest.ktodf(f2, df2train, df2test)
#gnb.predict(k2train, k2test)
#6, 0.88, 0.5, 6,3,3,0, no smoother
#37, 0.81, 0.58, 7,2,3,0, no smoother
#38, 0.79, 0.58, 7,2,3,0, no smoother
#48, 0.81, 0.58, 7,2,3,0, no smoother
#52, 0.79, 0.58, 7,2,3,0, no smoother
#53, 0.79, 0.58, 7,2,3,0, no smoother
#74, 0.79, 0.58, 7,2,3,0, no smoother
#100, 0.82, 0.58, 7,2,3,0, no smoother

# fage = kbest.findkbest(df10, df2, feats=18)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# gnb.predict(kagetrain, kagetest) 
#15, 0.8, 0.57, 29,11,13,3, smoother
#17, 0.8, 0.55, 27,13,11,5,no smoother
#18, 0.81, 0.625, 31,9,12,4, no smoother
#24, 0.78, 0.55, 27,13,12,4, no smoother
#41, 0.78, 0.6, 30,10,12,4, smoother
#54, 0.78, 0.6, 30,10,12,4, smoother
###############################################################################
#For ExtraTrees get best params and training and testing results
# f10 = kbest.findkbest(df10train, df10test, feats=36)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
#extra.numtrees(k10train, k10test) 
#extra.numsplits(k10train, k10test, ne=50) 
# extra.predict(k10train, k10test, ne=50, ms=6)
#8, 0.88, 0.5, 4,5,1,2, 10,5
#20, 0.92, 0.58, 5,4,1,2, 10,5
#34, 0.87, 0.66, 6,4,1,2, 50,2
#36  0.93, 0.75, 6,3,0,3, 50,6
#44, 0.87, 0.58, 6,3,2,1, 100,4
#47, 0.85, 0.75, 8,1,2,1, 50,2

# f2 = kbest.findkbest(df2train, df2test, feats=100)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
#extra.numtrees(k2train, k2test) 
#extra.numsplits(k2train, k2test, ne=100)
# extra.predict(k2train, k2test, ne=100, ms=5)
#6, 0.85, 0.5, 7,2,3,0 50,11
#37, 0.87, 0.66, 8,1,3,0, 50,3
#38, 0.84, 0.66, 8,1,3,0, 100,2
#48, 0.84, 0.58, 7,2,3,0, 100,3
#52, 0.84, 0.66, 8,1,3,0, 10,4
#53, 0.85, 0.66, 8,1,3,0, 100,4
#74, 0.80, 0.58, 7,2,3,0, 100,5
#100, 0.79, 0.66, 8,1,3,0, 100,5

# fage = kbest.findkbest(df10, df2, feats=41)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# extra.numtrees(kagetrain, kagetest) 
# extra.numsplits(kagetrain, kagetest, ne=50) 
# extra.predict(kagetrain, kagetest, ne=100, ms=11)
#15, 0.85, 0.62, 33,7,14,2, 100,6
#17, 0.78, 0.60, 32,8,14,2, 100,9
#18, 0.81, 0.625, 33,7,13,2, 50,10
#24, 0.78, 0.58, 32,8,15,2, 100,9
#41, 0.78, 0.69, 37,3,14,2, 100,11
#54, 0.81, 0.69, 38,2,15,1, 100,9
###############################################################################
#For RandomForest get best params and training and testing results
# f10 = kbest.findkbest(df10train, df10test, feats=20)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
#rf.paramfind(k10train, k10test)
# rf.predict(k10train, k10test, ne=50, ml=1, ms=2, md=5)
#8, 0.88, 0.83, 7,2,0,3, 500
#20, 0.93, 0.75, 6,3,0,3, 50
#34, 0.88, 0.66, 6,3,1,2, 10
#36, 0.92, 0.66, 6,3,1,2, 50
#44, 0.92, 0.66, 6,3,1,2, 500
#47, 0.92, 0.83, 8,1,2,1, 50

# f2 = kbest.findkbest(df2train, df2test, feats=6)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
# rf.paramfind(k2train, k2test)
# rf.predict(k2train, k2test, ne=50, ml=1, ms=2, md=5)
# #6, 0.85, 0.66, 8,1,3,0, 50
#37, 0.84, 0.58, 7,2,3,0, 50
#38, 0.80, 0.41, 4,5,2,1, 10
#48, 0.79, 0.5, 6,3,3,0, 10
#52, 0.77, 0.58, 7,2,3,0, 10
#53, 0.81, 0.66, 8,1,3,0, 10
#74, 0.84, 0.66, 7,2,2,1, 10
#100, 0.80, 0.58, 7,2,3,0, 10

# fage = kbest.findkbest(df10, df2, feats=54)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# rf.paramfind(kagetrain, kagetest)
# rf.predict(kagetrain, kagetest, ne=10, ml=1, ms=2, md=5)
#15, 0.8, 0.57, 30,10,14,2, 10
#17, 0.78, 0.55, 29,11,14,2, 10
#18, 0.81, 0.58, 30,10,13,3, 50
#24, 0.82, 0.5, 26,14,14,2, 50
#41, 0.78, 0.64, 30,10,12,4, 500
#54, 0.73, 0.66, 35,5,14,2 10
###############################################################################
#best cases
#10m
# f10 = kbest.findkbest(df10train, df10test, feats=20)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# m, X_train, y_train, X_test, y_test, y_pred = rf.predict2(k10train, k10test,
#                                                   ne=50, ml=1, ms=2, md=5)
# plotbest.plotROC(m, X_test, y_test, "miRNA Aged Samples")

#2m
# f2 = kbest.findkbest(df2train, df2test, feats=6)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
# m, X_train, y_train, X_test, y_test, y_pred =  rf.predict2(k2train, k2test,
#                                                    ne=50, ml=1, ms=2, md=5)
# plotbest.plotROC(m, X_test, y_test, "miRNA Young Samples")
#age
# fage = kbest.findkbest(df10, df2, feats=41)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# m, X_train, y_train, X_test, y_test, y_pred =  extra.predict2(kagetrain, kagetest,
#                                                       ne=100, ms=11)
# plotbest.plotROC(m, X_test, y_test, "miRNA Predisposed Samples")
# ###############################################################################
#make data for R
# f10 = kbest.findkbest(df10train, df10test, feats=20)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# m, X_train, y_train, X_test, y_test, y_pred = rf.predict2(k10train, k10test,
#                                                   ne=50, ml=1, ms=2, md=5)
# plotbest.forR(X_test, y_test, y_pred, "miRNA_Aged.csv", f10)

# f2 = kbest.findkbest(df2train, df2test, feats=6)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
# m, X_train, y_train, X_test, y_test, y_pred =  rf.predict2(k2train, k2test,
#                                                     ne=50, ml=1, ms=2, md=5)
# plotbest.forR(X_test, y_test, y_pred, "miRNA_Young.csv", f2)

# fage = kbest.findkbest(df10, df2, feats=41)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# m, X_train, y_train, X_test, y_test, y_pred =  extra.predict2(kagetrain, kagetest,
#                                                       ne=100, ms=11)
# plotbest.forR(X_test, y_test, y_pred, "miRNA_Predisposed.csv", fage)