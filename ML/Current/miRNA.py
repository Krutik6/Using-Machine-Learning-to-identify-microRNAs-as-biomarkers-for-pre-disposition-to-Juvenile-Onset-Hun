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
#274 features found

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
# trainclasses2, testclasses2 = multiclasstest.mulitanalysis2(rfe2train,rf2etest, 5)
#multiclasstest.makeconfusions(rfe2train,rf2etest)
#ADA, NB and Extratees are best
################################################################################
#recursive feature elimination 
#use old data to predict young data - better than kbest
##Feature Selection
# supp = recfs.reccySVC(df10, df2, minfeat=5)
# rfeagetrain,rfeagetest = recfs.retreiveFeat(df10, df2, supp)
#127 features found

##Train Multiple Classes
# trainclassesage, testclassesage = multiclasstest.mulitanalysis2(rfeagetrain,rfeagetest, 5)
#multiclasstest.makeconfusions(rfeagetrain,rfeagetest)
#ADA, Extratrees and RF are best
################################################################################
#kbest only using ADA, Extratrees, NGB and RF
###############################################################################
#try only with 10m data
##Feature Selection
# ada10, gnb10, et10, rf10 = kbest.findbestvalk(df10train, df10test, reps=100)
###############################################################################
#try only with 2m data
##Feature Selection
#ada2, gnb2, et2, rf2 = kbest.findbestvalk(df2train, df2test, reps=100)
################################################################################
#use old data to predict young data
##Feature Selection
# ada, gnb, et, rf = kbest.findbestvalk(df10, df2, reps=100)
###############################################################################
#Hypertuning parameters
#10 - 6, 11, 45, 71
#2 - 6, 37, 38, 52, 100  
#age - 11, 21, 25, 53, 80
###############################################################################
#For ADAboost get best params and training and testing results
# f10 = kbest.findkbest(df10train, df10test, feats=71)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# ada.gridne(k10train, k10test) 
# ada.gridlr(k10train, k10test, ne=1000) 
# ada.predict(k10train, k10test, lr=1, ne=1000)

#f2 = kbest.findkbest(df2train, df2test, feats=100)
#k2train, k2test = kbest.ktodf(f2, df2train, df2test)
#ada.gridne(k2train, k2test)
#ada.gridlr(k2train, k2test, ne=500)
#ada.predict(k2train, k2test, lr=0.01, ne=500)

# fage = kbest.findkbest(df10, df2, feats=80)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# ada.gridne(kagetrain, kagetest) 
# ada.gridlr(kagetrain, kagetest, ne=100)
# ada.predict(kagetrain, kagetest, lr=1, ne=100)
###############################################################################
#For GaussianNegativeBinomial get best params and training and testing results
# f10 = kbest.findkbest(df10train, df10test, feats=71)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# gnb.predict(k10train, k10test)

#f2 = kbest.findkbest(df2train, df2test, feats=100)
#k2train, k2test = kbest.ktodf(f2, df2train, df2test)
#gnb.predict(k2train, k2test)

# fage = kbest.findkbest(df10, df2, feats=80)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# gnb.predict(kagetrain, kagetest) 
###############################################################################
#For ExtraTrees get best params and training and testing results
# f10 = kbest.findkbest(df10train, df10test, feats=71)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# extra.numtrees(k10train, k10test) 
# extra.numsplits(k10train, k10test, ne=1000) 
# extra.predict(k10train, k10test, ne=1000, ms=10)

# f2 = kbest.findkbest(df2train, df2test, feats=100)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
#extra.numtrees(k2train, k2test) 
#extra.numsplits(k2train, k2test, ne=100)
# extra.predict(k2train, k2test, ne=100, ms=5)

# fage = kbest.findkbest(df10, df2, feats=80)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# extra.numtrees(kagetrain, kagetest) 
# extra.numsplits(kagetrain, kagetest, ne=500) 
# extra.predict(kagetrain, kagetest, ne=500, ms=4)
###############################################################################
#For RandomForest get best params and training and testing results
# f10 = kbest.findkbest(df10train, df10test, feats=71)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# rf.paramfind(k10train, k10test)
# rf.predict(k10train, k10test, ne=100, ml=1, ms=2, md=5)

# f2 = kbest.findkbest(df2train, df2test, feats=6)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
# rf.paramfind(k2train, k2test)
# rf.predict(k2train, k2test, ne=50, ml=1, ms=2, md=5)

# fage = kbest.findkbest(df10, df2, feats=80)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# rf.paramfind(kagetrain, kagetest)
# rf.predict(kagetrain, kagetest, ne=10, ml=1, ms=2, md=5)
###############################################################################
#best cases
#10m
# f10 = kbest.findkbest(df10train, df10test, feats=11)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# m, X_train, y_train, X_test, y_test, y_pred = ada.predict2(k10train, k10test,
#                                                   ne=50, lr=1)
# plotbest.plotROC(m, X_test, y_test, "miRNA Aged Samples")

#2m
# f2 = kbest.findkbest(df2train, df2test, feats=6)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
# m, X_train, y_train, X_test, y_test, y_pred =  rf.predict2(k2train, k2test,
#                                                     ne=50, ml=1, ms=2, md=5)
# plotbest.plotROC(m, X_test, y_test, "miRNA Young Samples")
#age
# fage = kbest.findkbest(df10, df2, feats=80)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# m, X_train, y_train, X_test, y_test, y_pred =  rf.predict2(kagetrain, kagetest,
#                                                       ne=10)
# plotbest.plotROC(m, X_test, y_test, "miRNA Predisposed Samples")
# ###############################################################################
#make data for R
# f10 = kbest.findkbest(df10train, df10test, feats=11)
# k10train, k10test = kbest.ktodf(f10, df10train, df10test)
# m, X_train, y_train, X_test, y_test, y_pred = ada.predict2(k10train, k10test,
#                                                   ne=50, lr=1)
# plotbest.forR(X_test, y_test, y_pred, "miRNA_Aged.csv", f10)

# f2 = kbest.findkbest(df2train, df2test, feats=6)
# k2train, k2test = kbest.ktodf(f2, df2train, df2test)
# m, X_train, y_train, X_test, y_test, y_pred =  rf.predict2(k2train, k2test,
#                                                     ne=50, ml=1, ms=2, md=5)
# plotbest.forR(X_test, y_test, y_pred, "miRNA_Young.csv", f2)

# fage = kbest.findkbest(df10, df2, feats=80)
# kagetrain, kagetest = kbest.ktodf(fage, df10, df2)
# m, X_train, y_train, X_test, y_test, y_pred =  rf.predict2(kagetrain, kagetest,
#                                                       ne=10)
# plotbest.forR(X_test, y_test, y_pred, "miRNA_Predisposed.csv", fage)