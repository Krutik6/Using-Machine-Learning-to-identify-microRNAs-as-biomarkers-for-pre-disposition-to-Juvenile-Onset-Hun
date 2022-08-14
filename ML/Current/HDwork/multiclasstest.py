# -*- coding: utf-8 -*-
import os, sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import structuring

from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
import pandas as pd
import plotly
import plotly.graph_objs as go
import plotly.offline as ply
plotly.offline.init_notebook_mode()
from sklearn.model_selection import cross_val_score
import numpy as np
import matplotlib.pyplot as plt 
from sklearn.metrics import plot_confusion_matrix

#global parameters
PLOT_DIR = 'C:\\Users\\nkp68\\OneDrive - Newcastle University\\Writing_Papers\\HD\\ML\\Current\\Info\\Plot'

def multiclasstestfunc():
  print("multiclasstest module loaded")
  
  
def mulitanalysis(train, test, C):
    names = ["Linear_SVM", "Polynomial_SVM",
              "Gaussian_Process",
              "Extra_Trees", "Random_Forest", "Neural_Net", "AdaBoost",
              "Naive_Bayes"]

    classifiers = [
        SVC(kernel="linear", C=0.025),
        SVC(kernel="poly", degree=3, C=0.025),
        GaussianProcessClassifier(1.0 * RBF(1.0)),
        ExtraTreesClassifier(n_estimators=10, min_samples_split=2),
        RandomForestClassifier(max_depth=10, n_estimators=10),
        MLPClassifier(alpha=1, max_iter=1000),
        AdaBoostClassifier(n_estimators=10),
        GaussianNB()]
    
    X_train, y_train, X_test, y_test = structuring.toArrays(train, test)
    
    TrainDF = []
    TestDF = []
    cv = StratifiedKFold(n_splits=C, random_state=42, shuffle=True)
    
    for train_idx, test_idx, in cv.split(X_train, y_train):
      X_train, y_train = structuring.smotey(X_train, y_train)
      X_train, X_test = structuring.scale1to0(X_train, X_test)
      
      TrainScores = []
      TestScores = []
      
      for name, clf in zip(names, classifiers):
         clf.fit(X_train, y_train)
         train_score = clf.score(X_train, y_train)
         test_score = clf.score(X_test, y_test)

         TrainScores.append(train_score)
         TestScores.append(test_score)

      TrainDF.append(TrainScores)
      TestDF.append(TestScores)
    
    TrainMultiClassScores = pd.DataFrame(TrainDF, columns=names)
    TrainMultiClassScores.loc['mean'] = TrainMultiClassScores.mean()

    TestMultiClassScores = pd.DataFrame(TestDF, columns=names)
    TestMultiClassScores.loc['mean'] = TestMultiClassScores.mean()
      
    return TrainMultiClassScores, TestMultiClassScores



def algorithmtest(traindf, testdf):
    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
    
    names = ["Linear_SVM", "Polynomial_SVM","Gaussian_Process","Extra_Trees",
             "Random_Forest", "Neural_Net", "AdaBoost","Naive_Bayes"]

    classifiers = [
      SVC(kernel="linear", C=0.025),
      SVC(kernel="poly", degree=3, C=0.025),
      GaussianProcessClassifier(1.0 * RBF(1.0)),
      ExtraTreesClassifier(n_estimators=10, min_samples_split=2),
      RandomForestClassifier(max_depth=10, n_estimators=100),
      MLPClassifier(alpha=1, max_iter=1000),
      AdaBoostClassifier(n_estimators=100),
      GaussianNB()]

    plot_data_train=[]

    clf_data=zip(names, classifiers)

    for clf_name, clf in clf_data:
        print('Running '+clf_name)
        kf=StratifiedKFold(n_splits=10, shuffle=True)
        train_scores=cross_val_score(clf, X_train, y_train, cv=kf)
        test_scores=cross_val_score(clf, X_test, y_test, cv=kf)
        print(train_scores)
        plot_data_train.append(
            go.Scatter(
                x=[i+1 for i in range(10)],
                y=train_scores,
                mode='lines',
                name=clf_name)
        )
        
    layout = go.Layout(
        xaxis=dict(
            title='Fold no.'
        ),
        yaxis=dict(
            range=[np.min([i['y'] for i in plot_data_train]), 1],
            title='Accuracy'
        )
    )

    fig=go.Figure(data=plot_data_train, layout=layout)
    ply.iplot(fig)
    fig.show(renderer="svg")
    fig.write_image(os.path.join(PLOT_DIR, 'algttrain.png'))
    
    
def makeconfusions(traindf, testdf):

    X_train, y_train, X_test, y_test = structuring.toArrays(traindf, testdf)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
    
    classifiers = [
        SVC(kernel="linear", C=0.025),
        SVC(kernel="poly", degree=3, C=0.025),
        GaussianProcessClassifier(1.0 * RBF(1.0)),
        ExtraTreesClassifier(n_estimators=10, min_samples_split=2),
        RandomForestClassifier(max_depth=10, n_estimators=100),
        MLPClassifier(alpha=1, max_iter=1000),
        AdaBoostClassifier(n_estimators=100),
        GaussianNB()]
    for cls in classifiers:
        cls.fit(X_train, y_train)
        
    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(15,10))
    
    for cls, ax in zip(classifiers, axes.flatten()):
        plot_confusion_matrix(cls, 
                              X_test, 
                              y_test, 
                              ax=ax, 
                              cmap='Blues',
                              display_labels=['HD', 'WT'])
        ax.title.set_text(type(cls).__name__)
    plt.tight_layout()  
    plt.savefig(os.path.join(PLOT_DIR, 'confusions.png')) 
    
    
    
def mulitanalysis2(train, test, C):
    names = ["Linear_SVM", "Polynomial_SVM",
              "Gaussian_Process",
              "Extra_Trees", "Random_Forest", "Neural_Net", "AdaBoost",
              "Naive_Bayes"]

    classifiers = [
        SVC(kernel="linear", C=0.025),
        SVC(kernel="poly", degree=3, C=0.025),
        GaussianProcessClassifier(1.0 * RBF(1.0)),
        ExtraTreesClassifier(n_estimators=10, min_samples_split=2),
        RandomForestClassifier(max_depth=10, n_estimators=10),
        MLPClassifier(alpha=1, max_iter=1000),
        AdaBoostClassifier(n_estimators=10),
        GaussianNB()]
    
    X_train, y_train, X_test, y_test = structuring.toArrays(train, test)
    X_train, y_train = structuring.smotey(X_train, y_train)
    X_train, X_test = structuring.scale1to0(X_train, X_test)
    
    TrainDF = []
    TestDF = []
    cv = StratifiedKFold(n_splits=C, random_state=42, shuffle=True)
    
    for train_idx, test_idx, in cv.split(X_train, y_train):
      
      TrainScores = []
      TestScores = []
      
      for name, clf in zip(names, classifiers):
         clf.fit(X_train, y_train)
         train_score = clf.score(X_train, y_train)
         test_score = clf.score(X_test, y_test)

         TrainScores.append(train_score)
         TestScores.append(test_score)

      TrainDF.append(TrainScores)
      TestDF.append(TestScores)
    
    TrainMultiClassScores = pd.DataFrame(TrainDF, columns=names)
    TrainMultiClassScores.loc['mean'] = TrainMultiClassScores.mean()

    TestMultiClassScores = pd.DataFrame(TestDF, columns=names)
    TestMultiClassScores.loc['mean'] = TestMultiClassScores.mean()
      
    return TrainMultiClassScores, TestMultiClassScores