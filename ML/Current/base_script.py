import os
import site
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
import seaborn as sns
print(sklearn.__version__)

# Add path to sources root to Python's PATH variable
site.addsitedir(os.path.dirname(os.path.dirname(os.path.abspath(''))))
from ML import *

if Import:
    Data = pd.read_csv(Data)
    Data = Data.drop(Data.columns[[0]], axis=1)
    print(Data.shape)

if ImportVal:
    Val = pd.read_csv(ValData)
    Val = Val.drop(Val.columns[[0]], axis=1)
    print(Val.shape)

if Hmap:
    from copy import deepcopy

    Data_corr = deepcopy(Data)
    Data_corr.Samples.replace(('HD', 'WT'), (1, 0), inplace=True)
    cor = Data_corr.corr(method='spearman')

    # f, ax = plt.subplots(figsize=(18, 18))
    # cmap = sns.diverging_palette(220, 10, as_cmap=True)
    # heatmap = sns.heatmap(cor, cmap=cmap, center=0.0, vmax=1,
    #                       linewidths=1, ax=ax)
    # plt.savefig(os.path.join(PLOT_DIR, 'hmap_train.png'))

if FE:
    
    def correlation(dataset, threshold):
        col_corr = set()  # Set of all the names of correlated columns
        corr_matrix = dataset.corr()
        for i in range(len(corr_matrix.columns)):
            for j in range(i):
                if abs(corr_matrix.iloc[i, j]) > threshold:  # we are interested in absolute coeff value
                    colname = corr_matrix.columns[i]  # getting the name of column
                    col_corr.add(colname)
        return col_corr

    corr_features = correlation(Data_corr, 0.8)
    len(set(corr_features))
    print(corr_features)

    Data_FE = Data.drop(corr_features, axis=1)
    Val_FE = Val.drop(corr_features, axis=1)

    # print(Data_FE)
    # print(Val_FE.shape)

if KFold:
    from sklearn import datasets, metrics, model_selection, svm
    from sklearn.neural_network import MLPClassifier
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.svm import SVC
    from sklearn.gaussian_process import GaussianProcessClassifier
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.gaussian_process.kernels import RBF
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import ExtraTreesClassifier
    from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
    from sklearn.naive_bayes import GaussianNB
    from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
    from sklearn.linear_model import SGDClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from imblearn.over_sampling import RandomOverSampler
    from sklearn import preprocessing
    from sklearn.preprocessing import MinMaxScaler

    names = ["Nearest_Neighbors", "Linear_SVM", "Polynomial_SVM", "RBF_SVM", "Gaussian_Process",
             "Gradient_Boosting", "Decision_Tree", "Extra_Trees", "Random_Forest", "Neural_Net", "AdaBoost",
             "Naive_Bayes", "QDA", "SGD", "LogisticRegression"]

    classifiers = [
        KNeighborsClassifier(8),
        SVC(kernel="linear", C=0.025),
        SVC(kernel="poly", degree=3, C=0.025),
        SVC(kernel="rbf", C=1, gamma=2),
        GaussianProcessClassifier(1.0 * RBF(1.0)),
        GradientBoostingClassifier(n_estimators=100, learning_rate=1),
        DecisionTreeClassifier(max_depth=5),
        ExtraTreesClassifier(n_estimators=10, min_samples_split=2),
        RandomForestClassifier(max_depth=5, n_estimators=100),
        MLPClassifier(alpha=1, max_iter=1000),
        AdaBoostClassifier(n_estimators=100),
        GaussianNB(),
        QuadraticDiscriminantAnalysis(),
        SGDClassifier(loss="hinge", penalty="l2"),
        LogisticRegression(solver='liblinear', multi_class='ovr')]

    TrainDF = []
    TestDF = []

    cv = StratifiedKFold(n_splits=5, random_state=42, shuffle=True)

    array = Data_FE.values
    X_train = array[:, 0:26]
    y_train = array[:, 27]

    array2 = Val_FE.values
    X_test = array2[:, 0:26]
    y_test = array2[:, 27]

    print(y_train)
    print(y_test)

    for train_idx, test_idx, in cv.split(X_train, y_train):

        ROS = RandomOverSampler(random_state=42)
        X_train_oversampled, y_train_oversampled = ROS.fit_resample(X_train, y_train)
        X_test_oversampled, y_test_oversampled = ROS.fit_resample(X_test, y_test)

        scaler = min_max_scaler = preprocessing.MinMaxScaler().fit(X_train_oversampled)
        X_train_scaled = scaler.transform(X_train_oversampled)
        X_test_scaled = scaler.transform(X_test_oversampled)

        Train = []
        Test = []
        for name, clf in zip(names, classifiers):
            clf.fit(X_train_scaled, y_train_oversampled)
            train_score = clf.score(X_train_scaled, y_train_oversampled)
            test_score = clf.score(X_test_scaled, y_test_oversampled)

            Train.append(train_score)
            Test.append(test_score)

        TrainDF.append(Train)
        TestDF.append(Test)

if ClassPlots:
    Training = pd.DataFrame(TrainDF, columns=names)
    Train_mean = Training.mean(axis=0)
    df = pd.DataFrame({'Classifiers':Train_mean.index, 'Prediction Scores':Train_mean.values})
    print(df)
    f, ax = plt.subplots(figsize=(12, 8))
    sns.set(style="whitegrid")
    ax = sns.barplot(y="Classifiers", x="Prediction Scores", data=df).set_title("Classfier performance on Training data")
    plt.savefig(os.path.join(PLOT_DIR, 'ClassCompare_Train.png'))
    plt.close()

    Testing = pd.DataFrame(TestDF, columns=names)
    Test_mean = Testing.mean(axis=0)
    df2 = pd.DataFrame({'Classifiers': Test_mean.index, 'Prediction Scores': Test_mean.values})
    print(df2)
    f, ax = plt.subplots(figsize=(12, 8))
    sns.set(style="whitegrid")
    ax = sns.barplot(y="Classifiers", x="Prediction Scores", data=df2).set_title("Classfier performance on Training data")
    plt.savefig(os.path.join(PLOT_DIR, 'ClassCompare_Test.png'))
    plt.close()

if LogisticRegression:
    from sklearn.model_selection import StratifiedKFold
    from imblearn.over_sampling import RandomOverSampler
    from sklearn import preprocessing
    from sklearn.preprocessing import MinMaxScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn import metrics
    from sklearn.metrics import confusion_matrix

    cv = StratifiedKFold(n_splits=5, random_state=42, shuffle=True)

    Acc = []
    ROC = []
    CM = []

    array = Data_FE.values
    X_train = array[:, 0:26]
    y_train = array[:, 27]

    array2 = Val_FE.values
    X_test = array2[:, 0:26]
    y_test = array2[:, 27]
    for train_idx, test_idx, in cv.split(X_train, y_train):
        ROS = RandomOverSampler(random_state=42)
        X_train_oversampled, y_train_oversampled = ROS.fit_resample(X_train, y_train)
        X_test_oversampled, y_test_oversampled = ROS.fit_resample(X_test, y_test)

        scaler = min_max_scaler = preprocessing.MinMaxScaler().fit(X_train_oversampled)
        X_train_scaled = scaler.transform(X_train_oversampled)
        X_test_scaled = scaler.transform(X_test_oversampled)

        model = LogisticRegression()
        model.fit(X_train_scaled, y_train_oversampled)

        y_pred = model.predict(X_train_scaled)
        Acc_score = model.score(X_test_scaled, y_test_oversampled)
        roc_score = roc_auc_score(y_test_oversampled, model.decision_function(X_test_scaled))
        confmat = metrics.confusion_matrix(y_test_oversampled, model.predict(X_test_scaled))

        Acc.append(Acc_score)
        ROC.append(roc_score)
        CM.append(confmat)

    # print(f'Validation Accuracy: {Acc}')
    # print(f'ROC scores: {ROC}')
    # print(f'Confusion matrices: {CM}')

if ConfusionMatrix:
    from sklearn.metrics import plot_confusion_matrix
    # Plot normalized confusion matrix
    titles_options = [("Confusion matrix", None),
                      ("Normalised confusion matrix", 'true')]
    for title, normalize in titles_options:
        disp = plot_confusion_matrix(model, X_test_scaled, y_test_oversampled,
                                     display_labels=['HD', 'WT'],
                                     cmap=plt.cm.Blues,
                                     normalize='true')
        disp.ax_.set_title(title)

        # print(title)
        # print(disp.confusion_matrix)
        plt.savefig(os.path.join(PLOT_DIR, 'ConfuMat.png'))

if ROC:
    import sklearn.metrics as metrics
    y_pred = model.predict(X_test_scaled)
    y_test = np.where(y_test_oversampled =="HD", 1, 0)
    y_proba = np.where(y_pred =="HD", 1, 0)
    fpr, tpr, _ = metrics.roc_curve(y_test, y_proba)
    auc = metrics.roc_auc_score(y_test, y_proba)
    plt.plot(fpr, tpr, label="HD sample score, auc=" + str(auc))
    plt.legend(loc=4)
    plt.plot([0, 1], [0, 1], 'r--')
    plt.ylabel('True Positive Rate (Positive label = HD)')
    plt.xlabel('False Positive Rate (Negative label = WT)')
    plt.savefig(os.path.join(PLOT_DIR, 'ROCAUC.png'))

