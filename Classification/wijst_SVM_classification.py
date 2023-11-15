import os
import pandas as pd
import numpy as np
from sklearn.svm import LinearSVC
from sklearn.model_selection import KFold
import time as tm

class Wijst_SVM:

    def __init__(self, data, labels, outputdir):
        self.data = data
        self.labels = labels
        self.outputdir = outputdir

    def run_SVM_crossval(self):
        os.chdir(self.outputdir)

        # normalize data
        data = np.log1p(self.data)

        Classifier = LinearSVC()

        tr_time=[]
        ts_time=[]
        truelab = []
        pred = []

        nfold = KFold(n_splits=5)
        for train_index, test_index in nfold.split(data):
            train_data, test_data = data.iloc[train_index], data.iloc[test_index]
            train_labels, test_labels = self.labels.iloc[train_index], self.labels.iloc[test_index]

            start=tm.time()
            Classifier.fit(train_data, train_labels.values.ravel())
            tr_time.append(tm.time()-start)

            start=tm.time()
            predicted = Classifier.predict(test_data)
            ts_time.append(tm.time()-start)

            truelab.extend(test_labels.values)
            pred.extend(predicted)


        truelab = pd.DataFrame(truelab)
        pred = pd.DataFrame(pred)

        tr_time = pd.DataFrame(tr_time)
        ts_time = pd.DataFrame(ts_time)

        truelab.to_csv("crossval_SVM_True_Labels.csv", index = False)
        pred.to_csv("crossval_SVM_Pred_Labels.csv", index = False)
        tr_time.to_csv("crossval_SVM_Training_Time.csv", index = False)
        ts_time.to_csv("crossval_SVM_Testing_Time.csv", index = False)

    def run_SVM_trainset(self, trainset, trainlabels):
        os.chdir(self.outputdir)

        # normalize data
        testdata = np.log1p(self.data)
        trainset = np.log1p(trainset)

        Classifier = LinearSVC()

        tr_time=[]
        ts_time=[]
        truelab = []
        pred = []

        start=tm.time()
        Classifier.fit(trainset, trainlabels.values.ravel())
        tr_time.append(tm.time()-start)

        start=tm.time()
        predicted = Classifier.predict(testdata)
        ts_time.append(tm.time()-start)

        truelab.extend(self.labels.values)
        pred.extend(predicted)

        truelab = pd.DataFrame(truelab)
        pred = pd.DataFrame(pred)

        tr_time = pd.DataFrame(tr_time)
        ts_time = pd.DataFrame(ts_time)

        truelab.to_csv("trained_SVM_True_Labels.csv", index = False)
        pred.to_csv("trained_SVM_Pred_Labels.csv", index = False)
        tr_time.to_csv("trained_SVM_Training_Time.csv", index = False)
        ts_time.to_csv("trained_SVM_Testing_Time.csv", index = False)
