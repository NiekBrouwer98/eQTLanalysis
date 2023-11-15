from run_SVM_rejection import run_SVM
import os
import pandas as pd
import numpy as np
from sklearn.svm import LinearSVC
from sklearn.model_selection import KFold
from wijst_data_preprocessing import get_dataframe, get_labels, get_barcodes_to_labels
from sklearn.calibration import CalibratedClassifierCV
import time as tm

class Wijst_SVMrejection:

    def __init__(self, data, labels, outputdir):
        self.data = data
        self.labels = labels
        self.outputdir = outputdir

    def run_SVM_rejection_crossval(self):
        os.chdir(self.outputdir)

        # normalize data
        data = np.log1p(self.data)

        Classifier = LinearSVC()
        clf = CalibratedClassifierCV(Classifier)

        tr_time=[]
        ts_time=[]
        truelab = []
        pred = []

        nfold = KFold(n_splits=5)
        for train_index, test_index in nfold.split(data):
            train_data, test_data = data.iloc[train_index], data.iloc[test_index]
            train_labels, test_labels = self.labels.iloc[train_index], self.labels.iloc[test_index]

            start=tm.time()
            clf.fit(train_data, train_labels.values.ravel())
            tr_time.append(tm.time()-start)

            start=tm.time()
            predicted = clf.predict(test_data)
            prob = np.max(clf.predict_proba(test_data), axis = 1)
            unlabeled = np.where(prob < 0.7)
            predicted[unlabeled] = 'Unknown'
            ts_time.append(tm.time()-start)

            truelab.extend(test_labels.values)
            pred.extend(predicted)

        truelab = pd.DataFrame(truelab)
        pred = pd.DataFrame(pred)

        tr_time = pd.DataFrame(tr_time)
        ts_time = pd.DataFrame(ts_time)

        truelab.to_csv("crossval_SVMrej_True_Labels.csv", index = False)
        pred.to_csv("crossval_SVMrej_Pred_Labels.csv", index = False)
        tr_time.to_csv("crossval_SVMrej_Training_Time.csv", index = False)
        ts_time.to_csv("crossval_SVMrej_Testing_Time.csv", index = False)

    def run_SVM_rejection_trainset(self, trainset, trainlabels):
        os.chdir(self.outputdir)

        # normalize data
        testdata = np.log1p(self.data)
        trainset = np.log1p(trainset)

        Classifier = LinearSVC()
        clf = CalibratedClassifierCV(Classifier)

        tr_time=[]
        ts_time=[]
        truelab = []
        pred = []

        start=tm.time()
        clf.fit(trainset, trainlabels.values.ravel())
        tr_time.append(tm.time()-start)

        start=tm.time()
        predicted = clf.predict(testdata)
        prob = np.max(clf.predict_proba(testdata), axis = 1)
        unlabeled = np.where(prob < 0.7)
        predicted[unlabeled] = 'Unknown'
        ts_time.append(tm.time()-start)

        truelab.extend(self.labels.values)
        pred.extend(predicted)

        truelab = pd.DataFrame(truelab)
        pred = pd.DataFrame(pred)

        tr_time = pd.DataFrame(tr_time)
        ts_time = pd.DataFrame(ts_time)

        truelab.to_csv("trained_SVMrej_True_Labels.csv", index = False)
        pred.to_csv("trained_SVMrej_Pred_Labels.csv", index = False)
        tr_time.to_csv("trained_SVMrej_Training_Time.csv", index = False)
        ts_time.to_csv("trained_SVMrej_Testing_Time.csv", index = False)
