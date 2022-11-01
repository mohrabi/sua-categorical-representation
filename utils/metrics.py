from random import random
import numpy as np
from numpy import array
import numpy_indexed as npi

from sklearn.metrics import confusion_matrix
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.covariance import EmpiricalCovariance as EC

def sepind(trials, labels, estimator=EC, order=2):
    b, w, t, m = [], [], [], []

    nanind = np.where(np.bitwise_or(np.isnan(labels).squeeze(), np.isnan(trials[:,0,0])))
    trials = np.delete(trials, nanind, 0)
    y = np.delete(labels, nanind)
    for i_bin in range(trials.shape[2]):
        X = trials[:, :, i_bin]
        
        s_total = estimator().fit(X).covariance_ * X.shape[0]
        s_within = np.zeros_like(s_total)
        for group in npi.group_by(y).split(X):
            s_within += estimator().fit(group).covariance_ * group.shape[0]
        s_between = s_total - s_within

        m.append(np.linalg.norm(s_between, order) / np.linalg.norm(s_within, order))
        b.append(np.linalg.norm(s_between, order))
        w.append(np.linalg.norm(s_within, order))
        t.append(np.linalg.norm(s_total, order))
    
    stats = {
        "scatter_t": array(t),
        "scatter_b": array(b),
        "scatter_w": array(w)
    }
    return array(m), stats


class DescriminationConfidenceEstimator:
    def __init__(self, random_state=0):
        self.mdl = []
        self.random_state = random_state

    def fit(self, X_train, y_train):
        self.mdl = make_pipeline(StandardScaler(), SVC(C=2.4, probability=True, kernel='linear', 
            random_state=self.random_state)).fit(X_train, y_train)
        return self

    def score(self, X_test, y_test):
        y_pred = self.mdl.predict(X_test)
        dist = self.mdl['svc'].decision_function(X_test) / np.linalg.norm(self.mdl['svc'].coef_)
        d0, d1 = dist[y_test==0], dist[y_test==1]
        
        n   = X_test.shape[1]
        cfn = confusion_matrix(y_test, y_pred, labels=[0, 1], sample_weight=None, normalize=None)
        dth = np.abs(d0.mean() - d1.mean())
        dpr = np.sqrt(2 / n) * np.abs(d0.mean() - d1.mean()) / np.sqrt(d0.var() + d1.var())
        return cfn, dth, dpr, d0.var(), d1.var()