import numpy as np
from numpy import array
import numpy_indexed as npi

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