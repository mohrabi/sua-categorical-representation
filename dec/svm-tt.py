"""
Initial code used to implement time-time classifiers on data
Key notes:
    * In this code, at each time-bin, a seperate pca transform is applied to data
Discussion:
    * In the IT cortex, the model trained on first information peaks has a 20 percent accuracy on the second peak
TO DO:
    * Omit the effects of pca, so that the results be are only affected by the amount of information
"""

import os

import numpy as np
from imblearn.under_sampling import RandomUnderSampler as rus
from sklearn.decomposition import PCA
from sklearn.metrics import confusion_matrix, roc_auc_score
from sklearn.model_selection import train_test_split as tts
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from tqdm import tqdm

from src.ipm import NeuronLoader, grablabels


def mkeqdata(data, groups, n, seed=0):
    np.random.seed(seed)
    ind0 = np.random.choice(np.argwhere(groups==0).squeeze(), size=n, replace=False)
    ind1 = np.random.choice(np.argwhere(groups==1).squeeze(), size=n, replace=False)
    
    X = np.concatenate([data[ind0], data[ind1]], axis=0)
    y = np.concatenate([np.zeros_like(ind0), np.ones_like(ind1)])
    return X, y

def initmetric(nmb_bin, nmb_rep):
    cfn = np.nan * np.ones([2, 2, nmb_bin, nmb_bin, nmb_rep])
    dth = np.nan * np.ones([nmb_bin, nmb_bin, nmb_rep])
    dpr = np.nan * np.ones([nmb_bin, nmb_bin, nmb_rep])
    auc = np.nan * np.ones([nmb_bin, nmb_bin, nmb_rep])
    return cfn, dth, dpr, auc

def calcmetric(mdl, X_test, y_test):
    y_pred = mdl.predict(X_test)
    y_prob = mdl.predict_proba(X_test)
    dist = mdl['svc'].decision_function(X_test) / np.linalg.norm(mdl['svc'].coef_)
    d0, d1 = dist[y_test==0], dist[y_test==1]
    
    cfn = confusion_matrix(y_test, y_pred, labels=[0, 1], sample_weight=None, normalize=None)
    dth = np.abs(d0.mean() - d1.mean())
    dpr = np.sqrt(2) * np.abs(d0.mean() - d1.mean()) / np.sqrt(d0.var() + d1.var())
    auc = roc_auc_score(y_test, y_prob[:, 1])
    
    return cfn, dth, dpr, auc

def savemetric(path: str, prefix: str, cfn: np.array, dth: np.array, dpr: np.array, auc: np.array, tim: np.array) -> None:
    np.save(os.path.join(path, prefix + '-cfn.npy'), cfn)
    np.save(os.path.join(path, prefix + '-dth.npy'), dth)
    np.save(os.path.join(path, prefix + '-dpr.npy'), dpr)
    np.save(os.path.join(path, prefix + '-auc.npy'), auc)
    np.save(os.path.join(path, prefix + '-tim.npy'), tim)

nmb_component = .95
nmb_rep = 50

for monkey in ["both", "jenab", "zebel"]:
    for selectivity in ["fast", "slow"]:
        data = NeuronLoader(f'G:/Data/{selectivity.capitalize()}/{monkey.capitalize()}')
        out_path  = f'G:\\Codes\\Processing\\out\\svm-tt\\{monkey.lower()}-{selectivity.lower()}'
        os.makedirs(out_path, exist_ok=True)

        time = np.arange(0, data.it.shape[2], 20)
        
        for region in ['it', 'pfc']:
            region_data = data.it if region=='it' else data.pfc

            for categ in ['super-ordinate', 'mid-level']:
                nmb_bin = time.size
                cfn, dth, dpr, auc = initmetric(nmb_bin, nmb_rep)
                
                yc = grablabels(hierarchy=categ)
                nmb_sample = np.min([(yc==0).sum(), (yc==1).sum()])

                for seed in (pbar := tqdm(range(nmb_rep))):
                    pbar.set_description(monkey+"-"+selectivity+"-"+region+"-"+categ)
                    X, y = mkeqdata(region_data, yc, nmb_sample, seed=seed)

                    for itr, ibin_train in enumerate(time):
                        # pca = PCA(n_components=None if nmb_component==-1 else nmb_component, 
                        #           copy=True, whiten=False, svd_solver='auto', random_state=seed)
                        # Xr = pca.fit_transform(X[:,:,ibin_train])
                        Xr = X[:,:,ibin_train]
                        X_train, X_test, y_train, y_test = tts(Xr, y, test_size=.3,
                                                               random_state=seed, shuffle=True, stratify=y)
                        mdl = make_pipeline(StandardScaler(), SVC(C=2.4, probability=True, kernel='linear', 
                                                                  random_state=seed)).fit(X_train, y_train)
                        
                        ind_train, ind_test = tts(np.arange(X.shape[0]), test_size=.3, 
                                                  random_state=seed, shuffle=True, stratify=y)
                        for ite, ibin_test in enumerate(time):
                            # Xr = pca.transform(X[:,:,ibin_test])
                            Xr = X[:,:,ibin_test]
                            X_train, X_test, y_train, y_test = Xr[ind_train], Xr[ind_test], y[ind_train], y[ind_test]
                            
                            cfn[:,:,itr,ite,seed], dth[itr,ite,seed], \
                            dpr[itr,ite,seed], auc[itr,ite,seed] = calcmetric(mdl, X_test, y_test)
                
                savemetric(out_path, region.lower() + "-" + categ, cfn, dth, dpr, auc, time-175)

