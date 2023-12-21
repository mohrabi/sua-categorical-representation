import json
import os
import pickle
import sys

import numpy as np
import pandas as pd
from scipy.io import loadmat
from tqdm import tqdm

from sklearn.model_selection import train_test_split as tts

with open('../dirs.json', 'r') as f:
    dirs = json.load(f)
sys.path.append(dirs['root'])

from utils.metrics import DescriminationConfidenceEstimator as DCE

def read_neuron(path):
    cm = loadmat(os.path.join(os.path.dirname(path), 'cm.mat'))['cm']
    nstim = 165
        
    ua = loadmat(path)['ua']
    cm = cm[:ua.shape[0]]
    
    exit_flag = 0
    data = []
    for istim in np.arange(1, nstim+1):
        X = ua[(cm==istim).flatten(), :]
        nrep = X.shape[0]
        
        if nrep < 3:
            exit_flag = 1
            break

        pIndex = np.random.permutation(nrep)
        X = X[pIndex]

        data.append(X.mean(0))
    
    return exit_flag, np.array(data)

def movavg(inp, bl, ss):
    nbin = int(np.floor((inp.shape[2] - bl) / ss))
    out = np.nan * np.zeros((inp.shape[0], inp.shape[1], nbin))
    for ibin in range(nbin):
        out[:, :, ibin] = np.sum(inp[:, :, int(ibin*ss):int(ibin*ss)+bl], 2) / bl * 1000
    return out

## Main 
nmb_rep = 30
selectivity = 'fast'
monkey = 'both'

# Create output directories
_outPath = os.path.join(dirs['out']['dec'], 'c-ovr-final')
outPath  = os.path.join(_outPath, 'time-time-stimulus-split')
os.makedirs(outPath, exist_ok=True)

# Load the data
print("Loading the Data")
itc = []
for monkey in ['jenab', 'zebel']:
    info = pd.read_csv(f'/Data/{selectivity.capitalize()}/{monkey.capitalize()}/itcNeuralInfo.csv')
    for row in info.iterrows():
        ef, ua = read_neuron(row[1][0])
        if ef == 0:
            itc.append(ua)
itc = np.transpose(itc, [1, 0, 2])
itc = movavg(itc, 25, 10)

pfc = []
for monkey in ['jenab', 'zebel']:
    info = pd.read_csv(f'/Data/{selectivity.capitalize()}/{monkey.capitalize()}/pfcNeuralInfo.csv')
    for row in info.iterrows():
        ef, ua = read_neuron(row[1][0])
        if ef == 0:
            pfc.append(ua)
pfc = np.transpose(pfc, [1, 0, 2])
pfc = movavg(pfc, 25, 10)

_data = {
    'itc': itc,
    'pfc': pfc
}

print("Data Loaded Successfully")

for general_seed in range(2, 100):
    print(f"iteration: {general_seed}")

    # Initialize result dictionaries
    cfn = {'itc': {'fac': [], 'bod': [], 'nat': [], 'art': []},
        'pfc': {'fac': [], 'bod': [], 'nat': [], 'art': []}}
    dpr = {'itc': {'fac': [], 'bod': [], 'nat': [], 'art': []},
        'pfc': {'fac': [], 'bod': [], 'nat': [], 'art': []}}
    dth = {'itc': {'fac': [], 'bod': [], 'nat': [], 'art': []},
        'pfc': {'fac': [], 'bod': [], 'nat': [], 'art': []}}
    d0v = {'itc': {'fac': [], 'bod': [], 'nat': [], 'art': []},
        'pfc': {'fac': [], 'bod': [], 'nat': [], 'art': []}}
    d1v = {'itc': {'fac': [], 'bod': [], 'nat': [], 'art': []},
        'pfc': {'fac': [], 'bod': [], 'nat': [], 'art': []}}

    # Model training
    for region in ['itc', 'pfc']:
        for category in ['fac', 'bod', 'nat', 'art']:
            seed = general_seed
            
            print(f"{region} - {category}")
            
            with open("../utils/info.pkl", "rb") as handler:
                info = pickle.load(handler)
            info = info[:165]
            validStimuliIndex = ((info.cat != 'none') & ((info.sfr == "A") | (info.sfr == "BI"))).to_numpy()
            X = _data[region][validStimuliIndex]
            info = info[validStimuliIndex].reset_index(drop=True)

            assert(X.shape[0] == len(info))

            y = info[category].to_numpy()
            y0 = np.argwhere(y).flatten()
            y1 = np.argwhere(~y)
            np.random.seed(seed=seed)
            y1 = np.random.choice(y1.flatten(), size=y0.size, replace=False).flatten()
            y = np.concatenate([y0, y1])
            info = info.loc[y, :]
            X = X[y]
            y = info[category].to_numpy()

            np.random.seed(seed=seed)
            y = y[np.random.permutation(len(y)).flatten()]

            ind_tr, ind_te = tts(np.arange(len(y)), train_size=.5, random_state=seed, shuffle=True, stratify=y)

            X_train, X_test, y_train, y_test = X[ind_tr], X[ind_te], y[ind_tr], y[ind_te]

            _cfn, _dpr, _dth, _d0v, _d1v = [], [], [], [], []
            for itime in np.arange(X_train.shape[2]):
                __cfn, __dpr, __dth, __d0v, __d1v = [], [], [], [], []
                for jtime in np.arange(X_train.shape[2]):
                    mdl = DCE(random_state=seed).fit(X_train=X_train[:, :, itime], y_train=y_train)
                    out = mdl.score(X_test[:, :, jtime], y_test=y_test)
                    __cfn.append(out[0])
                    __dth.append(out[1])
                    __dpr.append(out[2])
                    __d0v.append(out[3])
                    __d1v.append(out[4])
                _cfn.append(__cfn)
                _dth.append(__dth)
                _dpr.append(__dpr)
                _d0v.append(__d0v)
                _d1v.append(__d1v)

            cfn[region][category].append(_cfn)
            dth[region][category].append(_dth)
            dpr[region][category].append(_dpr)
            d0v[region][category].append(_d0v)
            d1v[region][category].append(_d1v)

    # Save output of repetition
    for var, name in zip([cfn, dth, dpr, d0v, d1v], ['cfn', 'dth', 'dpr', 'd0v', 'd1v']):
        with open(os.path.join(outPath, f'c-ovr-p-{general_seed}-{name}.pickle'), 'wb') as handle:
            pickle.dump(var, handle, protocol=pickle.HIGHEST_PROTOCOL)