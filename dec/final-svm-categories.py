import json
import os
import pickle
import sys

import numpy as np
import pandas as pd
from scipy.io import loadmat
from tqdm.notebook import tqdm

with open('../dirs.json', 'r') as f:
    dirs = json.load(f)
sys.path.append(dirs['root'])

from utils.metrics import DescriminationConfidenceEstimator as DCE


def read_neuron(path, random_state=0):
    cm = loadmat(os.path.join(os.path.dirname(path), 'cm.mat'))['cm']
    nstim = np.unique(cm).max()
    
    nsample = 900
    data = np.nan * np.zeros([nstim, 1, nsample])
    
    ua = loadmat(path)['ua']
    cm = cm[:ua.shape[0]]
    
    exit_flag = 0
    split_1, split_2 = [], []
    for istim in np.arange(1, 166):
        np.random.seed(165 * random_state + istim)
        X = ua[(cm==istim).flatten(), :]
        nrep = X.shape[0]
        
        if nrep < 3:
            exit_flag = 1
            break

        pIndex = np.random.permutation(nrep)
        X = X[pIndex]

        first_half = X[:int(nrep/2), :]
        second_half = X[int(nrep/2):, :]

        split_1.append(first_half.mean(0))
        split_2.append(second_half.mean(0))
    
    return exit_flag, np.array(split_1), np.array(split_2)

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
outPath  = os.path.join(_outPath, f'{monkey.lower()}-{selectivity.lower()}')
os.makedirs(outPath, exist_ok=True)


for general_seed in (pbar := tqdm(range(150))):
    pbar.set_description(f"Repetition: {general_seed}")

    # Load the data
    itc_t, itc_v = [], []
    for monkey in ['jenab', 'zebel']:
        info = pd.read_csv(f'/Data/{selectivity.capitalize()}/{monkey.capitalize()}/itcNeuralInfo.csv')
        for row in info.iterrows():
            ef, spl1, spl2 = read_neuron(row[1][0], random_state=general_seed)
            if ef == 0:
                itc_t.append(spl1)
                itc_v.append(spl2)

    itc_t = np.transpose(itc_t, [1, 0, 2])
    itc_v = np.transpose(itc_v, [1, 0, 2])

    pfc_t, pfc_v = [], []
    for monkey in ['jenab', 'zebel']:
        info = pd.read_csv(f'/Data/{selectivity.capitalize()}/{monkey.capitalize()}/pfcNeuralInfo.csv')
        for row in info.iterrows():
            ef, spl1, spl2 = read_neuron(row[1][0], random_state=general_seed)
            if ef == 0:
                pfc_t.append(spl1)
                pfc_v.append(spl2)

    pfc_t = np.transpose(pfc_t, [1, 0, 2])
    pfc_v = np.transpose(pfc_v, [1, 0, 2])

    _data = {
        'itc_t': itc_t,
        'itc_v': itc_v,
        'pfc_t': pfc_t,
        'pfc_v': pfc_v
    }

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
            for seed in (range(nmb_rep)):
                with open("../utils/info.pkl", "rb") as handler:
                    info = pickle.load(handler)
                info = info[:165]
                validStimuliIndex = ((info.cat != 'none') & ((info.sfr == "A") | (info.sfr == "BI"))).to_numpy()
                data_t = movavg(_data[f'{region}_t'][validStimuliIndex], 25, 1)
                data_v = movavg(_data[f'{region}_v'][validStimuliIndex], 25, 1)
                info = info[validStimuliIndex].reset_index(drop=True)
            
                _cfn, _dpr, _dth, _d0v, _d1v = [], [], [], [], []
                
                y = info[category].to_numpy()
                y0 = np.argwhere(y).flatten()
                y1 = np.argwhere(~y)
                np.random.seed(seed=seed)
                y1 = np.random.choice(y1.flatten(), size=y0.size, replace=False).flatten()
                y = np.concatenate([y0, y1])
                info = info.loc[y, :]
                data_t, data_v = data_t[y], data_v[y]
                y = info[category].to_numpy()
                
                for itime in np.arange(data_t.shape[2]):
                    out = DCE(random_state=seed).fit(X_train=data_t[:, :, itime], y_train=y).score(data_v[:, :, itime], y_test=y)
                    _cfn.append(out[0])
                    _dth.append(out[1])
                    _dpr.append(out[2])
                    _d0v.append(out[3])
                    _d1v.append(out[4])

                
                cfn[region][category].append(_cfn)
                dth[region][category].append(_dth)
                dpr[region][category].append(_dpr)
                d0v[region][category].append(_d0v)
                d1v[region][category].append(_d1v)
    
    # Save output of repetition
    for var, name in zip([cfn, dth, dpr, d0v, d1v], ['cfn', 'dth', 'dpr', 'd0v', 'd1v']):
        with open(os.path.join(dirs['out']['dec'], f'c-ovr-r-{general_seed}-{name}.pickle'), 'wb') as handle:
            pickle.dump(var, handle, protocol=pickle.HIGHEST_PROTOCOL)
