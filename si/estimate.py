
import os
import sys
import pickle
import json

with open('dirs.json', 'r') as f:
    dirs = json.load(f)
sys.path.append(dirs['root'])

import numpy as np
from sklearn.decomposition import PCA
from utils.ipm import NeuronLoader
from utils.metrics import sepind
from tqdm import tqdm

out_path = dirs['out']['sep']
os.makedirs(out_path, exist_ok=True)

with open("./utils/info.pkl", "rb") as handler:
    info = pickle.load(handler)
info = info[:165]

def mkeqdata(data, groups, n, seed=0):
    np.random.seed(seed)
    
    inds, X, y, l = [], [], [], []
    for i, g in enumerate(np.unique(groups)):
        if 'none'==g:
            continue

        inds.append(np.random.choice(np.argwhere(groups==g).squeeze(), size=n, replace=False))
        X.extend(data[inds[-1]])
        y.extend(i * np.ones_like(inds[-1]))
        l.append(g)
    
    X, y = np.array(X), np.array(y)
    rind = np.random.permutation(X.shape[0])
    X, y = X[rind], y[rind]
    return X, y, l

nmb_rep = 50
monkey ="both"
selectivity = "fast"

data = NeuronLoader(f'G:/Data/{selectivity.capitalize()}/{monkey.capitalize()}')
time_bins = np.arange(data.it.shape[2], step=5)

for region in ['pfc']:#['it', 'pfc']:
    region_data = data.it if region=='it' else data.pfc
    
    for ordination in ['sup']:#['sup', 'mid']:
        nmb_smp = info[ordination].groupby(info[ordination]).count().drop('none').values.min()
        nmb_cls = info[ordination].groupby(info[ordination]).count().drop('none').size
        nmb_bin = time_bins.size

        separabilityIndex = [];
        for seed in (pbar := tqdm(range(nmb_rep))):
            X, y, l = mkeqdata(region_data, info[ordination].to_numpy(), nmb_smp, seed=seed)

            si, _ = sepind(X, y)
            separabilityIndex.append(si)
        separabilityIndex = np.array(separabilityIndex)
        np.save(os.path.join(out_path, f"{monkey}-{selectivity}-{region}-{ordination}.npy"), np.array(separabilityIndex))
