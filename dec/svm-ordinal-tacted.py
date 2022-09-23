import json
import os
import pickle
import sys

import numpy as np
from tqdm import tqdm

with open('./dirs.json', 'r') as f:
    dirs = json.load(f)
sys.path.append(dirs['root'])

from sklearn.model_selection import train_test_split as tts
from utils.ipm import NeuronLoader
from utils.metrics import DescriminationConfidenceEstimator as DCE
from utils.utils import mkeqdata

nmb_rep = 50

out_root = os.path.join(dirs['out']['dec'], 'tacted')
os.makedirs(out_root, exist_ok=True)

for monkey in ["zebel"]:
    for selectivity in ["fast"]:
        data = NeuronLoader(f'G:/Data/{selectivity.capitalize()}/{monkey.capitalize()}')
        outPath  = os.path.join(out_root, f'{monkey.lower()}-{selectivity.lower()}')
        os.makedirs(outPath, exist_ok=True)
            
        for region in ['it', 'pfc']:
            for sfr in ["BI", "BL", "BH"]:
                with open("./utils/info.pkl", "rb") as handler:
                    info = pickle.load(handler)
                info = info[:165]
                validStimuliIndex = ((info.sfr == sfr)).to_numpy()
                region_data = data.it[:165] if region=='it' else data.pfc[:165]
                region_data = region_data[validStimuliIndex]
                info = info[validStimuliIndex]
                
                for ordination in ['sup']:
                    nmb_smp = info[ordination].groupby(info[ordination]).count().drop('none', errors='ignore').values.min()
                    nmb_cls = info[ordination].groupby(info[ordination]).count().drop('none', errors='ignore').size
                
                    cfn, dpr = [], []

                    for seed in (pbar := tqdm(range(nmb_rep))):
                        _cfn, _dpr = [], []
                        pbar.set_description(sfr+"-"+ monkey+"-"+selectivity+"-"+region+"-"+ordination)

                        X, y, l = mkeqdata(region_data, info[ordination].to_numpy(), nmb_smp, seed=seed)
                        trainIndex, testIndex = tts(np.arange(len(y)), test_size=.3, random_state=seed, shuffle=True, stratify=y)
                        XTrain, XTest = X[trainIndex], X[testIndex]
                        yTrain, yTest = y[trainIndex], y[testIndex]

                        for itime in np.arange(len(data.time)):
                            out = DCE(random_state=seed).fit(X_train=XTrain[:, :, itime], y_train=yTrain).score(XTest[:, :, itime], yTest)
                            _cfn.append(out[0])
                            _dpr.append(out[2])
                        
                        cfn.append(_cfn)
                        dpr.append(_dpr)

                    prefix = f"{sfr.lower()}-{region}-{ordination}-"
                    np.save(os.path.join(outPath, prefix + 'cfn.npy'), np.array(cfn))
                    np.save(os.path.join(outPath, prefix + 'dpr.npy'), np.array(dpr))
                    with open(os.path.join(outPath, prefix + 'meta.npy'), "wb") as handler:
                        pickle.dump({'labels': l, 'time': data.time, 'nmb_cls': nmb_cls, 'nmb_smp': nmb_smp, 'info': info}, handler)
