import os
import sys
import pickle
import json

with open('dirs.json', 'r') as f:
    dirs = json.load(f)
sys.path.append(dirs['root'])

import numpy as np
from utils.ipm import NeuronLoader
from utils.metrics import sepind
from utils.utils import mkeqdata
from tqdm import tqdm

outPath = os.path.join(dirs['out']['sep'], "tacted")
os.makedirs(outPath, exist_ok=True)

nmb_rep = 50
monkey ="both"
selectivity = "fast"

data = NeuronLoader(f'G:/Data/{selectivity.capitalize()}/{monkey.capitalize()}')
time_bins = np.arange(data.it.shape[2], step=5)

for region in ['it', 'pfc']:
    region_data = data.it if region=='it' else data.pfc
    
    for sfr in ["BI", "BL", "BH"]:
        with open("./utils/info.pkl", "rb") as handler:
            info = pickle.load(handler)
        info = info[:165]
        validStimuliIndex = ((info.sfr == sfr)).to_numpy()
        region_data = data.it if region=='it' else data.pfc
        region_data = region_data[validStimuliIndex]
        info = info[validStimuliIndex]

        for ordination in ['sup', 'mid']:
            nmb_smp = info[ordination].groupby(info[ordination]).count().drop('none', errors='ignore').values.min()
            nmb_cls = info[ordination].groupby(info[ordination]).count().drop('none', errors='ignore').size
            nmb_bin = time_bins.size

            separabilityIndex = []
            for seed in (pbar := tqdm(range(nmb_rep))):
                X, y, l = mkeqdata(region_data, info[ordination].to_numpy(), nmb_smp, seed=seed)

                si, _ = sepind(X, y)
                separabilityIndex.append(si)
            separabilityIndex = np.array(separabilityIndex)
            
            prefix = f"{sfr}-{region}-{ordination}-"
            np.save(os.path.join(outPath, f"{prefix}si.npy"), separabilityIndex)
            with open(os.path.join(outPath, prefix + 'meta.npy'), "wb") as handler:
                pickle.dump({'labels': l, 'time': data.time, 'nmb_cls': nmb_cls, 'nmb_smp': nmb_smp, 'info': info}, handler)