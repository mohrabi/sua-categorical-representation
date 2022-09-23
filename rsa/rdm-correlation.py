import json
import os
import pickle
import sys

import numpy as np
from tqdm import tqdm

with open('./dirs.json', 'r') as f:
    dirs = json.load(f)
sys.path.append(dirs['root'])

from utils.ipm import NeuronLoader
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import kendalltau

nmb_rep = 100

outPath = os.path.join(dirs['out']['rsa'])
os.makedirs(outPath, exist_ok=True)

def hist_equalize(matrix):
    matrix = np.copy(matrix)
    shape = matrix.shape
    matrix = matrix.reshape(-1)
    index = matrix.argsort()
    values = np.linspace(0, 1, index.shape[0])
    matrix[index] = values
    return matrix.reshape(shape)

def mkrdm(data):
    rdmat = 1 - cosine_similarity(data)
    return rdmat

def mkherdm(data):
    rdmat = 1 - cosine_similarity(data)
    return hist_equalize(rdmat)

monkey = "both"
selectivity = "fast"
data = NeuronLoader(f'G:/Data/{selectivity.capitalize()}/{monkey.capitalize()}')

with open("./utils/info.pkl", "rb") as handler:
    info = pickle.load(handler)
info = info[:165]

noneIndex = (info.cat == "none")
info = info[~noneIndex].reset_index(drop=True)
data.it = data.it[~noneIndex]
data.pfc = data.pfc[~noneIndex]

argsort = np.argsort(info.cat)
rdmItc = np.nan * np.ones([argsort.size, argsort.size, data.time.size])
rdmPfc = np.nan * np.ones([argsort.size, argsort.size, data.time.size])

for itime in range(data.time.size):
    rdmItc[:,:,itime] = mkrdm(data.it[argsort, :, itime])
    rdmPfc[:,:,itime] = mkrdm(data.pfc[argsort, :, itime])

tt = np.nan * np.zeros((850, 850))
for itime in tqdm(range(data.time.size)):
    for jtime in range(data.time.size):
        r, _ = kendalltau(rdmItc[:, :, itime], rdmPfc[:, :, jtime])
        tt[itime, jtime] = r

np.save(os.path.join(outPath, "rdm-correlation-between-it-and-pfc.npy"), tt)