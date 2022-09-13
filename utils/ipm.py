from ctypes import sizeof
import glob
import os
from scipy.io import loadmat
import numpy as np
import numpy_indexed as npi
import pandas as pd

def movavg(inp: np.array, bl: int, ss: int) -> np.array:
    nbin = int(np.floor((inp.shape[2] - bl) / ss))
    out = np.nan * np.zeros((inp.shape[0], inp.shape[1], nbin))
    for ibin in range(nbin):
        out[:, :, ibin] = np.mean(inp[:, :, int(ibin*ss):int(ibin*ss)+bl], 2)
    return out

class NeuronLoader:
    def __init__(self, root):
        self.it = np.load(os.path.join(root, 'ITC.npy'))
        self.pfc = np.load(os.path.join(root, 'PFC.npy'))
        nbin = self.it.shape[2]
        self.time = np.linspace(-175, 674 if nbin==850 else 974 , nbin)

class DataLoader:
    def __init__(self, root):
        self.root = root
        self.it = []
        self.pfc = []
        self.it_id = []
        self.pfc_id = []
        self.cm_it = []
        self.cm_pfc = []
    
    def load(self):
        for session in glob.iglob(f'{self.root}/20*'):
            data_dir = os.path.join(session, 'trials')
            cm = loadmat(os.path.join(data_dir, 'cm.mat'))['cm']
            
            self.it.append(loadmat(os.path.join(data_dir, 'mu_it.mat'))['ua'])
            self.it_id.append('mu')
            self.cm_it.append(cm)
            for u in glob.iglob(f'{data_dir}/su_it*.mat'):
                self.it.append(loadmat(u)['ua'])
                self.it_id.append('su')
                self.cm_it.append(cm)

            self.pfc.append(loadmat(os.path.join(data_dir, 'mu_pfc.mat'))['ua'])
            self.pfc_id.append('mu')
            self.cm_pfc.append(cm)
            for u in glob.iglob(f'{data_dir}/su_pfc*.mat'):
                self.pfc.append(loadmat(u)['ua'])
                self.pfc_id.append('su')
                self.cm_pfc.append(cm)
        return self
    
    def split_it(self, bl, ss, seed):
        train, test = [], []
        for uit in range(len(self.it)):
            ua = self.it[uit]
            cm = self.cm_it[uit]
            
            unit_train, unit_test = [], []
            for lbl, group in enumerate(npi.group_by(cm).split(ua)):
                np.random.seed(seed)
                
                ind = np.arange(group.shape[0])
                np.random.shuffle(ind)
                train_ind = ind[:int(np.ceil(len(ind)/2))]
                test_ind = ind[int(np.ceil(len(ind)/2)):]
                
                unit_train.append(group[train_ind].mean(0))
                unit_test.append(group[test_ind].mean(0))
            train.append(unit_train)
            test.append(unit_test)
        train, test = np.transpose(np.array(train), [1, 0, 2]), np.transpose(np.array(test), [1, 0, 2])
        return self.movavg(train, bl, ss), self.movavg(test, bl, ss)
    
    def split_pfc(self, bl, ss, seed):
        train, test = [], []
        for uit in range(len(self.it)):
            ua = self.pfc[uit]
            cm = self.cm_pfc[uit]
            
            unit_train, unit_test = [], []
            for lbl, group in enumerate(npi.group_by(cm).split(ua)):
                np.random.seed(seed)
                
                ind = np.arange(group.shape[0])
                np.random.shuffle(ind)
                train_ind = ind[:int(np.ceil(len(ind)/2))]
                test_ind = ind[int(np.ceil(len(ind)/2)):]
                
                unit_train.append(group[train_ind].mean(0))
                unit_test.append(group[test_ind].mean(0))
            train.append(unit_train)
            test.append(unit_test)
        train, test = np.transpose(np.array(train), [1, 0, 2]), np.transpose(np.array(test), [1, 0, 2])
        return self.movavg(train, bl, ss), self.movavg(test, bl, ss)
    
    def movavg(self, inp, bl, ss):
        nbin = int(np.floor((inp.shape[2] - bl) / ss))
        out = np.nan * np.zeros((inp.shape[0], inp.shape[1], nbin))
        for ibin in range(nbin):
            out[:, :, ibin] = np.mean(inp[:, :, int(ibin*ss):int(ibin*ss)+bl], 2)
        return out


def generate_neuron_file_path(root, row):
    folder_name = row['folder name']
    region = row['region']
    unit_type = row['unit type']
    
    if unit_type == 0:
        file_name = f"mu_{region.lower()}.mat"
    else:
        file_name = f"su_{region.lower()}_{unit_type}.mat"
    
    return os.path.join(root, folder_name, 'Trial', file_name)


class NeuronReader:
    def __init__(self, root, from_file=False):
        self.root = root
        self.sessioninfo = pd.read_csv(os.path.join(root, 'selective_sessioninfo.csv'))
        self.fnames = self.sessioninfo.apply(lambda row: generate_neuron_file_path(root, row), axis=1).to_numpy()
        self.from_file = from_file
        self.it = self.load('IT')
        self.pfc = self.load('PFC')
        self.time = np.linspace(-175, 674, self.it.shape[2])

    def load(self, region):
        nunit, nstim = (self.sessioninfo.region==region.upper()).sum(), 165
        nsample = loadmat(self.fnames[0])['ua'].shape[1]
        data = np.nan * np.zeros([nstim, nunit, nsample])

        if self.from_file or ~os.path.isfile(os.path.join(self.root, f'{region}.npy')):
            for iunit, fpath in enumerate(self.fnames[self.sessioninfo.region==region.upper()]):
                ua = loadmat(fpath)['ua']
                cm = loadmat(os.path.join(os.path.dirname(fpath), 'cm.mat'))['cm']
                if cm.size > ua.shape[0]:
                    cm = cm[:ua.shape[0]]
                for istim in np.arange(1, nstim + 1, 1):
                    data[istim-1, iunit, :] = ua[(cm == istim).squeeze(), :].mean(0)
            out = self.movavg(data, 50, 1)
            np.save(os.path.join(self.root, f'{region}.npy'), out)
        else:
            out = np.load(os.path.join(self.root, f'{region}.npy'))
        return out
    
    def movavg(self, inp, bl, ss):
        nbin = int(np.floor((inp.shape[2] - bl) / ss))
        out = np.nan * np.zeros((inp.shape[0], inp.shape[1], nbin))
        for ibin in range(nbin):
            out[:, :, ibin] = np.sum(inp[:, :, int(ibin*ss):int(ibin*ss)+bl], 2) / bl * 1000
        return out


def cvtlabels(hierarchy='super-ordinate'):
    human_face  = np.array([1, 3, *list(range(5,9)), *list(range(102,108)), *list(range(156,166))]) - 1
    monkey_face = np.array([13, 17, 18, 110, *list(range(166,171))]) - 1
    animal_face = np.array([*list(range(10,13)), *list(range(14,17)), 108, 109]) - 1

    human_body  = np.array([19, 20, 22, 23,26, 28, *list(range(111,114))]) - 1
    monkey_body = np.array([*list(range(29,32)), 114]) - 1
    animal_body = np.array([*list(range(32,38)), 115, 116]) - 1

    natural     = np.array([*list(range(38,55)), 72, *list(range(117,123))]) - 1
    artificial  = np.array([*list(range(55,72)), 73, *list(range(123,129))]) - 1

    lbls = np.nan * np.ones([170])
    if (hierarchy=='super-ordinate'):
        lbls[human_face] = 0
        lbls[monkey_face] = 0
        lbls[animal_face] = 0
        lbls[human_body] = 0
        lbls[monkey_body] = 0
        lbls[animal_body] = 0
        lbls[natural] = 1
        lbls[artificial] = 1
    elif (hierarchy=='mid-level'):
        lbls[human_face] = 1
        lbls[monkey_face] = 1
        lbls[animal_face] = 1
        lbls[human_body] = 0
        lbls[monkey_body] = 0
        lbls[animal_body] = 0
        lbls[natural] = np.nan
        lbls[artificial] = np.nan
    
    return lbls

def grablabels(hierarchy='super-ordinate'):
    human_face  = np.array([1, 3, *list(range(5,9)), *list(range(102,108)), *list(range(156,166))]) - 1
    monkey_face = np.array([13, 17, 18, 110, *list(range(166,171))]) - 1
    animal_face = np.array([*list(range(10,13)), *list(range(14,17)), 108, 109]) - 1

    human_body  = np.array([19, 20, 22, 23,26, 28, *list(range(111,114))]) - 1
    monkey_body = np.array([*list(range(29,32)), 114, 36]) - 1
    animal_body = np.array([*list(range(32,36)), 37, 115, 116]) - 1

    natural     = np.array([*list(range(38,55)), 72, *list(range(117,123))]) - 1
    artificial  = np.array([*list(range(55,72)), 73, *list(range(123,129))]) - 1

    lbls = np.nan * np.ones([170])
    if (hierarchy=='super-ordinate'):
        lbls[human_face] = 0
        lbls[monkey_face] = 0
        lbls[animal_face] = 0
        lbls[human_body] = 0
        lbls[monkey_body] = 0
        lbls[animal_body] = 0
        lbls[natural] = 1
        lbls[artificial] = 1
    elif (hierarchy=='mid-level'):
        lbls[human_face] = 1
        lbls[monkey_face] = 1
        lbls[animal_face] = np.nan
        lbls[human_body] = 0
        lbls[monkey_body] = 0
        lbls[animal_body] = np.nan
        lbls[natural] = np.nan
        lbls[artificial] = np.nan
        lbls[36] = 0
        lbls[114] = 0
        lbls[115] = 0
    elif (hierarchy=='face-body'):
        lbls[human_face] = 1
        lbls[monkey_face] = 1
        lbls[animal_face] = np.nan
        lbls[human_body] = 0
        lbls[monkey_body] = 0
        lbls[animal_body] = np.nan
        lbls[natural] = np.nan
        lbls[artificial] = np.nan
        lbls[36] = 0
        lbls[114] = 0
        lbls[115] = 0
    elif (hierarchy=='natural-artificial'):
        lbls[human_face] = np.nan
        lbls[monkey_face] = np.nan
        lbls[animal_face] = np.nan
        lbls[human_body] = np.nan
        lbls[monkey_body] = np.nan
        lbls[animal_body] = np.nan
        lbls[natural] = 0
        lbls[artificial] = 1
    if (hierarchy=='super-ordinate-old'):
        lbls[human_face] = 0
        lbls[monkey_face] = 0
        lbls[animal_face] = 0
        lbls[human_body] = 0
        lbls[monkey_body] = 0
        lbls[animal_body] = 0
        lbls[natural] = 1
        lbls[artificial] = 1
    elif (hierarchy=='mid-level-old'):
        lbls[human_face] = 1
        lbls[monkey_face] = 1
        lbls[animal_face] = 1
        lbls[human_body] = 0
        lbls[monkey_body] = 0
        lbls[animal_body] = 0
        lbls[natural] = np.nan
        lbls[artificial] = np.nan
    elif (hierarchy=='super-ordinate-b'):
        lbls[np.array(list(set(b_series).intersection(set(human_face))))] = 0
        lbls[np.array(list(set(b_series).intersection(set(monkey_face))))] = 0
        lbls[np.array(list(set(b_series).intersection(set(animal_face))))] = 0
        lbls[np.array(list(set(b_series).intersection(set(human_body))))] = 0
        lbls[np.array(list(set(b_series).intersection(set(monkey_body))))] = 0
        lbls[np.array(list(set(b_series).intersection(set(animal_body))))] = 0
        lbls[np.array(list(set(b_series).intersection(set(natural))))] = 1
        lbls[np.array(list(set(b_series).intersection(set(artificial))))] = 1
    elif (hierarchy=='mid-level-b'):
        lbls[np.array(list(set(b_series).intersection(set(human_face))))] = 1
        lbls[np.array(list(set(b_series).intersection(set(monkey_face))))] = 1
        lbls[np.array(list(set(b_series).intersection(set(animal_face))))] = 1
        lbls[np.array(list(set(b_series).intersection(set(human_body))))] = 0
        lbls[np.array(list(set(b_series).intersection(set(monkey_body))))] = 0
        lbls[np.array(list(set(b_series).intersection(set(animal_body))))] = 0
        lbls[np.array(list(set(b_series).intersection(set(natural))))] = np.nan
        lbls[np.array(list(set(b_series).intersection(set(artificial))))] = np.nan
    elif (hierarchy=='categories'):
        lbls[animal_face] = 0
        lbls[human_face] = 1
        lbls[monkey_face] = 2
        lbls[human_body] = 3
        lbls[monkey_body] = 4
        lbls[animal_body] = 5
        lbls[natural] = 6
        lbls[artificial] = 7
    
    # lbls[33] = np.nan # Fish
    lbls[19] = np.nan # Human with pumpkin as head
    lbls[28] = np.nan # Monkey Body with large head: Decoded as both face and body
    lbls[43] = np.nan # Banana
    lbls[69] = np.nan # Face-Liked Jar
    lbls[126]= np.nan # Face-Liked Cellphone
    lbls[63] = np.nan # Face-Liked Kettle
    lbls[61] = np.nan # Face-Liked Milk Bottle
    lbls[44] = np.nan # Pineapple
    return lbls[:165]