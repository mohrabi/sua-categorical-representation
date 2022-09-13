from logging import getLogger
import pickle
import os
import numpy as np
from sklearn.model_selection import train_test_split
# from torch.utils.data import Subset

from .logger import create_logger, PD_Stats

def init_experience(params, *args, dump_params=True):
    """
    Initialize the experience:
    - dump parameters
    - create checkpoint repo
    - create a logger
    - create a panda object to keep track of the training statistics
    - Taken from SwAV source code. Changed the function name
    """

    if not os.path.isdir(params.logPath):
        os.mkdir(params.logPath)

    # dump parameters
    if dump_params:
        pickle.dump(params, open(os.path.join(params.logPath, "params.pkl"), "wb"))

    # create repo to store checkpoints
    # params.checkpoints = os.path.join(params.dumppath, "checkpoints")
    # if not os.path.isdir(params.checkpoints):
    #     os.mkdir(params.checkpoints)

    # create a panda object to log loss and acc
    # training_stats = PD_Stats(
    #     os.path.join(params.dumppath, "stats" + ".pkl"), args
    # )

    # create a logger
    logger = create_logger(
        os.path.join(params.logPath, "train.log"), rank=0
    )
    logger.info("============ Initialized logger ============")
    logger.info(
        "\n".join("%s: %s" % (k, str(v)) for k, v in sorted(dict(vars(params)).items()))
    )
    logger.info("The experiment will be stored in %s\n" % params.logPath)
    logger.info("")
    return logger

# def train_val_dataset(dataset, val_split=0.25):
#     train_idx, val_idx = train_test_split(list(range(len(dataset))), test_size=val_split)
#     datasets = {}
#     datasets['train'] = Subset(dataset, train_idx)
#     datasets['val'] = Subset(dataset, val_idx)
#     return datasets