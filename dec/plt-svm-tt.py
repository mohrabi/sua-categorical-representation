import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker


def myfmt(x, pos):
    return '{0:.1f}'.format(x)

titlefontsize = 32
fontsize = 24

monkey = "both"
selectivity = "fast"
out_path = f"G:/Codes/Processing/out/svm-tt/{monkey}-{selectivity}"

metric = "dpr"
for region in ["it", "pfc"]:
    fig, axs = plt.subplots(1, 2, figsize=np.array([30, 13]), dpi=30)#/2.54)
    time = np.load(os.path.join(out_path, f"{region.lower()}-super-ordinate-tim.npy"))
    for ax, ordination in zip(axs, ["super-ordinate", "mid-level"]):
        m = np.load(os.path.join(out_path, f"{region.lower()}-{ordination}-{metric}.npy"))
        h = ax.contourf(time, time, m.mean(2), 400, cmap='inferno')
        cbar = fig.colorbar(h, ax=ax, shrink=.6, format=ticker.FuncFormatter(myfmt))
        cbar.ax.tick_params(labelsize=fontsize)
        ax.set_title(ordination, fontsize=fontsize)

        ax.axis("square")
        ax.set_xlabel('train time (ms)', fontsize=fontsize)
        ax.set_ylabel('test time (ms)', fontsize=fontsize)
        ax.tick_params(axis='both', which='major', labelsize=fontsize)

    fig.suptitle(f"{region.upper()}: d-prime of Distance of Observations from Hyper-plane", fontsize=titlefontsize)
    fig.tight_layout()
    plt.savefig(os.path.join(out_path, f"{region}-{metric}.jpg"), dpi=1000, bbox_inches='tight')
    plt.close(fig=fig);

metric = "cfn"
for region in ["it", "pfc"]:
    fig, axs = plt.subplots(1, 2, figsize=np.array([30, 13]), dpi=30)#/2.54)
    time = np.load(os.path.join(out_path, f"{region.lower()}-super-ordinate-tim.npy"))
    for ax, ordination in zip(axs, ["super-ordinate", "mid-level"]):
        m = np.load(os.path.join(out_path, f"{region.lower()}-{ordination}-{metric}.npy"))
        m = np.trace(m, axis1=0, axis2=1) / m.sum((0, 1))
        h = ax.contourf(time, time, m.mean(2), 400, cmap='seismic')
        cbar = fig.colorbar(h, ax=ax, shrink=.6, format=ticker.FuncFormatter(myfmt))
        cbar.ax.tick_params(labelsize=fontsize)
        ax.set_title(ordination, fontsize=titlefontsize)

        ax.axis("image")
        ax.set_xlabel('train time (ms)', fontsize=fontsize)
        ax.tick_params(axis='both', which='major', labelsize=fontsize)    

    axs[0].set_ylabel('test time (ms)', fontsize=fontsize)

    fig.suptitle(f"{region.upper()}: Decoder Accuracy", fontsize=titlefontsize)
    fig.tight_layout()
    plt.savefig(os.path.join(out_path, f"{region}-{metric}.jpg"), dpi=1000, bbox_inches='tight')
    plt.close(fig=fig);
