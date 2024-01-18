import numpy as np

import matplotlib.pyplot as plt

from scipy.io import loadmat
import h5py

from PIL import Image


def load_track(matlab_file):
    # https://github.com/jstaf/fly_tracker


    try:
        # matlab file pre 7.3
        m = loadmat(matlab_file, squeeze_me=True, struct_as_record=False)
        data = vars(m['trk'])['data'][:, :, [0, 1]]
    except NotImplementedError:
        # matlab file since 7.3
        with h5py.File(matlab_file, 'r') as f:
            data = f['trk']['data'][()].T[:, :, [0, 1]]

    return data
    
def load_separation(separation_file):
     
    sep = np.loadtxt(separation_file).astype(int)

    return sep

def get_line(sep):

    # define line separating left and right
    line = dict()
    for p1, p2 in zip(sep, sep[1:]):
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        pxl = np.max([dx, dy])

        x = np.linspace(p1[0], p2[0], pxl + 1).astype(int)
        y = np.linspace(p1[1], p2[1], pxl + 1).astype(int)

        line = {**line, **{ j: i for i, j in zip(x, y)}}

    return line




def get_sides(data, line):

    results = dict()

    for fly, xy in enumerate(data):

        # drop nan
        fnan = np.isnan(xy).any(axis=1)
        xy = xy[~fnan]

        # get x and y
        x = xy[:, 0]
        y = xy[:, 1]

        # masks for left and right of line
        bl = np.array([line[int(j)] >= i for i, j in zip(x, y)])
        br = ~bl

        results[fly] = {
            'x_pxl': x,
            'y_pxl': y,
            'left_mask': bl,
            'right_mask': br,
            'nan_frames': fnan.sum(),
        }

    return results

def plot_sides(p_png, sep, line, res, path=''):
        
        # load first frame
        img = Image.open(p_png)

        fig, ax = plt.subplots()

        # plot photo
        ax.imshow(img)

        # plot line defining points
        ax.scatter(sep[:,0], sep[:,1], zorder=99, color='k')

        # plot line definition
        y, x = line.keys(), line.values()
        ax.plot(x, y, color='C3', zorder=98)
        
        # plot trajectories
        cmap_paired = plt.cm.tab20.colors
        x, y = res['x_pxl'], res['y_pxl']
        l, r = res['left_mask'], res['right_mask']

        ax.scatter(x[l], y[l], marker=',', s=1, color=cmap_paired[2*0])
        ax.scatter(x[r], y[r], marker=',', s=1, color=cmap_paired[2*0 + 1])

        if path:
            fig.savefig(path)
            plt.close(fig)

