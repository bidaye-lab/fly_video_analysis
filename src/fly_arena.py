import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d

import matplotlib.pyplot as plt

from scipy.io import loadmat
import h5py

from PIL import Image


def load_track(matlab_file):
    '''Load fly tracker data from matlab file

    Matlab file generated with 
    https://github.com/jstaf/fly_tracker

    Parameters
    ----------
    matlab_file : path-like
        Path to matlab file

    Returns
    -------
    data : np.ndarray
        Array of shape (n_flies, n_frames, 2) with x and y coordinates
    '''


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
    '''File with separation line between left and right side of chamber

    The line is defined by a series of points drawn for example with ImageJ.
    It does not need to span the whole chamber, but it need to cover the 
    part of the y axis where the flies are walking.

    Parameters
    ----------
    separation_file : path-like
        Path to file with separation line

    Returns
    -------
    sep : np.ndarray
        Array of shape (n_points, 2) with x and y coordinates
    '''
     
    sep = np.loadtxt(separation_file).astype(int)

    return sep

def get_line(sep):
    '''Get line separating left and right side of chamber based on separation points

    Line is defines a dict with ints for values and keys mapping y to x pixels

    Parameters
    ----------
    sep : np.ndarray
        Array of shape (n_points, 2) with x and y coordinates

    Returns
    -------
    line : dict
        Dictionary with y as keys and x as values
    '''

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
    '''Analyze fly trajectories and assign frames to left or right side of chamber

    Results are returned as a dictionary with fly as keys (0-indexed) and the
    following values:
    - x_pxl: x coordinates in pixels after removing nan
    - y_pxl: y coordinates in pixels after removing nan
    - left_mask: boolean mask with True for frames on left side of chamber
    - right_mask: boolean mask with True for frames on right side of chamber
    - nan_frames: number of frames with nan coordinates

    Parameters
    ----------
    data : np.ndarray
        Array of shape (n_flies, n_frames, 2) with x and y coordinates
    line : dict
        Dictionary with y as keys and x as values defining line separating left and right

    Returns
    -------
    results : dict
        Dictionary with fly as keys and dict with results as values
    '''

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


def add_velocity(results, sigma):
    '''Add velocity to results dict

    The following keys are added to the results dict:
    - velocity: velocity in pixels per frame
    - velocity_smoothed: velocity in pixels per frame after smoothing with Gaussian filter

    Parameters
    ----------
    results : dict
        Dictionary with fly as keys and dict with results as values
    sigma : float
        Sigma for Gaussian filter to smooth velocity
    '''

    for _, res in results.items():

        # get x and y
        x = res['x_pxl']
        y = res['y_pxl']

        # get velocity
        vel = np.linalg.norm(np.diff(np.array([x, y]).T, axis=0), axis=1)
        vel_smth = np.linalg.norm(np.diff(uniform_filter1d(np.array([x, y]).T, sigma, axis=0), axis=0), axis=1)

        # add to results
        res['velocity'] = vel
        res['velocity_smoothed'] = vel_smth



def count_stops(results, thresh_vel, thresh_walk, thresh_stop):
    '''Add stop counts to results dict

    Smoothed velocity is used to define stopping and walking bouts
    and can be controlled with the thresh_* parameters.
     
    The following keys are added to the results dict:
    - n_stop_left: number of stopping bouts on left side of chamber
    - n_stop_right: number of stopping bouts on right side of chamber
    - n_stop_frames_left: number of stopping frames in left side of chamber
    - n_stop_frames_right: number of stopping frames in right side of chamber
    - stop_frames: boolean mask with True for stopping frames

    Parameters
    ----------
    results : dict
        Dictionary with fly as keys and dict with results as values
    thresh_vel : float
        Threshold for velocity to define stopping
    thresh_walk : int
        Ignore walking bouts shorter than this
    thresh_stop : int
        Ignore stopping bouts shorter than this
    '''
    
    for _, res in results.items():
        
        stop = pd.Series(res['velocity_smoothed'] < thresh_vel)

        # cycle through stop and walk periods
        split = np.split(stop, np.flatnonzero(np.diff(stop)) + 1)
        for s in split:
            # if walk periods are shorter than thresh_walk, set them to stop
            if not s.sum() and (len(s) < thresh_walk):
                stop.loc[s.index] = True

        # redefine periods and cylce again
        split = np.split(stop, np.flatnonzero(np.diff(stop)) + 1)
        for s in split:
            # if stop intervals are shorter than thresh_stop, set them to walk
            if s.sum() and (len(s) < thresh_stop):
                stop.loc[s.index] = False

        # count stops and write to video
        split = np.split(stop, np.flatnonzero(np.diff(stop)) + 1)
        n_r, n_l = 0, 0
        f_r, f_l = 0, 0
        for s in split:
            if s.sum():
                f_i = s.index[0]
                if res['left_mask'][f_i]:
                    n_r += 1
                    f_r += s.sum()
                else:
                    n_l += 1
                    f_l += s.sum()

                for f in s.index:
                    x = res['x_pxl'][f].astype(int)
                    y = res['y_pxl'][f].astype(int)

        res['n_stop_left'] = n_l
        res['n_stop_right'] = n_r
        res['n_stop_frames_left'] = f_l
        res['n_stop_frames_right'] = f_r
        res['stop_frames'] = stop



def plot_trajectory(p_png, sep, line, res, path=''):
    '''Plot fly trajectory and separation line

    Parameters
    ----------
    p_png : path-like
        Example frame to plot ontop of
    sep : np.ndarray
        Array of shape (n_points, 2) with x and y coordinates
    line : dict
        Dictionary with y as keys and x as values defining line separating left and right
    res : dict
        Dictionary with results for a single fly
    path : path-like, optional
        If not '', save plot to disk and close figure, by default ''
    '''
        
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

def plot_velocity_and_stops(res, thresh_vel, xscale=1, path=''):
    '''Plot velocity and stops

    This is useful to check if stop detection thresholds are set correctly.

    Parameters
    ----------
    res : dict
        Dictionary with results for a single fly
    thresh_vel : float
        Threshold for velocity used when counting stops
    xscale : int, optional
        Scale x axis by this factor, by default 1
    path : path-like, optional
        If not '', save plot to disk and close figure, by default ''
    '''
  
    y = res['velocity']
    fig, ax = plt.subplots(figsize=(5 * len(y) / 1000 * xscale, 5))

    ax.plot(y, c='C0', label='raw trace')

    y = res['velocity_smoothed']
    ax.plot(y, c='C1', label='smoothed')
    a, b = y.min(), y.max()
    ymin, ymax = a - 0.1 * (b - a), b + 0.1 * (b - a)

    stop = res['stop_frames']
    x = stop.index[stop]
    y = np.zeros_like(x) - 0.25
    ax.scatter(x, y, marker='.', color='k')

    ax.axhline(thresh_vel, c='gray', ls='--')

    ax.axhline(0, c='k', lw=0.5)
    ax.set_xlabel('frame')
    ax.set_ylabel('velocity')
    ax.set_ylim(ymin, ymax)
    ax.set_title('stops L: {} R: {}'.format(res['n_stop_left'], res['n_stop_right']))
    ax.margins(x=0)
    ax.legend(loc='upper right')

    if path:
        fig.savefig(path)
        plt.close(fig)

def summary_df(results):
    '''Generate summary dataframe from results dict

       Dataframe has fly as index and the following columns:
        - n_frames_left: number of frames on left side of chamber
        - n_frames_right: number of frames on right side of chamber
        - ratio_frames_right_left: ratio of frames on right side of chamber to left side
        - n_frames_dropped: number of frames with nan coordinates
        - n_stops_left: number of stopping bouts on left side of chamber
        - n_stops_right: number of stopping bouts on right side of chamber
        - n_stop_frames_left: number of stopping frames in left side of chamber
        - n_stop_frames_right: number of stopping frames in right side of chamber
        - ratio_stop_frames_left: ratio of stopping frames in left side of chamber to all frames
        - ratio_stop_frames_right: ratio of stopping frames in right side of chamber to all frames
        - stops_per_frame_left: number of stopping bouts per frame on left side of chamber
        - stops_per_frame_right: number of stopping bouts per frame on right side of chamber

    Parameters
    ----------
    results : dict
        Dictionary with fly as keys and dict with results as values

    Returns
    -------
    df : pd.DataFrame
        Summary dataframe
    '''

    df = pd.DataFrame(columns=[
        'n_frames_left',
        'n_frames_right',
        'ratio_frames_right_left',
        'n_frames_dropped',
        'n_stops_left',
        'n_stops_right',
        'n_stop_frames_left',
        'n_stop_frames_right',
        'ratio_stop_frames_left',
        'ratio_stop_frames_right',
        'stops_per_frame_left',
        'stops_per_frame_right'
        ])
    df.index.name = 'fly'

    for i, res in results.items():
        
        fly = i + 1

        nfl, nfr = res['left_mask'].sum(), res['right_mask'].sum()
        rfrl = nfr / nfl if nfl > 0 else 0
        nan = res['nan_frames']

        nsl, nsr, nsfl, nsfr = res['n_stop_left'], res['n_stop_right'], res['n_stop_frames_left'], res['n_stop_frames_right']
        rsfl = nsfl / nfl if nfl > 0 else 0
        rsfr = nsfr / nfr if nfr > 0 else 0
        spfl = nsl / nfl if nfl > 0 else 0
        spfr = nsr / nfr if nfr > 0 else 0
        df.loc[fly, :] = [ nfl, nfr, rfrl, nan, nsl, nsr, nsfl, nsfr, rsfl, rsfr, spfl, spfr ]

    return df