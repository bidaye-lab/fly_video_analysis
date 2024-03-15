import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d

import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.gridspec as gridspec

from scipy.io import loadmat
import h5py

from skvideo.io import vread, vwrite
from PIL import ImageDraw, Image
import imageio

from io import BytesIO
from joblib import Parallel, delayed, parallel_backend

def load_track(matlab_file):
    """Load fly tracker data from matlab file

    Matlab file generated with
    https://github.com/kristinbranson/FlyTracker

    Parameters
    ----------
    matlab_file : path-like
        Path to matlab file

    Returns
    -------
    data : np.ndarray
        Array of shape (n_flies, n_frames, 3) with x and y coordinates and orientation
    """

    try:
        # matlab file pre 7.3
        m = loadmat(matlab_file, squeeze_me=True, struct_as_record=False)
        data = vars(m["trk"])["data"][:, :, [0, 1, 2]]
    except NotImplementedError:
        # matlab file since 7.3
        with h5py.File(matlab_file, "r") as f:
            data = f["trk"]["data"][()].T[:, :, [0, 1, 2]]

    return data


def load_separation(separation_file):
    """File with separation line between left and right side of chamber

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
    """

    sep = np.loadtxt(separation_file).astype(int)

    return sep


def get_line(sep):
    """Get line separating left and right side of chamber based on separation points

    Line is defines a dict with ints for values and keys mapping y to x pixels

    Parameters
    ----------
    sep : np.ndarray
        Array of shape (n_points, 2) with x and y coordinates

    Returns
    -------
    line : dict
        Dictionary with y as keys and x as values
    """

    # define line separating left and right
    line = dict()
    for p1, p2 in zip(sep, sep[1:]):
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        pxl = np.max([dx, dy])

        x = np.linspace(p1[0], p2[0], pxl + 1).astype(int)
        y = np.linspace(p1[1], p2[1], pxl + 1).astype(int)

        line = {**line, **{j: i for i, j in zip(x, y)}}

    return line


def get_sides(data, line):
    """Analyze fly trajectories and assign frames to left or right side of chamber

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
        Array of shape (n_flies, n_frames, 3) with x and y coordinates and orientation
    line : dict
        Dictionary with y as keys and x as values defining line separating left and right

    Returns
    -------
    results : dict
        Dictionary with fly as keys and dict with results as values
    """

    results = dict()

    for i_fly, xyo in enumerate(data):
        # drop nan
        fnan = np.isnan(xyo[:, :2]).any(axis=1)

        # replace nan with previous value
        xyo = pd.DataFrame(xyo).ffill(axis=0).values

        # get x and y
        x = xyo[:, 0]
        y = xyo[:, 1]
        ori = xyo[:, 2]
        ori = -ori # flip angle so right: 0, left: 180, up: 90, down: -90
        ori = np.unwrap(ori, discont=np.pi) # remove jumps in angle

        # masks for left and right of line
        bl = np.array([line[int(j)] >= i for i, j in zip(x, y)])
        br = ~bl

        results[i_fly] = {
            "x_pxl": x,
            "y_pxl": y,
            "ori": ori,
            "left_mask": bl,
            "right_mask": br,
            "nan_frames": fnan.sum(),
        }

    return results


def add_velocity(results, sigma):
    """Add velocity to results dict

    The following keys are added to the results dict:
    - velocity: velocity in pixels per frame
    - velocity_smoothed: velocity in pixels per frame after smoothing with Gaussian filter

    Parameters
    ----------
    results : dict
        Dictionary with fly as keys and dict with results as values
    sigma : float
        Sigma for Gaussian filter to smooth velocity
    """

    for _, res in results.items():
        # get x and y
        x = res["x_pxl"]
        y = res["y_pxl"]

        # get velocity
        vel = np.linalg.norm(np.diff(np.array([x, y]).T, axis=0), axis=1)
        vel_smth = np.linalg.norm(
            np.diff(uniform_filter1d(np.array([x, y]).T, sigma, axis=0), axis=0), axis=1
        )
        res["velocity"] = np.insert(vel, 0, np.nan)
        res["velocity_smoothed"] = np.insert(vel_smth, 0, np.nan)

def _angle_diff_rad(arr):
    '''Calculate change in angle between frames

    Parameters
    ----------
    arr : np.ndarray
        Array of shape (n_frames,) with angle in radians

    Returns
    -------
    diff : np.ndarray
        Array of shape (n_frames,) with change in angle in radians
    '''

    a, b = arr[:-1], arr[1:]
    # make sure that diff is between -180 and +180
    diff = (b - a + np.pi) % (2 * np.pi) - np.pi
    # insert nan for first frame so that arr and diff have same shape
    diff = np.insert(diff, 0, np.nan)

    return diff

def add_angle_change(results):
    '''Add change in angle to results dict

    Parameters
    ----------
    results : dict
        Dictionary with fly as keys and dict with results as values
    '''

    for _, res in results.items():
        ori = res["ori"]
        dori = _angle_diff_rad(ori)
        res["dori"] = dori
  

def count_stops(results, thresh_vel, thresh_walk, thresh_stop):
    """Add stop counts to results dict

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
    """

    for _, res in results.items():
        stop = pd.Series(res["velocity_smoothed"] < thresh_vel)

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
                if res["left_mask"][f_i]:
                    n_r += 1
                    f_r += s.sum()
                else:
                    n_l += 1
                    f_l += s.sum()

                for f in s.index:
                    x = res["x_pxl"][f].astype(int)
                    y = res["y_pxl"][f].astype(int)

        res["n_stop_left"] = n_l
        res["n_stop_right"] = n_r
        res["n_stop_frames_left"] = f_l
        res["n_stop_frames_right"] = f_r
        res["stop_frames"] = stop

def load_video(video_file):
    "Load .mp4 file as numpy array of shape (n_frames, height, width, 3)"
    return vread(str(video_file))

def annotate_video(video, p_video_out, results):
    """Annotate video with stop numbers

    This will annotate each frame in the the chamber recording where a
    fly is stopping. The number of the stopping bout is written on the fly.

    Parameters
    ----------
    video : np.array
        Array of shape (n_frames, height, width, 3) with RGB video
    p_video_out : path-like
        Path to output video
    results : dict
        Dictionary with fly as keys and dict with results as values
    """

    for i, res in results.items():
        # full trajectory
        x, y = res["x_pxl"].astype(int), res["y_pxl"].astype(int)

        # cycle through tab10 cmap
        rgb = tuple([int(f * 255) for f in colors.to_rgb("C{}".format(i))])

        stop = res["stop_frames"]
        split = np.split(stop, np.flatnonzero(np.diff(stop)) + 1)
        n_stops = 0
        for s in split:
            if s.sum():  # only stops
                n_stops += 1
                for frame in s.index:  # add marker to all frames in stop
                    img = Image.fromarray(video[frame])
                    draw = ImageDraw.Draw(img)
                    draw.text((x[frame], y[frame]), str(n_stops), rgb)
                    video[frame] = np.array(img)

    vwrite(str(p_video_out), video)

def _convert_figure_to_array(fig, dpi):
    'Convert matplotlib figure to image array'
    buffer = BytesIO()
    fig.savefig(buffer, dpi=dpi)
    buffer.seek(0)
    return imageio.imread(buffer)

def _pad_nan(arr, d):
    'Add d nan values to the beginning of arr'
    return np.pad(arr, (d,0), mode='constant', constant_values=np.nan)


def _format_angle(angle, _):
    'Ensure that axis labels are between -180 and 180 degrees'
    angle = angle % 360
    if angle > 180:
        angle = angle - 360
    return f'{angle:1.0f}'

def plot_angle_frame(vid, res, frame, prev_frames, borders=50, path=''):
    '''Plot orientation trace, velocity, angle and angle change for a single frame

    Parameters
    ----------
    vid : np.array
        Array of shape (n_frames, height, width, 3) with RGB video
    res : dict
        Dictionary with results for a single fly
    frame : int
        Frame to plot
    prev_frames : int
        How many preceeding frames to plot for each frame
    borders : int, optional
        Number of pixels to add around trace, by default 50
    path : path-like, optional
        If not '', save plot to file, by default ''

    Returns
    -------
    fig : plt.Figure
        Matplotlib figure object of plot
    '''

    # current frame
    img = vid[frame]

    # original traces
    x, y, ori, dori, v, v_smth = res['x_pxl'], res['y_pxl'], res['ori'], res['dori'], res['velocity'], res['velocity_smoothed']

    # pad all arays with nans to allow plotting of the first frames
    x, y, ori, dori, v, v_smth = [ _pad_nan(i, prev_frames) for i in [x, y, ori, dori, v, v_smth]]
    
    # new index of current frame
    frame = frame + prev_frames 

    # indices of preceeding frames
    t = np.arange(-prev_frames, 0)

    # slices
    s = slice(frame-prev_frames+1, frame+1)

    # select data
    x, y, ori, dori, v, v_smth = [i[s] for i in [x, y, ori, dori, v, v_smth]]

    # initialize figure
    fig = plt.figure(figsize=(20,10))
    gs = gridspec.GridSpec(nrows=3, ncols=2, width_ratios=[2, 1])

    # plot frame
    ax = fig.add_subplot(gs[:, 0])
    ax.imshow(img)
    ax.set_title(f'frame {frame:5}')

    # plot original trace
    ax.scatter(x, y, s=10, zorder=10, c=t)
    ax.set_axis_off()

    # add head direction
    def _add_arrow(ax, x, y, angle, color='.2', scale=4):
        dx, dy = np.cos(angle), np.sin(angle)
        ax.arrow(x, y, scale*dx, scale*dy, fc=color, ec=color, head_width=scale*.5, head_length=scale*.5)

    for xi, yi, oi in zip(x, y, ori):
        _add_arrow(ax, xi, yi, oi)

    # trim axis
    x_med = (np.nanmax(x) + np.nanmin(x)) / 2
    y_med = (np.nanmax(y) + np.nanmin(y)) / 2
    width_xy = np.nanmax([
        np.nanmax(x) - np.nanmin(x),
        np.nanmax(y) - np.nanmin(y)
        ])
    dxy = width_xy / 2 + borders
    ax.set_xlim(x_med - dxy, x_med + dxy)
    ax.set_ylim(y_med - dxy, y_med + dxy)

    # next plot: velocity
    ax = fig.add_subplot(gs[0, 1])
    ax.set_title('Velocity (raw and smoothed)')
    ax.scatter(t, v, c=t)
    ax.plot(t, v_smth, c='gray', zorder=10)

    # next plot: angle
    ax = fig.add_subplot(gs[1, 1])
    ax.set_title('Angle')
    ax.scatter(t, np.rad2deg(ori), c=t)
    ax.yaxis.set_major_formatter(_format_angle)

    # next plot: change in angle
    ax = fig.add_subplot(gs[2, 1])
    ax.set_title('Change in angle')
    ax.scatter(t, np.rad2deg(dori), c=t)

    fig.tight_layout()

    if path:
        fig.savefig(path)
    plt.close(fig)
    return fig

def make_angle_video(video, p_video_out, res, prev_frames=50, fps=30, n_jobs=-1):
    '''Make video with angle and velocity plots

    Parameters
    ----------
    video : np.array
        Array of shape (n_frames, height, width, 3) with RGB video
    p_video_out : path-like
        Path to output video
    res : dict
        Dictionary with results for a single fly
    prev_frames : int, optional
        How many preceeding frames to plot for each frame, by default 50
    fps : float, optional
        Frames per second of output video, by default 30
    n_jobs : int, optional
        Number of CPUs to use to generate video frames simultaneously, by default -1
    '''

    def get_image(frame, dpi=100):
        fig = plot_angle_frame(video, res, frame, prev_frames)
        img = _convert_figure_to_array(fig, dpi)
        plt.close(fig)
        del fig
        return img

    with parallel_backend("loky", n_jobs=n_jobs):
        frames = Parallel()(
            delayed(get_image)(frame) for frame in range(len(video))
        )

    with imageio.get_writer(
        p_video_out, format="FFMPEG", fps=fps
    ) as writer:
        for frame in frames:
            writer.append_data(frame)

def plot_trajectory(p_png, sep, line, res, path=""):
    """Plot fly trajectory and separation line

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
    """

    # load first frame
    img = Image.open(p_png)

    fig, ax = plt.subplots()

    # plot photo
    ax.imshow(img)

    # plot line defining points
    ax.scatter(sep[:, 0], sep[:, 1], zorder=99, color="k")

    # plot line definition
    y, x = line.keys(), line.values()
    ax.plot(x, y, color="C3", zorder=98)

    # plot trajectories
    cmap_paired = plt.cm.tab20.colors
    x, y = res["x_pxl"], res["y_pxl"]
    l, r = res["left_mask"], res["right_mask"]

    ax.scatter(x[l], y[l], marker=",", s=1, color=cmap_paired[2 * 0])
    ax.scatter(x[r], y[r], marker=",", s=1, color=cmap_paired[2 * 0 + 1])

    if path:
        fig.savefig(path)
        plt.close(fig)


def plot_velocity_and_stops(res, thresh_vel, xscale=1, path=""):
    """Plot velocity and stops

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
    """

    y = res["velocity"]
    fig, ax = plt.subplots(figsize=(5 * len(y) / 1000 * xscale, 5))

    ax.plot(y, c="C0", label="raw trace")

    y = res["velocity_smoothed"]
    ax.plot(y, c="C1", label="smoothed")
    a, b = np.nanmin(y), np.nanmax(y)
    ymin, ymax = a - 0.1 * (b - a), b + 0.1 * (b - a)

    stop = res["stop_frames"]
    x = stop.index[stop]
    y = np.zeros_like(x) - 0.25
    ax.scatter(x, y, marker=".", color="k")

    ax.axhline(thresh_vel, c="gray", ls="--")

    ax.axhline(0, c="k", lw=0.5)
    ax.set_xlabel("frame")
    ax.set_ylabel("velocity")
    ax.set_ylim(ymin, ymax)
    ax.set_title("stops L: {} R: {}".format(res["n_stop_left"], res["n_stop_right"]))
    ax.margins(x=0)
    ax.legend(loc="upper right")

    if path:
        fig.savefig(path)
        plt.close(fig)


def summary_df(results):
    """Generate summary dataframe from results dict

       Dataframe has fly as index and the following columns:
        - n_frames: number of frames
        - n_frames_left: number of frames on left side of chamber
        - n_frames_right: number of frames on right side of chamber
        - ratio_frames_right_left: ratio of frames on right side of chamber to left side
        - nan_frames: number of frames with nan coordinates
        - n_stops_left: number of stopping bouts on left side of chamber
        - n_stops_right: number of stopping bouts on right side of chamber
        - n_stop_frames_left: number of stopping frames in left side of chamber
        - n_stop_frames_right: number of stopping frames in right side of chamber
        - ratio_stop_frames_left: ratio of stopping frames in left side of chamber to all frames
        - ratio_stop_frames_right: ratio of stopping frames in right side of chamber to all frames
        - stops_per_frame_left: number of stopping bouts per frame on left side of chamber
        - stops_per_frame_right: number of stopping bouts per frame on right side of chamber
        - avg_velocity_raw_left: average velocity in pixels per frame on left side of chamber
        - avg_velocity_raw_right: average velocity in pixels per frame on right side of chamber
        - avg_velocity_smoothed_left: average smoothed velocity in pixels per frame on left side of chamber
        - avg_velocity_smoothed_right: average smoothed velocity in pixels per frame on right side of chamber
        - avg_angle_change: average change in angle in degrees per frame
        - avg_angle_change_left: average change in angle in degrees per frame on left side of chamber
        - avg_angle_change_right: average change in angle in degrees per frame on right side of chamber

    Parameters
    ----------
    results : dict
        Dictionary with fly as keys and dict with results as values

    Returns
    -------
    df : pd.DataFrame
        Summary dataframe
    """

    df = pd.DataFrame(
        columns=[
            "n_frames",
            "n_frames_left",
            "n_frames_right",
            "ratio_frames_right_left",
            "nan_frames",
            "n_stops_left",
            "n_stops_right",
            "n_stop_frames_left",
            "n_stop_frames_right",
            "ratio_stop_frames_left",
            "ratio_stop_frames_right",
            "stops_per_frame_left",
            "stops_per_frame_right",
            "avg_velocity_raw_left",
            "avg_velocity_raw_right",
            "avg_velocity_smoothed_left",
            "avg_velocity_smoothed_right",
            "avg_angle_change",
            "avg_angle_change_left",
            "avg_angle_change_right",
        ]
    )
    df.index.name = "fly"

    for i, res in results.items():
        fly = i + 1

        n_frames = len(res["x_pxl"])
        ml, mr = res["left_mask"], res["right_mask"]

        nfl, nfr = ml.sum(), mr.sum()
        rfrl = nfr / nfl if nfl > 0 else 0
        nan = res["nan_frames"]

        nsl, nsr, nsfl, nsfr = (
            res["n_stop_left"],
            res["n_stop_right"],
            res["n_stop_frames_left"],
            res["n_stop_frames_right"],
        )
        rsfl = nsfl / nfl if nfl > 0 else 0
        rsfr = nsfr / nfr if nfr > 0 else 0
        spfl = nsl / nfl if nfl > 0 else 0
        spfr = nsr / nfr if nfr > 0 else 0

        vl, vr = np.nanmean(res["velocity"][ml]), np.nanmean(res["velocity"][mr])
        vlsm, vrsm = np.nanmean(res["velocity_smoothed"][ml]), np.nanmean(res["velocity_smoothed"][mr])

        dori = res["dori"]
        dori, doril, dorir = np.nanmean(dori), np.nanmean(dori[ml]), np.nanmean(dori[mr])
        dori, doril, dorir = np.rad2deg(dori), np.rad2deg(doril), np.rad2deg(dorir)
        
        df.loc[fly, :] = [
            n_frames,
            nfl,
            nfr,
            rfrl,
            nan,
            nsl,
            nsr,
            nsfl,
            nsfr,
            rsfl,
            rsfr,
            spfl,
            spfr,
            vl, vr,
            vlsm, vrsm,
            dori, doril, dorir,
        ]

    return df
