# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: free_walking
#     language: python
#     name: python3
# ---

# %%
# %load_ext autoreload
# %autoreload 2

from pathlib import Path
from src import fly_arena as fa


# %% [markdown]
# # Load data
# Here we choose which data files to use for the analysis and 
# create an output directory to store the results.
#
# We need to define the following:
# - `p_mat`: matlab file with tracking data
# - `p_sep`: text file with separation line
# - `p_video`: video recording of chamber
# - `p_frame`: example frame for plotting
# - `p_out`: output directory 
#
#
# Then, we load the data and split the trajectories into left and right side of the chamber.

# %%
# select files necessary for analysis
data_root = Path('../data/side_assay/')

p_mat = data_root / 'starved-bb-track.mat' 
p_sep = data_root / 'coordinates.txt' 
p_video = data_root / 'starved-bb.avi' 
p_frame = data_root / 'frame.png'

# create output directory
p_out = data_root / f'{p_video.stem}_side_analysis'
p_out.mkdir(exist_ok=True)

# load tracking data and separation line
data = fa.load_track(p_mat)
sep = fa.load_separation(p_sep)
line = fa.get_line(sep)

# use separation line analyze trajectories
results = fa.get_sides(data, line)

# %% [markdown]
# # Analyze behavior
#
# The trajectories will tell us how much time the flies spend in each side of the chamber
# and how often they stop.
# We can fine-tune the analysis by changing the following parameters:
# - `thresh_walk`: Minimum number of frames for walking bout
# - `thresh_stop`: Minimum number of frames for stopping bout
# - `thresh_vel`: Threshold for velocity to define stopping in pxl/frame
# - `smooth_vel`: Gaussian kernel width for smoothing velocity
#
# The following outputs are geneated in the `p_out` directory:
# - `fly{n}_trajectories.png`: left/right separation of trajectory for fly n 
# - `fly{n}_stops.png`: velocity traces and stops for fly n
# - `summary.csv`: summary for all flies
# - `video_annotated.avi`: chamber recording showing the stopping bouts

# %%
# thresholds for stop analysis
thresh_stop, thresh_walk = 5, 5
thresh_vel = 0.25
smooth_vel = 15

# extract velocity and count stops
fa.add_velocity(results, sigma=smooth_vel)
fa.count_stops(results, thresh_vel, thresh_walk, thresh_stop)

for fly, res in results.items():
    # plot trajectories
    fa.plot_trajectory(p_frame, sep, line, res, path=p_out / f'fly{fly}_trajectories.png')
    # plot velocity and stops
    fa.plot_velocity_and_stops(res, thresh_vel, path=p_out / f'fly{fly}_stops.png')

# summary statistics
df = fa.summary_df(results)
df.to_csv(p_out / 'summary.csv', index=False)
df


# %% [markdown]
# The next cell generates the annotated video, which takes some time to run.

# %%
# optional: generate annotated video
fa.annotate_video(p_video, p_out / 'video_annotated.avi', results)
