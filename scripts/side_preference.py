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
import pandas as pd

from src import (
    fly_arena as fa
)


# %%
p_mat = Path('./right_left/starved-bb/starved-bb-track.mat')
p_sep = './right_left/coordinates.txt'
p_frame = './right_left/frame.png'

p_out = p_mat.parent / f'{p_mat.stem}_side_analysis'
p_out.mkdir(exist_ok=True)


# %%

data = fa.load_track(p_mat)

sep = fa.load_separation(p_sep)

line = fa.get_line(sep)

results = fa.get_sides(data, line)


# %%

df = pd.DataFrame({
    'frames_left' : pd.Series(dtype=int),
    'frames_right': pd.Series(dtype=int),
    'ratio' : pd.Series(dtype=float),
    'dropped_frames' : pd.Series(dtype=int),
})


df.index.name = 'fly'

for i, res in results.items():
    
    fly = i + 1

    fa.plot_sides(p_frame, sep, line, res, p_out / f'fly_{fly}.png')

    fl, fr = res['left_mask'].sum(), res['right_mask'].sum()
    ratio = fr / fl if fl > 0 else 0
    nan = res['nan_frames']

    df.loc[fly, :] = [ fl, fr, ratio, nan ]

df

# %%
