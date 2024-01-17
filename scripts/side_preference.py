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
from pathlib import Path

from src import (
    fly_arena as fa
)


# %%
ps = [Path('./right_left/starved-bb/starved-bb-track.mat')]
print(ps)

# %%
p = ps[0]

data = fa.load_track(p)


# %%
