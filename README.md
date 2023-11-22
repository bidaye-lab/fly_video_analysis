# Fly Video Analysis

Collection of workflows, scripts, pipelines etc that process and analyze recordings of flies.


# Installation

```
# create conda environment with necessary dependencies
conda create -n video_analysis -f environment.yml

# get source code
git clone https://github.com/bidaye-lab/fly_video_analysis

# install code as local local python module
cd fly_video_analysis
pip install -e .
```


# jupyter notebooks and git
Jupyter notebooks are saved as `JSON` file,
which do not work well with git version control.
Therefore, no `.ipynb` notebooks is stored in this repo.
Instead, they are stored as `.py` python script files in the [`scripts/`](./scripts/) folder.
Note that the `.py` files do not contain any output cells.

To work with notebooks, one has to convert the `.py` files back to `.ipynb`.
This is done via the [`jupytext`](https://jupytext.readthedocs.io/en/latest/index.html) python module.
To approaches are outlined below.

## option 1: jupytext vscode extension (recommended for vscode)
The vscode extension [Jupytext for Notebooks (congyiwu)](https://marketplace.visualstudio.com/items?itemName=congyiwu.vscode-jupytext)
adds the option "Open as Jupyter Notebook" to the `.py` files (right click on the file in the explorer view).
The `.ipynb` file is never actually stored on disc,
but any changes made to the `.ipynb` file are directly written to the `.py` file.

## option 2: jupytext CLI
`jupytext` can be called from the command line to sync between `.ipynb` and `.py` file.
The [configuration file](./pyproject.toml) defines `scripts` as the folder to store `.py` files,
and `notebooks` folder to create the corresponding `.ipynb` files.

Calling `jupytext` on a `.py` file will create the corresponding `.ipynb` file in the `notebooks` folder:
```
jupytext --sync scripts/example.py
[jupytext] Reading scripts/example.py in format py
[jupytext] Updating notebooks/example.ipynb
[jupytext] Updating the timestamp of scripts/example.py
```
Calling `jupytext` again on either the `.py` or the `.ipynb` file will update the contents of the older file with the newer ones.
Note that the `notebooks` folder is included in the [`.gitignore`](./.gitignore),
because notebooks should not be included in version control.
It is the users responsibility to make sure that the changes in the `.ipynb` files are reflected in the `.py` file before comits or merges.

