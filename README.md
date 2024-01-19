Collection of workflows, scripts, pipelines etc that process and analyze recordings of flies.

# How to use this repo
For more information on the structure of this repo, 
see this [template repo](https://github.com/bidaye-lab/template_data_pipelines).

## Analysis pipelines
The script files `scripts/*.py` are workflows for the individual steps in the analysis pipeline.

|script file|use case|
|---|---|
|`side_preference.py`| analyze time and stopping behavior in left and right side of chamber|

## old scripts
These old scripts need to become part of the new code structure.
- `scripts/still_frames/`


## Installation
```
# create conda environment with necessary dependencies
conda env create -n free_walking_analysis -f environment.yml

# get source code
git clone https://github.com/bidaye-lab/free_walking_analysis

# install code as local local python module
cd free_walking_analysis
pip install -e .
```
