# Diffusion-Accessible-Domains: Simulation

This repository contains all the files required to replicate the figure 3 results from Singh and Chakrabarti (`https://www.biorxiv.org/content/10.1101/2022.08.31.505992v1.full`). 

## Installation

See this on how to install Conda -> `https://conda.io/projects/conda/en/latest/user-guide/install/index.html`

To create an environment using the .yml file in this repository:
```bash
conda env create --file=environment.yml
```

## Usage
After installing `conda`, follow the steps below to run `script.r`:
```bash
# Clone this repository to your local computer
git clone https://github.com/Shaonlab/Simulation-DADs.git

# Create a conda environment using the provided yml file
conda env create --file=/path/to/environment.yml

# Activate the conda environment
conda activate hic_analysis
```
To run the `simulation.r` script:
```bash
source("simulation.r")
```
