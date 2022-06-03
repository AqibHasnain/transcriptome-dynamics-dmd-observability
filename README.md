# Code and data for "Learning transcriptome dynamics for discovery of optimal genetic reporters of novel compounds".

This code has been tested on operations systems: `macOS 10.14, 10.15, 11`

Python 3.8 has been used to develop the code in this repo using the following packages and their versions: 

- `numpy (1.21.4)`
- `pandas (1.3.4)`
- `matplotlib (3.3.3)`
- `seaborn (0.11.0)`
- `scipy (1.7.2)`
- `sklearn (0.0)`
- `biopython (1.7.8)`
- `bioservices (1.7.11)`
- `cvxpy (1.1.7)`

All code can be run on standard CPUs and has not been tested or adapted to run on GPUs. 

To install this repository run: 

`git clone https://github.com/AqibHasnain/transcriptome-dynamics-dmd-observability.git` 

The typical install time is between 1 and 10 minutes on a "normal" desktop computer.

All necessary data to generate the figures from the manuscript are provided in the `data` directory. 

In the python notebook `main.ipynb`, you can run dynamic mode decomposition and observability analysis on a time-series RNA-seq dataset of 624 genes from the organism Pseudomonas fluorescens SBW25. The code to reproduce Figure 2 is found in `main.ipynb`. The code to reproduce Figure 3 is found in `main.ipynb`. The code to reproduce Figure 4 is found in `reporters.ipynb`. The code to reproduce Figure 5 is found in `virtual_sensors.ipynb`. Finally, the code to reproduce Figure 6 is found in `reporters.ipynb`. All python notebooks call on supplemental python scripts (e.g. `dmd.py` or 'sensor_placement.py') to facilitate readable code and debugging. 

The expected run-time for `main.ipynb` is between 10 and 30 minutes if generating all results from scratch (see doRun and doSave flags in the notebooks). If importing pre-saved results, `main,ipynb` run-time is approximately 1 minute. 

The expected run-time for `reporters.ipynb` is 20 minutes if fitting hill functions to the fluorescence data. If loading pre-saved results, run-time is expected to be 1 minute. 

The expected run-time for `virtual_sensors.ipynb` is 1 minute.  
























