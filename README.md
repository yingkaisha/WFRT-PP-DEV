# The GEFS total precipitation post-processing system

This project aims to post-process the 0.25 degree gridded GEFS ensemble mean; produce bias-corrected and calibrated total precipitation forecasts in British Columbia.

The system is operated daily at tombstone.eos.ubc.ca; its contains four sections: (1) GEFS raw forecast download; (2) Analog Ensemble (AnEn) and Minimum Divergence Schaake Shuffle (MDSS) post-processing; (3) Convolutional Neural Network (CNN) based post-processing; (4) Data visualization and file cleaning. 

Technical details of the system are documented in:

* Sha, Y., D. J. Gagne II, G. West, and R. Stull, 2021 (submitted): A hybrid analog-ensemble, convolutional-neural-network method for post-processing precipitation forecasts. Mon Wea. Rev.

## The main folder

**`rotation_ubc.sh`**

The bash script that performs the entire post-processing routine

**`main_nonDL.py`**

The Python script that performs AnEn and MDSS post-processing.

**`main_DL.py`**

The Python script that perfroms CNN-based post-processing; it takes the output of `main_nonDL.py` as input.

**`main_plot.py`**

The data visualization script.

**`namelist_ubc.py`** and **`namelist_plot`**

Configuration files that contain the hyperparameter, directory, and metadata information of the post-processing system.

**`requirements.txt`**

The script that contains the major Python dependencies of the system.

### Folder `GEFS`

Bash and Python scripts for downloading the GEFS files from NOMADS.

### Folder `libs`

Functions and modules.

