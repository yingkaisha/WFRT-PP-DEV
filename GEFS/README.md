# GEFS

This folder contains Python and Bash scripts that access and pre-process the operational GEFS.

**`GEFS_download.py`**

Download the GEFS ensemble mean fields from the NOAA Operational Model Archive and Distribution System (NOMADS).

**`main.sh`**

The main routine. Proposed as a daily cron job at the WFRT servers. 

# Usage

1. Modify the `namelist.py` to specify the target driectory, NOMADS urls, and file prefix.
3. Execute `main.sh`

# Dependency
* Bash
* Python 3.6

# Contributors

* Kyle Sha <yingkai@eoas.ubc.ca>
* Roland Schigas

# Overview

