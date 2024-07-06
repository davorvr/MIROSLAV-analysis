<img src="logo/hexsticker.png" height="170" align="left" hspace=15></img>
<br> MIROSLAV analysis <br>
[![Full paper](https://img.shields.io/badge/Full%20paper-10.1101%2F2024.06.25.600592-%23B12929)](https://doi.org/10.1101/2024.06.25.600592) [![Zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.12191589-%230F81C2)](https://doi.org/10.5281/zenodo.12191589) <br clear="left"/>
========

### What is it?

This repo holds a complete, user-friendly MIROSLAV software toolkit for analysis of MIROSLAV home cage activity and habitat environmental data. Here, you can also find real MIROSLAV data from our lab as an example that the tools operate on.

Everything you need to construct the MIROSLAV device can be found in the [`MIROSLAV-hardware`](https://github.com/davorvr/MIROSLAV-hardware) and [`MIROSLAV-firmware`](https://github.com/davorvr/MIROSLAV-firmware) repositories.

### What is MIROSLAV, anyway?

**MIROSLAV (_Multicage InfraRed Open Source Locomotor Activity eValuator_)** is a platform for non-invasive monitoring of circadian locomotor activity in laboratory rodents. MIROSLAV is fully open source and scalable to hundreds of cages. All of its hardware and software components are described in the paper: https://doi.org/10.1101/2024.06.25.600592

This repository holds the software for all stages of MIROSLAV data processing:

| **Step** | **Processing tool**                                                                                                                                                                                            | **What does it do?**                                                                                                                                                      | **What does it output?**                                                                                                                |
|:--------:|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------|
| **1**    | **_Prepare-a-SLAV_** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/davorvr/MIROSLAV-analysis/blob/main/1_Prepare-a-SLAV.ipynb)         | Processes raw textual logs into a data structure, labels the data, resamples the data to reduce size and facilitate subsequent analyses.                                  | A tabular dataframe containing labelled sensor readings.                                                                                |
| **2**    | **_TidySLAV_** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/davorvr/MIROSLAV-analysis/blob/main/2_TidySLAV.ipynb)               | Melts Prepare-a-SLAV’s table into a tall format, performs standardisation against baseline values, detects periods when the sensors were disconnected, and adds metadata. | A tall dataframe with experimental metadata, suitable for plotting and statistical analysis.                                            |
| **3**    | **_MIROSine_** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/davorvr/MIROSLAV-analysis/blob/main/3_MIROSine.ipynb)                     | Operates on TidySLAV’s dataframe, reduces each sensor’s readings to daily 24-hour rhythm amplitudes, mean activities (MESORs), and times of peak activity (acrophases).   | A dataframe containing three parameters for every day of each sensor’s recordings, describing the 24-hour rhythm of circadian activity. |
| **3.1**  | **_MIRO The Explorer_** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/davorvr/MIROSLAV-analysis/blob/main/3-1_MIRO_The_Explorer.ipynb) | Generate exploratory plots showing temporal dynamics of the three MIROSine parameters over the course of the experiment                                                   | Exploratory plots.                                                                                                                      |
| **3.2**  | **_StatistiSLAV_** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/davorvr/MIROSLAV-analysis/blob/main/3-2_StatistiSLAV.ipynb)           | Statistical comparisons of treated groups in specified timepoints                                                                                                         | Plots and an Excel document containing results.                                                                                         |

<sup>_Table adapted from our MIROSLAV paper._</sup>

#### Environmental monitoring data

Additionally, the [`E_EnviroSLAV` directory](https://github.com/davorvr/MIROSLAV-analysis/tree/main/E_EnviroSLAV) contains tools for parsing and plotting environmental data. Two logs are supplied:

1. [`miroslav-env-rack_M.2024-06-27T14-04-44-916883.gz`](https://github.com/davorvr/MIROSLAV-analysis/tree/main/E_EnviroSLAV/0_raw_env/miroslav-env-rack_M.2024-06-27T14-04-44-916883.gz) - A short test environmental monitoring log from the current MIROSLAV v0.4, as a sample of the current log format which is now stable and will be used for future MIROSLAV versions as well. The log follows current MIROSLAV file naming conventions, detailed in [the Prepare-a-SLAV config file](https://github.com/davorvr/MIROSLAV-analysis/blob/main/1_Prepare-a-SLAV_config.toml).

2. [`0_proto_env.csv.gz`](https://github.com/davorvr/MIROSLAV-analysis/tree/main/E_EnviroSLAV/0_raw_env/0_proto_env.csv.gz) - A log recorded with a MIROSLAV prototype during a real experiment, recorded in an older format, but kept in its original form for posterity.

The parser, [`1_EnviroSLAV_parse.py`](https://github.com/davorvr/MIROSLAV-analysis/tree/main/E_EnviroSLAV/1_EnviroSLAV_parse.py), converts both to `.parquet` files of an equivalent format.

The visualiser, [`2_EnviroSLAV_plot.R`](https://github.com/davorvr/MIROSLAV-analysis/tree/main/E_EnviroSLAV/2_EnviroSLAV_plot.R), can read either, but the example code works on prototype logs as they contain nearly 3 weeks of real data recorded in an animal habitat during an experiment.

### How to utilise this workflow?

#### Google Colab
 
The entire workflow can be run in your web browser. Simply click on the ![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg) buttons, which will lead you to an interactive environment where you can run the notebooks with our data to see a working example, or upload your own and see what you got with MIROSLAV!

#### Locally, on your computer

To run the analysis on your computer, you need Python, R, some dependencies, and this repository. Here are a few quick tips:

* **Windows:**
  * We personally prefer the [WinPython](https://winpython.github.io/) to get an easy-to-setup Python distribution with a lot of useful preinstalled software.
  * We install R from [the official website](https://cran.r-project.org/bin/windows/base/).
  * Finally, we install [VSCode](https://code.visualstudio.com/download) with the necessary extensions:
    * [Python](https://marketplace.visualstudio.com/items?itemName=ms-python.python)
    * [Jupyter](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter) 
    * [R](https://marketplace.visualstudio.com/items?itemName=REditorSupport.r)

* **Linux:**

  * Most Linux distributions come with Python preinstalled. R and VSCode can usually be installed using the distro's package manager.
  * We open up VSCode and install the aforementioned extension.

Finally, we open up the local GitHub repository's directory and get to work! :\)

### Dependencies

The notebooks are adapted to install these dependencies automatically, but they're listed here so you know what you'll be installing.

#### Python

Listed in the `Python_requirements.txt` file:

* `mirofile`
* `pandas>=2.2,<3.0.0`
* `fastparquet`
* `plotly>=5.19`
* `ipywidgets>=8.1.2`

#### R

Listed in the `R_requirements.txt` file:

* `dplyr`
* `lubridate`
* `progress`
* `arrow`
* `tzdb`
* `ggplot2`
* `patchwork`
* `emmeans`
* `glmmTMB`
* `DHARMa`

**Note:** The current version of `glmmTMB` (at the time of writing), `1.1.9`, has a bug when modelling Pearson type VII distributions `t_family`. [We have reported this bug and it was fixed promptly](https://github.com/glmmTMB/glmmTMB/issues/1024), but it will probably we released with the next version of `glmmTMB`. If you're not sure if you have a recent enough version of `glmmTMB`, you can run this code, which will check it for you and install the GitHub version if necessary:

  ```R
  if (packageVersion("glmmTMB") <= "1.1.9") {
    install.packages("remotes")
    remotes::install_github("glmmTMB/glmmTMB")
  }
  ```

### Related repositories

* [`mirofile`](https://github.com/davorvr/mirofile) - Your buddy for dealing with raw MIROSLAV data. Our Python library used by Prepare-a-SLAV to parse raw MIROSLAV logs into a dataframe.
* [`MIROSLAV-hardware`](https://github.com/davorvr/MIROSLAV-hardware) - Everything you need to construct the MIROSLAV device and start monitoring hundreds of rodents' locomotor activity patterns 24/7 in their home cages.
* [`MIROSLAV-firmware`](https://github.com/davorvr/MIROSLAV-firmware) - Contains MIROSLAVino (Arduino firmware) and Record-a-SLAV (Python data acquisition script), everything you need to breathe life into your MIROSLAV device.

### License

You can modify any part of MIROSLAV freely under the **GPLv3** license - if you have any questions, problems, or ideas on how to improve MIROSLAV, feel free to reach out to us, submit a GitHub issue, or a pull request.

The MIROSLAV logo is created by Petra Šoštarić Mužić and licensed under the **CC BY-NC-SA 4.0** license - more details can be found in the repository's [`logo/` directory](logo/).
