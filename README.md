# MIROSLAV analysis

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/weiji14/deepbedmap/)

### What is it?

This repo holds a complete, user-friendly MIROSLAV software toolkit for analysis of MIROSLAV data.

You can run the entire Python/R workflow from your browser using Google Colab:

`url`

### What is MIROSLAV, anyway?

**MIROSLAV (_Multicage InfraRed Open Source Locomotor Activity eValuator_)** is a platform for non-invasive monitoring of circadian locomotor activity in laboratory rodents. All of its hardware and software components are described in the paper: #url

This repository holds the software for all stages of MIROSLAV data processing:

| **Step** 	|   **Processing tool**   	| **What does it do?**                                                                                                                                                      	| **What does it output?**                                                                                                                	|
|:--------:	|:-----------------------:	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|-----------------------------------------------------------------------------------------------------------------------------------------	|
|   **1**  	|   **_Prepare-a-SLAV_**  	| Processes raw textual logs into a data structure, labels the data, resamples the data to reduce size and facilitate subsequent analyses.                                  	| A tabular dataframe containing labelled sensor readings.                                                                                	|
|   **2**  	|      **_TidySLAV_**     	| Melts Prepare-a-SLAV’s table into a tall format, performs standardisation against baseline values, detects periods when the sensors were disconnected, and adds metadata. 	| A tall dataframe with experimental metadata, suitable for plotting and statistical analysis.                                            	|
|   **3**  	|      **_MIROSine_**     	| Operates on TidySLAV’s dataframe, reduces each sensor’s readings to daily 24-hour rhythm amplitudes, mean activities (MESORs), and times of peak activity (acrophases).   	| A dataframe containing three parameters for every day of each sensor’s recordings, describing the 24-hour rhythm of circadian activity. 	|
|  **3.1** 	| **_MIRO The Explorer_** 	| Generate exploratory plots showing temporal dynamics of the three MIROSine parameters over the course of the experiment                                                   	| Exploratory plots                                                                                                                       	|
|  **3.2** 	|    **_StatistiSLAV_**   	| Statistical comparisons of treated groups in specified timepoints                                                                                                         	| Plots and a results Excel document                                                                                                      	|

<sup><sup>_Table adapted from our MIROSLAV paper_</sup></sup>

### How to utilise this workflow?

#### Google Colab

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/weiji14/deepbedmap/) 
The entire workflow can be run in your web browser. Simply click the Colab button, which will lead you to an interactive environment where you can run the notebooks, and clone the entire workflow to adapt for your own logs' analysis.

#### Locally, on your computer

Getting the analysis running locally is a bit more complicated. There are many ways to do it:

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

**Note:** The current version of `glmmTMB` (at the time of writing), `1.1.9`, has a bug when modelling Pearson type VII distributions `t_family`. [We have reported this bug and it was fixed promptly](https://github.com/glmmTMB/glmmTMB/issues/1024), but it will probably we released with the next version of `glmmTMB`.
In the meanwhile, users can install the GitHub version:

1. Install the `remotes` package:

    ```R
    install.packages("remotes")
    ```

2. Install the development version of `glmmTMB` directly from the GitHub repo:

    ```R
    remotes::install_github("glmmTMB/glmmTMB")
    ```

### License

All code is distributed under the GPLv3 license.
