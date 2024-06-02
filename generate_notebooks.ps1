$scriptpath = $MyInvocation.MyCommand.Path
$dir = Split-Path $scriptpath

Set-Location "$dir"
& ".\.venv\Scripts\activate.ps1"
jupytext --set-formats ipynb,py:percent .\1_Prepare-a-SLAV.py
jupytext --set-formats ipynb,py:percent .\2_TidySLAV.py
jupytext --set-formats ipynb,R:percent .\3_MIROSine.R
jupytext --set-formats ipynb,R:percent .\3-1_MIRO_The_Explorer.R
jupytext --set-formats ipynb,R:percent .\3-2_StatistiSLAV.R
