#!/bin/bash

scriptpath=$(realpath "$0")
dir=$(dirname "$scriptpath")

cd "$dir" || exit

source ./.venv/Scripts/activate
jupytext --set-formats ipynb,py:percent ./1_Prepare-a-SLAV.py
jupytext --set-formats ipynb,py:percent ./2_TidySLAV.py
jupytext --set-formats ipynb,R:percent ./3_MIROSine.R
jupytext --set-formats ipynb,R:percent ./3-1_MIRO_The_Explorer.R
jupytext --set-formats ipynb,R:percent ./3-2_StatistiSLAV.R
