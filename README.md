# tele
Code for 12.860, Assignment 4
Processes CESM, ERSST, and Ocean2k data

## Data Dependencies
- CESM - historical and past1000 (downloaded from WHOI CMIP5 server)
	- tos
- ERSST (downloaded from https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html)
	- sst
- Ocean2k (downloaded from https://data.noaa.gov/dataset/dataset/noaa-wds-paleoclimatology-pages-ocean2k-synthesis-data-set)
	- "Ocean2kLR2015.sst.xlsx" file, Southern0483aShevenell2011_ODP-1098B time series

## Packages
- tele: subset .nc files, compute statistics, FFT, trendline, extreme events, and visualize

## Scripts
- tsa.jl: timeseries analysis for CESM (1850-2005) and ERSST. Both are gridded datasets, and all tele.jl functionality is used
- oc2k.jl: rudimentary timeseries analysis of statistics and extreme events for CESM (850-1850) and Shevenell time series

## Julia Dependencies
- Revise
- Dates
- NCDatasets
- PyPlot
- PyCall
- DataFrames
- Statistics
- CFTime
- DSP
- FFTW
- NaNMath
- DelimitedFiles

## Python Dependencies
- cartopy
- cartopy.features
- matplotlib.patches
- matplotlib.ticer

## Notes
- Discovered Julia can't handle doing statistics on a multidimensional array with sparse NaN values. "skipmissing" does not work for this situation. Wrote a bunch of stupid functions (spatial_mean, spatial_std, tempmean, tempstd) to take care of this. Didn't cry. Definitely developed better practices for handling this than I had in monsoon.jl
- Additionally, Julia/DataFrames won't compute an annual average. Wrote some naive code to do this, could be made better in future.
- Anomaly, FFT code is definitely reusable and works well

## To make package
Run in Julia (anywhere)
```
t = Template(; user="b-r-hamilton", dir="~/Code", authors="Brynnydd Hamilton", julia=v"1.6")
```
and then ```t("tele.jl")```
This will create a package template in Code with the name tele. It will also create a GitHub repo with the name "tele". I had to go onto Github website and manually change the name to "tele.jl" to be able to push. Do a force push to overwrite the current data
**Do not have a hyphen in your package name!!!**
