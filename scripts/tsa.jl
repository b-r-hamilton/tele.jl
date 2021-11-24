using tele
using NCDatasets, Revise, Statistics, PyPlot
using CFTime, DSP, FFTW
using DataFrames, NaNMath
#analyzing data from https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html

#FILE PATHS
data_dir = "/home/brynn/Documents/JP/860/tsa_data/"
ersst_p = data_dir * "sst.mnmean.nc"
# cesm_p2 = data_dir * "tos_Omon_CCSM4_past1000_r1i1p1_135001-185012.nc"
cesm_p2 = data_dir * "tos_Omon_CCSM4_historical_r1i1p1_185001-200512.nc"


#READ NETCDFS TO BOUNDING BOX
bbox = [71, 55, 290, 350]
(late, lone, timee, datae) = read_netcdf(ersst_p, "sst", bbox)
(latc, lonc, timec, datac) = read_netcdf(cesm_p2, "tos", bbox)

#FILTER CHARACTERISTICS
responsetype = Lowpass(0.1*10^(-7),fs = 3.80265176 * 10^(-7))
designmethod = Butterworth(4)

#COMPUTE TIME SERIES (MEAN OVER BBOX, MEAN ANOMALY OVER BBOX, FILTER EACH)
#ersst
anome = anomaly(datae, timee)
tavge_anom = spatial_mean(anome)
tavge_raw = spatial_mean(datae)
tavge_anom_filt = filt(digitalfilter(responsetype, designmethod), tavge_anom)
tavge_raw_filt = filt(digitalfilter(responsetype, designmethod), tavge_raw)

#cesm
timec = reinterpret(typeof(timee[begin]), timec) #convert from DatetimeNoLeap to DateTime
# datac = Float32.(datac) #convert to Float32 (used to be Union{Missing, Float32} but our missing are gone anyways)
datac = datac .- 273.15 #convert to celsius

anomc = anomaly(datac, timec)
tavgc_anom = spatial_mean(anomc)
tavgc_raw = spatial_mean(datac)
tavgc_anom_filt = filt(digitalfilter(responsetype, designmethod), tavgc_anom)
tavgc_raw_filt = filt(digitalfilter(responsetype, designmethod), tavgc_raw)

#Make contour plots
savetype = ".png"
contplot(latc,lonc,temp_mean(datac), "images/datac"*savetype)
contplot(late,lone,temp_mean(datae), "images/datae"*savetype)
contplot(latc,lonc,temp_mean(anomc), "images/anomc"*savetype)
contplot(late,lone,temp_mean(anome), "images/anome"*savetype)

#Make time series plots
filtseriesplot(tavgc_raw, tavgc_raw_filt, timec, tavge_raw, tavge_raw_filt, timee, "", "images/raw"*savetype)
filtseriesplot(tavgc_anom, tavgc_anom_filt, timec, tavge_anom, tavge_anom_filt, timee, "", "images/anom"*savetype)

#COMPUTE FFTs
ffte_raw = fft(tavge_raw)
ffte_anom = fft(tavge_anom)
ffte_freq = FFTW.fftfreq(length(timee), 3.80265176 * 10^(-7))

fftc_raw = fft(tavgc_raw)
fftc_anom = fft(tavgc_anom)
fftc_freq = FFTW.fftfreq(length(timec), 3.80265176 * 10^(-7))

#PLot FFT
fftplot(ffte_freq[:], [ffte_raw, ffte_anom], "images/ffte"*savetype, "ERSST")
fftplot(fftc_freq[:], [fftc_raw, fftc_anom], "images/fftc"*savetype, "CESM")

#Plot FFT
figure()
halfway = 1006
plot(ffte_freq[begin+1:halfway], abs.(ffte_raw)[begin+1:halfway], alpha = 0.5)
halfway = 936
plot(fftc_freq[begin+1:halfway], abs.(fftc_raw)[begin+1:halfway], alpha = 0.5)
ylabel("Frequency [Hz]")
xlabel("Intensity")
legend(["ERSST", "CESM"])
savefig("images/fftec"*savetype)
close()


#histogram
#fig2b - boxplotGPCC Total # of Precips Obs
figure(figsize = (6, 8))
boxdict = boxplot([tavge_raw, tavgc_raw], notch = true, labels =  ["ERSST", "CESM (1850 - 2005)"])
grid(true)
ylabel("SST [°C]")
savefig("images/boxes_ersst"*savetype)
close()


#Extreme values analysis
df = DataFrame(min = [], max = [], mean = [], median = [], std = [])
push!(df, character(tavge_raw))
push!(df, character(tavgc_raw))

quant = quantile(tavgc_raw, [0.25, 0.5, 0.75])
thresh = (quant[3] - quant[1])* 1 + quant[3]
ind = findall(>(thresh), tavgc_raw)

extreme1 = temp_mean(datac[:,:, ind])
extreme2 = temp_mean(anomc[:,:, ind])

#Pot extreme values
contplot(latc,lonc,extreme1, "images/ext_raw"*savetype)
contplot(latc,lonc,extreme2, "images/ext_anom"*savetype)
