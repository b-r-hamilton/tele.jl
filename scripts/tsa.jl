using tele
using NCDatasets, Revise, Statistics, PyPlot
using CFTime, DSP, NaNMath, FFTW
using DataFrames
#analyzing data from https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html

#FILE PATHS
data_dir = "/home/brynn/Documents/JP/860/tsa_data/"
ersst_p = data_dir * "sst.mnmean.nc"
cesm_p2 = data_dir * "tos_Omon_CCSM4_past1000_r1i1p1_135001-185012.nc"

#READ NETCDFS TO BOUNDING BOX
bbox = [-72, -60, 300, 315]
(late, lone, timee, datae) = read_netcdf(ersst_p, "sst", bbox)
(latc, lonc, timec, datac) = read_netcdf(cesm_p2, "tos", bbox)

#FILTER CHARACTERISTICS
responsetype = Lowpass(0.1*10^(-7),fs = 3.80265176 * 10^(-7))
designmethod = Butterworth(4)

#COMPUTE TIME SERIES (MEAN OVER BBOX, MEAN ANOMALY OVER BBOX, FILTER EACH)
#ersst
anome = anomaly(datae, timee)
tavge_anom = mean(anome, dims = (1,2))[1,1,:]
tavge_raw = mean(datae, dims = (1,2))[1,1,:]
tavge_anom_filt = filt(digitalfilter(responsetype, designmethod), tavge_anom)
tavge_raw_filt = filt(digitalfilter(responsetype, designmethod), tavge_raw)

#cesm
timec = reinterpret(typeof(timee[begin]), timec) #convert from DatetimeNoLeap to DateTime
replace!(datac, missing=>NaN) #convert missing vals to NaN
datac = Float32.(datac) #convert to Float32 (used to be Union{Missing, Float32} but our missing are gone anyways)
datac = datac .- 273.15 #convert to celsius
anomc = anomaly(datac, timec)
tavgc_anom = mean(anomc, dims = (1,2))[1,1,:]
replace!(datac, NaN=>0) #convert missing vals to NaN
tavgc_raw = mean(datac, dims = (1,2))[1,1,:]
tavgc_anom_filt = filt(digitalfilter(responsetype, designmethod), tavgc_anom)
tavgc_raw_filt = filt(digitalfilter(responsetype, designmethod), tavgc_raw)

#Make contour plots
savetype = ".png"
contplot(latc,lonc,mean(datac, dims = 3)[:,:], "images/datac"*savetype)
contplot(late,lone,mean(datae, dims = 3)[:,:], "images/datae"*savetype)
contplot(latc,lonc,mean(anomc, dims = 3)[:,:], "images/anomc"*savetype)
contplot(late,lone,mean(anome, dims = 3)[:,:], "images/anome"*savetype)

#Make time series plots
filtseriesplot(tavgc_raw[6012 - 2013:end], tavgc_raw_filt[6012 - 2013:end], timec[6012 - 2013:end], tavge_raw, tavge_raw_filt, timee, "Average Time Series in ROI", "images/raw"*savetype)
filtseriesplot(tavgc_anom[6012 - 2013:end], tavgc_anom_filt[6012 - 2013:end], timec[6012 - 2013:end], tavge_anom, tavge_anom_filt, timee, "Average Anomaly Time Series in ROI", "images/anom"*savetype)

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

#Plot
figure()
halfway = 1006
plot(ffte_freq[begin+1:halfway], abs.(ffte_raw)[begin+1:halfway], alpha = 0.5)
halfway = 3006
plot(fftc_freq[begin+1:halfway], abs.(fftc_raw)[begin+1:halfway], alpha = 0.5)
ylabel("Frequency [Hz]")
xlabel("Intensity")
legend(["ERSST", "CESM"])
savefig("images/fftec"*savetype)
close()


#histogram
#fig2b - boxplotGPCC Total # of Precips Obs
figure()
subplot(1,2,1)
title("CESM")
boxplot([tavgc_raw, tavgc_anom], notch = true, labels =  ["Raw", "Anomaly"])
subplot(1,2,2)
title("ERSST")
boxplot([tavge_raw, tavge_anom], notch = true, labels = ["Raw", "Anomaly"])
savefig("images/boxes"*savetype)
close()
