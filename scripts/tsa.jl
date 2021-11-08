using tele
using NCDatasets, Revise, Statistics, PyPlot
using CFTime, DSP, FFTW
using DataFrames, NaNMath
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
boxplot([tavge_raw, tavgc_raw], notch = true, labels =  ["ERSST", "CESM"])
ylabel("SST [°C]")
savefig("images/boxes"*savetype)
close()


df = DataFrame(min = [], max = [], mean = [], median = [], std = [])
push!(df, character(tavge_raw))
push!(df, character(tavgc_raw))
