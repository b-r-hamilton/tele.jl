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

#Compute raw average time series over bounding box
datac = datac .- 273.15 #convert to celsius
tavge_raw = spatial_mean(datae)
tavgc_raw = spatial_mean(datac)

#ersst
dataep = zeros(length(timee), length(lone), length(late))
permutedims!(dataep, datae, [3,1,2])
filte = filt(digitalfilter(responsetype, designmethod), dataep)
permutedims!(datae, filte, [2,3,1])
anome = anomaly(datae, timee)
#filte holds filtered anomalized data

#cesm
timec = reinterpret(typeof(timee[begin]), timec) #convert from DatetimeNoLeap to DateTime

datacp = zeros(length(timec), length(lonc), length(latc))
permutedims!(datacp, datac, [3,1,2])
filtc = filt(digitalfilter(responsetype, designmethod), datacp)
permutedims!(datac, filtc, [2,3,1])
anomc = anomaly(datac, timec)
#filtc holds filtered, anomalized data

#filtered, anomalized time series
tavge_anom = spatial_mean(anome)
tavgc_anom = spatial_mean(anomc)



#Make contour plots
savetype = ".pdf"
contplot(latc,lonc,temp_mean(datac), "images/datac"*savetype)
contplot(late,lone,temp_mean(datae), "images/datae"*savetype)
contplot(latc,lonc,temp_mean(anomc), "images/anomc"*savetype)
contplot(late,lone,temp_mean(anome), "images/anome"*savetype)

#Make time series plots
filtseriesplot(tavgc_raw, tavgc_anom, timec, tavge_raw, tavge_anom, timee, "", "images/timeseries"*savetype)
# filtseriesplot(tavgc_anom, tavgc_anom_filt, timec, tavge_anom, tavge_anom_filt, timee, "", "images/anom"*savetype)

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
ax = subplot()
halfway = 1006
plot(ffte_freq[begin+1:halfway], abs.(ffte_anom)[begin+1:halfway], alpha = 0.5)
halfway = 936
plot(fftc_freq[begin+1:halfway], abs.(fftc_anom)[begin+1:halfway], alpha = 0.5)
text(7 * 10 ^-10, 400, "f = 5.5 × 10^(-10) Hz")
text(7 * 10 ^-10, 380, "T = 47 years")
ax.set_xlim(ffte_freq[begin+1], 2 * 0.1*10^(-7))
ax.set_ylim(0, 500)
xlabel("Frequency [Hz]")
ylabel("Intensity")
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
