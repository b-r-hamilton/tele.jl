using tele
using NCDatasets, Revise, Statistics, PyPlot, CFTime, DSP, NaNMath
#analyzing data from https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html

data_dir = "/home/brynn/Documents/JP/860/tsa_data/"
ersst_p = data_dir * "sst.mnmean.nc"
nc = NCDataset(ersst_p)

bbox = [-72, -60, 300, 315]
(late, lone, timee, datae) = read_netcdf(ersst_p, "sst", bbox)
anome = anomaly(datae, timee)

tavge = mean(anome, dims = (1,2))[1,1,:]

#cesm_p1 = data_dir * "tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc"
cesm_p2 = data_dir * "tos_Omon_CCSM4_past1000_r1i1p1_135001-185012.nc"
(latc, lonc, timec, datac) = read_netcdf(cesm_p2, "tos", bbox)
timec = reinterpret(typeof(timee[begin]), timec)
replace!(datac, missing=>NaN)
datac = Float32.(datac)
anomc = anomaly(datac, timec)

# contplot(latc,lonc,mean(anomc, dims = 3)[:,:])
# contplot(late,lone,mean(anome, dims = 3)[:,:])


tavgc = mean(anomc, dims = (1,2))[1,1,:]

figure()
subplot(1,2,2)
plot(timee, tavge)
title("ERSST")
xlabel("time")

responsetype = Lowpass(0.01*2, fs = 1)
designmethod = Butterworth(4)
tavge_filt = filt(digitalfilter(responsetype, designmethod), tavge)
plot(timee, tavge_filt)
legend(["ERSST", "ERSST smoothed"])

subplot(1,2,1)
plot(timec, tavgc)
tavgc_filt = filt(digitalfilter(responsetype, designmethod), tavgc)
plot(timec, tavgc_filt)
title("CESM")
legend(["CESM", "CESM smoothed"])
xlabel("time")
ylabel("SST Anomaly [°C]")

suptitle("Average Anomaly in ROI, Raw and Smoothed (Butterworth w/ 0.02 month^-1 Lowpass Filter)")

using FFTW
ffte = fft(tavge)
ffte_freq = fftfreq(length(timee), 3.80265176 * 10^(-7))
period = (1 ./ ffte_freq)/3600/24/365
figure()
scatter(ffte_freq[begin+1:end-1], abs.(ffte)[begin+1:end-1])
xlim([0, maximum(ffte_freq)])

function character(series)
    max = NaNMath.maximum(series)
    return max
end
