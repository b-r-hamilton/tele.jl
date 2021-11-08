using NCDatasets, Revise, Statistics, PyPlot, DelimitedFiles
using tele
using Dates, DSP

data_dir = "/home/brynn/Documents/JP/860/tsa_data/"
oc2k_p = data_dir * "Ocean2kLR2015/input_data/Ocean2kLR2015sst.csv"
cesm_p1 = data_dir * "tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc"
cesm_p2 = data_dir * "tos_Omon_CCSM4_past1000_r1i1p1_135001-185012.nc"

#READ OCEAN2K Data
mat = readdlm(oc2k_p, ',')
c1 = 109 #age column index Shevenell TS - got this from excel sheet
c2 = 110
age = mat[begin:17, c1]
d_age = [DateTime(x) for x in age]
val = mat[begin:17, c2]

savetype = ".png"

#READ NETCDFS TO BOUNDING BOX

#Shevenell = -64.87, -64.2

bbox = [-63, -65, 295, 297]
(latc, lonc, timec1, datac1) = read_netcdf(cesm_p1, "tos", bbox)
(latc, lonc, timec, datac) = read_netcdf(cesm_p2, "tos", bbox)
datac = datac .- 273.15 #convert to celsius
contplot(latc,lonc,temp_mean(datac), "images/oc2k"*savetype)

#FILTER CHARACTERISTICS
responsetype = Lowpass(0.1*10^(-7),fs = 3.80265176 * 10^(-7))
designmethod = Butterworth(4)

#cesm
timec = reinterpret(DateTime, timec) #convert from DatetimeNoLeap to DateTime

anomc = anomaly(datac, timec)
tavgc_anom = spatial_mean(anomc)
tavgc_raw = spatial_mean(datac)


using PyPlot
tavgc_raw_a, years = yearly_mean(timec, tavgc_raw)
years_d = [DateTime(x) for x in years]

figure()
plot(years_d, tavgc_raw_a, alpha = 0.5)
plot(d_age, val, ".-")
legend(["CESM Annual Average", "Shevenell et al"])
xlabel("Time (years)")
ylabel("SST [°C]")
title("Alkenone Reconstruction of SST vs CESM Annual Average")
savefig("images/oc2k_cesm"*savetype)
