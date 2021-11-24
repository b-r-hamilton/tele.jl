using NCDatasets, Revise, Statistics, PyPlot, DelimitedFiles
using tele
using Dates, DSP, DataFrames

data_dir = "/home/brynn/Documents/JP/860/tsa_data/"
oc2k_p = data_dir * "Ocean2kLR2015/input_data/Ocean2kLR2015sst.csv"
cesm_p1 = data_dir * "tos_Omon_CCSM4_past1000_r1i1p1_085001-134912.nc"
cesm_p2 = data_dir * "tos_Omon_CCSM4_past1000_r1i1p1_135001-185012.nc"

#READ OCEAN2K Data
mat = readdlm(oc2k_p, ',')
c1 = 109 #age column index Shevenell TS - got this from excel sheetlen
c2 = 110
age = mat[begin:20, c1]
d_age = [DateTime(x) for x in age]
val = mat[begin:20, c2]
savetype = ".pdf"

#READ NETCDFS TO BOUNDING BOX
#Shevenell = -64.87, -64.2
bbox = [71, 55, 290, 350]

(latc, lonc, timec1, datac1) = read_netcdf(cesm_p1, "tos", bbox)
(latc, lonc, timec, datac) = read_netcdf(cesm_p2, "tos", bbox)
#concatenate two sets of data
datac = cat(datac1, datac; dims = 3)
timec = cat(timec1, timec; dims = 1)
datac = datac .- 273.15 #convert to celsius
anomc = anomaly(datac, timec)

#Plot data - mostly just to make sure we're in the right spot
contplot(latc,lonc,temp_mean(datac), "images/oc2k"*savetype)

#FILTER CHARACTERISTICS
responsetype = Lowpass(0.1*10^(-7),fs = 3.80265176 * 10^(-7))
designmethod = Butterworth(4)

#Compute spatial average of CESM data
timec = reinterpret(DateTime, timec) #convert from DatetimeNoLeap to DateTime
tavgc_raw = spatial_mean(datac)

#Convert to annual average
tavgc_raw_a, years = yearly_mean(timec, tavgc_raw)
years_d = [DateTime(x) for x in years]

#Plot annual average CESM and Shevenell dataset
figure()
plot(years_d, tavgc_raw_a, alpha = 0.5)
plot(d_age, val, ".-")
legend(["CESM Annual Average", "Shevenell et al"])
xlabel("Time (years)")
ylabel("SST [°C]")
title("Alkenone Reconstruction of SST vs CESM Annual Average")
savefig("images/oc2k_cesm"*savetype)

#compute statistics
df = DataFrame(min = [], max = [], mean = [], median = [], std = [])
push!(df, character(tavgc_raw))
push!(df, character(val))

#Make boxplot of statistics
figure(figsize = (6, 8))
boxplot([val, tavgc_raw], notch = true, labels =  ["Shevenell et al", "CESM (850 - 1850)"])
ylabel("SST [°C]")
grid(true)
savefig("images/boxes_shev"*savetype)
close()

#Extreme value analysis
quant = quantile(tavgc_raw, [0.25, 0.5, 0.75])
thresh = (quant[3] - quant[1])* 1.5 + quant[3]
ind = findall(>(thresh), tavgc_raw)

extreme1 = temp_mean(datac[:,:, ind])
extreme2 = temp_mean(anomc[:,:, ind])

#Plot extreme value anomalies
contplot(latc,lonc,extreme2, "images/ext_anom_lm"*savetype)
