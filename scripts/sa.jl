using tele
using NCDatasets, Revise, Statistics, PyPlot, PyCall
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
#sampling = monthly
#cutoff all signals below 3 year periodicity
responsetype = Lowpass(0.1*10^(-7),fs = 3.80265176 * 10^(-7))
designmethod = Butterworth(4)

#CALCULATE ANOMALY AND FILTER
"""
all this permutedims! nonsense is because my anomaly code assumes lon x lon x time
and filter assumes time x lon x lat
"""

#ersst
dataep = zeros(length(timee), length(lone), length(late))
permutedims!(dataep, datae, [3,1,2])
filte = filt(digitalfilter(responsetype, designmethod), dataep)
permutedims!(datae, filte, [2,3,1])
anome = anomaly(datae, timee)
permutedims!(filte, anome, [3,1,2])
#filte holds filtered anomalized data

#cesm
timec = reinterpret(typeof(timee[begin]), timec) #convert from DatetimeNoLeap to DateTime
datac = datac .- 273.15 #convert to celsius

datacp = zeros(length(timec), length(lonc), length(latc))
permutedims!(datacp, datac, [3,1,2])
filtc = filt(digitalfilter(responsetype, designmethod), datacp)
permutedims!(datac, filtc, [2,3,1])
anomc = anomaly(datac, timec)
permutedims!(filtc, anomc, [3,1,2])
#filtc holds filtered, anomalized data

anomc = anomaly(datac, timec)
anomcp = zeros(length(timec), length(lonc), length(latc))
permutedims!(anomcp, anomc, [3,1,2])
filtc = filt(digitalfilter(responsetype, designmethod), anomcp)

#COMPUTE EOFS
eofs = pyimport("eofs")

#cesm
solver = eofs.standard.Eof(filtc)
savetype = ".pdf"
plot_eofs(solver, lonc, latc, timec, "images/eofc"*savetype)

#ersst
solver = eofs.standard.Eof(filte)
plot_eofs(solver, lone, late, timee, "images/eofe"*savetype)

#CORRELATION PLOT
#at every grid point, correlate the mean time series with the grid point's timeseries
tavgc = spatial_mean(anomc)
tavgc_filt = filt(digitalfilter(responsetype, designmethod), tavgc)
tavge = spatial_mean(anome)
tavge_filt = filt(digitalfilter(responsetype, designmethod), tavge)

xcc = zeros(length(timec)*2 - 1, length(lonc), length(latc))
xce = zeros(length(timee)*2 - 1, length(lone), length(late))
regc = zeros(length(lonc), length(latc))
rege = zeros(length(lone), length(late))

for i in 1:length(lonc)
    for j in 1:length(latc)
        ts = filtc[:, i, j]
        xcc[:, i, j] = xcorr(tavgc_filt, ts)
        regc[i, j] = normalequation(1:length(ts), ts)
    end
end

for i in 1:length(lone)
    for j in 1:length(late)
        ts = filte[:, i, j]
        xce[:, i, j] = xcorr(tavge_filt, ts)
        rege[i, j] = normalequation(1:length(ts), ts)
    end
end

plot_corr(xcc, latc, lonc, timec, "images/xcc"*savetype)
plot_corr(xce, late, lone, timee, "images/xce"*savetype)

plot_reg(regc, latc, lonc, "images/regc"*savetype)
plot_reg(rege, late, lone, "images/rege"*savetype)


#extreme value anlaysis
quant = quantile(tavgc, [0.25, 0.5, 0.75])
thresh1 = (quant[3] - quant[1])* 1 + quant[3]
thresh2 = quant[1] - (quant[3] - quant[1])* 1
ind1 = findall(>(thresh1), tavgc)
ind2 = findall(<(thresh2), tavgc)
extreme1 = temp_mean(anomc[:,:, ind1])
extreme2 = temp_mean(anomc[:,:, ind2])

#Pot extreme values
contplot(latc,lonc,extreme2, "images/ext_anomc_low"*savetype)
contplot(latc, lonc, extreme1, "images/ext_anomc_high"*savetype)


#extreme value anlaysis
quant = quantile(tavge, [0.25, 0.5, 0.75])
thresh1 = (quant[3] - quant[1])* 1 + quant[3]
thresh2 = quant[1] - (quant[3] - quant[1])* 1
ind1 = findall(>(thresh1), tavge)
ind2 = findall(<(thresh2), tavge)
extreme1 = temp_mean(anome[:,:, ind1])
extreme2 = temp_mean(anome[:,:, ind2])

#Pot extreme values
contplot(late,lone,extreme2, "images/ext_anome_low"*savetype)
contplot(late, lone, extreme1, "images/ext_anome_high"*savetype)
