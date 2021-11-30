module tele
using NCDatasets, Dates, PyCall, Revise, PyPlot, Statistics, NaNMath

export read_netcdf, contplot, anomaly, batch_open, character,filtseriesplot,
normalequation, fftplot, yearly_mean, spatial_mean, temp_mean, plot_eofs, plot_corr, plot_reg

const ccrs = PyNULL()
const cfeature = PyNULL()
const mpatches = PyNULL()
const ticker = PyNULL()

"""
   function nc_dims(nc)
   Get all important netcdf dims assume CF-conventions
# Arguments
- `nc`: NCDataset
# Output
- `lat`: all lat values
- `lon`: all lon values
- `time`: time values (already converted to Datetime)
"""
function nc_dims(nc)
   lat = getindex(nc, "lat")
   lon = getindex(nc, "lon")
   time = getindex(nc, "time")
   if length(size(lat)) > 1
      a = lat[1,:]
      o = lon[:,1]
   else
      a = lat
      o = lon
   end

   return a, o, time
end

"""
   function read_netcdf(path, varname, bbox)
   Read in NetCDF file and spatially subset into square bounding box
# Arguments#
- `path`: string, path to netcdf file
- `varname`: string, name of variable you want to read
- `bbox`: geographic bounding box in format [latmin, latmax, lonmin, lonmax]
# Outputlat_diff = pi*r_e / 180 #km/°latlat_diflat_diff = pi*r_e / 180 #km/°lat
f = pi*r_e / 180 #km/°lat


- `lat`: subsetted lat
- `lon`: subsetted lon
- `time`: subsetted time
- `data`: subsetted data
"""
function read_netcdf(path, varname, bbox)
   nc = NCDataset(path, "r")
   lat, lon, time = nc_dims(nc)
   #find lat lon box indices
   lat_flag = lat[begin] > lat[end] ? true : false #check if lat is reversed
   a1 = searchsortedfirst(lat, bbox[1], rev = lat_flag)
   a2 = searchsortedfirst(lat, bbox[2], rev = lat_flag)
   o1 = searchsortedfirst(lon, bbox[3])
   o2 = searchsortedfirst(lon, bbox[4])
   println([a1,a2,o1,o2])

   #is our longitude array continuous - not true for CESM
   if maximum(abs.(diff(lon))) > 100
      ind = findall(>(100), abs.(diff(lon)))[1]
      lon1 = lon[begin:ind-1]
      lon2 = lon[ind:end]


      o1 = searchsortedfirst(lon1, bbox[3])

      bool = zeros(2)
      if o1 == 1
         o1 = searchsortedfirst(lon2, bbox[3])
         bool[1] = 1
      end

      o2 = searchsortedfirst(lon1, bbox[4])
      if o2 == 1
         o2 = searchsortedfirst(lon2, bbox[4])
         bool[2] = 1
      end

      println("bool is " *string(bool))
      if bool == zeros(2)
         println("first bool : "* string([o1,o2]))
         lon_sub = lon1[minimum([o1,o2]):maximum([o1,o2])]
         subdata = nc[varname][minimum([o1,o2]):maximum([o1,o2]), minimum([a1, a2]):maximum([a1, a2]), :]

      elseif bool == [1,1]
         lon_sub = lon1[minimum([o1,o2]):maximum([o1,o2])]
         o1 = length(lon1) + o1
         o2 = length(lon1) + o2
         subdata = nc[varname][minimum([o1,o2]):maximum([o1,o2]), minimum([a1, a2]):maximum([a1, a2]), :]
      #FILLOWING CODE IS NOT TESTED - this is where I'd put my good code (if I had it)
      # elseif bool == [0, 1]
      #    lon_sub = vcat(lon1[o1:end], lon2[begin:o2])
      #    subdata1 = nc[varname][o1:length(lon1), minimum([a1, a2]):maximum([a1, a2]), :]
      #    subdata2 = nc[varname][length(lon1):o2+length(lon1), minimum([a1, a2]):maximum([a1, a2]), :]
      #    subdata = cat(subdata1, subdata2, dims = 2)

      #this works! (amazing)
      elseif bool == [1, 0]
         lon_sub = vcat(lon2[o1:end], lon1[begin:o2])
         subdata1 = nc[varname][length(lon1) + o1:end, minimum([a1, a2]):maximum([a1, a2]), :]
         subdata2 = nc[varname][begin:o2, minimum([a1, a2]):maximum([a1, a2]), :]
         subdata = cat(subdata1, subdata2, dims = 1)
      end

   else
      subdata = nc[varname][minimum([o1,o2]):maximum([o1,o2]), minimum([a1, a2]):maximum([a1, a2]), :]
      lon_sub = lon[minimum([o1,o2]):maximum([o1,o2])]
   end

   na = nc[varname].attrib["missing_value"]
   replace!(subdata,na => NaN)
   replace!(subdata,missing=> NaN)
   subdata = Float64.(subdata)
   return lat[minimum([a1, a2]):maximum([a1, a2])], lon_sub, time, subdata

end

#this is hacky!
function load_python()
   copy!(ccrs, pyimport("cartopy.crs"))
   copy!(cfeature, pyimport("cartopy.feature"))
   copy!(mpatches, pyimport("matplotlib.patches"))
   copy!(ticker, pyimport("matplotlib.ticker"))
end

"""
   function contplot(lat, lon, data, s)
   Make a contour plot of data (assumes SouthPolarStereo view) and save
# Arguments#
- `lat`
- `lon`
- `data`: data in (lat x lon x time) format
- `s`: file path to save to
"""
function contplot(lat,lon, data,s)
   load_python()
   figure()
   ax = subplot(projection = ccrs.PlateCarree())
   ax.gridlines()
   #ax.set_extent([-180, 0, 50, 90], ccrs.PlateCarree())
   gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=true,
                     linewidth=2, color="gray", alpha=0.5, linestyle="--")
   gl.right_labels = false
   ax.coastlines()
   max = NaNMath.maximum(vcat(data...))
   min = NaNMath.minimum(vcat(data...))
   println(min, max)
   val = maximum([abs(max), abs(min)])

   cf = contourf(lon,lat,transpose(data), cmap = "coolwarm", levels = range(-val, stop = val, length = 10), transform=ccrs.PlateCarree())
   cb = colorbar(cf, fraction = 0.015)
   cb.set_label("SST [σ]")
   tight_layout()
   savefig(s,bbox_inches="tight")
   close()
end

"""
   function anomaly(data, time)
   Remove seasonality by computing the anomaly (stdev away from mean) w.r.t. behavior in month
# Arguments#
- `data`: array of data in (lat x lon x time) format
- `time`: array of DateTime
# Output
- `anomaly`: array of same dimensions of `data` where each value is the anomaly in stdev
"""
function anomaly(data, time)
   anomaly = similar(data)

   #make list of months
   m = []
   for i in time
      push!(m, Dates.month(i))
   end

   #compute average for each month
   marray = Array{Float32}(undef, size(data)[1], size(data)[2],12)
   sarray = Array{Float32}(undef, size(data)[1], size(data)[2],12)
   for i in range(1, stop = 12, step = 1) #iterate over each month, collecting all vals

      #find indices of months of interest
      bool_t_ind = findfirst.(isequal.(m), (i,)).== nothing
      t_ind = findall(x->x==0, bool_t_ind)
      marray[:, :, i] .= temp_mean(data[:,:,t_ind])
      sarray[:,:, i] .= temp_std(data[:,:,t_ind])
   end

   for i in range(1, stop = length(time), step = 1)
      m_index = Dates.month(time[i])
      anomaly[:,:,i] = (data[:,:,i] .- marray[:,:, m_index]) ./ sarray[:,:,m_index]
   end

   return anomaly

end

"""
   function yearly_mean(time, vals)
   Compute annual average of 1d array (assumes data is ordered AND continuous)
# Arguments#
- `time`: array of time values (DateTime)
- `vals`: array of associated values
# Output
- `ann_avg`: annual average of `vals`
- `years`: list of years associated with `ann_avg`
"""
function yearly_mean(time, vals)
   initial_month = Dates.month(time[1])
   initial_year = Dates.year(time[1])
   stop_month = Dates.month(time[end])
   stop_year = Dates.year(time[end])

   years = range(initial_year, stop = stop_year, step = 1) #array of all possible years

   ann_avg = zeros(length(years))
   #compute first and last values, because we might not have 12 months in first or last year
   ann_avg[1] = mean(vals[1:13-initial_month])
   ann_avg[end] = mean(vals[end - stop_month+1:end])

   #iterate through, computing the average for every 12 months, starting from end of year 1
   for i in range(2, length(ann_avg)-1, step = 1)
      start_index = 12-initial_month + (i-1)*12
      stop_index = start_index + 12
      ann_avg[i] = mean(vals[start_index:stop_index])
   end

   return ann_avg, years
end

"""
   function character(series)
   Compute statistics for a series
# Arguments#
- `series`: array of floats, can have NaN values (will be ignored)
# Output
- `dict`: dictionary of statistics
"""
function character(series)
    min = NaNMath.minimum(series)
    max = NaNMath.maximum(series)
    mean = NaNMath.mean(series)
    median = NaNMath.median(series)
    std = NaNMath.std(series)
    dict = Dict(:min=>min, :max=> max, :mean=>mean, :median=>median, :std=>std)
    return dict
end

"""
   function filtseriesplot(series)
      This is a pretty specific function - don't need to go into everything
      Makes plot with 2 subplots, each has raw data, filtered data, and a trendline
"""
function filtseriesplot(series1, series1filt, series1time, series2, series2filt, series2time, t, s)
   figure(figsize = (10, 4))
   subplot(1,2,1)
   shift = 50
   plot(series1time, series1)
   plot(series1time[shift:end], series1filt[shift:end])
   x_ind1 = range(1, stop = length(series1time[shift:end]), step = 1)
   s1_avg = mean(series1filt[shift:end])
   series1filt_slope = normalequation(x_ind1, series1filt[shift:end] .- s1_avg)
   series1filt_ne = series1filt_slope*x_ind1 .+ s1_avg
   plot(series1time[shift:end], series1filt_ne, color = "r")
   xpos = series1time[1000]
   ypos = s1_avg
   legend(["CESM", "CESM smoothed", "CESM smoothed trend: "* string(round(series1filt_slope, digits =6))])
   xlabel("time")

   if s == "images/anom.pdf"
      ylabel("SST Anomaly [σ]")
   else
      ylabel("SST [°C]")
   end

   subplot(1,2,2)
   plot(series2time, series2)
   shift = 50
   plot(series2time[shift:end], series2filt[shift:end])
   x_ind2 = range(1, stop = length(series2time[shift:end]), step = 1)
   s2_avg = mean(series2filt[shift:end])
   print(s2_avg)
   series2filt_slope = normalequation(x_ind2, series2filt[shift:end] .- s2_avg)
   series2filt_ne = series2filt_slope*x_ind2 .+ s2_avg
   plot(series2time[shift:end], series2filt_ne, color = "r")
   xpos = series2time[1000]
   ypos = s2_avg
   xlabel("time")
   legend(["ERSST", "ERSST smoothed", "ERSST smoothed trend: " * string(round(series2filt_slope, digits =6))])
   suptitle(t)
   savefig(s)
   close()
end

"""
   function fftplot(ffte_freq, ffte, s, t)
   Plot frequency vs intensity for FFTW resuls
# Arguments#
- `ffte_freq`: output of FFTW.fft_freq
- `ffte`: ouptu of FFTW.fft
- `s`: save strng
- `t`:title string
"""
function fftplot(ffte_freq, ffte, s, t)
   load_python()
   figure()
   ax = subplot()
   N = length(ffte_freq)
   halfway = 0.0
   if N % 2 == 1
      halfway = convert(Int, N/2 - 0.5 )
   else
      halfway = convert(Int, N/2)
   end
   replace!(ffte_freq, 0.0 => 1.0)
   print(halfway)

   if typeof(ffte[1]) in [Vector{ComplexF32}, Vector{ComplexF64}]
      for i in ffte
         ax.plot(ffte_freq[begin+1:halfway], abs.(i)[begin+1:halfway], linewidth = 2, alpha = 0.5)
      end
      legend(["Raw", "Anomaly"])
   else
      ax.plot(ffte_freq[begin+1:halfway], abs.(ffte)[begin+1:halfway])
   end
   xt = range(ffte_freq[begin + 1], stop = ffte_freq[halfway], step = 2*10^(-10) * 10)
   ax.set_xticks(xt)
   ax.set_xlabel("Frequency [Hz]")
   ax.set_ylabel("Intensity")

   #ax.set_ylim([0, 750])
   ax.set_xlim(ffte_freq[begin+1], 2 * 0.1*10^(-7))
   ax.set_ylim(0, 500)
   ax.grid()
   title(t)
   savefig(s)
   close()

end

"""
   function normalequation(x,y)
   Linear regression with least squares fit
   Assumes y-intercept is already 0 (mean removed )
# Arguments#
- `x`: dependent variable
- `y`: independent variable
# Output
- `val`: slope of line of best fit
"""
function normalequation(x,y)
   val = inv(transpose(x)*x)*transpose(x)*y
   return val
end

"""
   function tempmean(array)
   Computes temporal mean of 3d array, IGNORING NANs
# Arguments#
- `array`: 3d float array in format (lat x lon x time)
# Output
- `avg`: 2d array in format (lat x lon), where every value is temporal average
"""
function temp_mean(array)
   avg = similar(array)[:,:,1]

   for i in range(1, stop = size(avg)[1], step = 1)
      for j in range(1, stop = size(avg)[2], step = 1)
         ts = array[i, j, :]
         avg[i, j] = NaNMath.mean(ts)
      end
   end
   return avg
end

"""
   function tempstd(array)
   Computes temporal standard deviation of 3d array, IGNORING NANs
# Arguments#
- `array`: 3d float array in format (lat x lon x time)
# Output
- `avg`: 2d array in format (lat x lon), where every value is temporal standard deviation
"""
function temp_std(array)
   avg = similar(array)[:,:,1]

   for i in range(1, stop = size(avg)[1], step = 1)
      for j in range(1, stop = size(avg)[2], step = 1)
         ts = array[i, j, :]
         avg[i, j] = NaNMath.std(ts)
      end
   end
   return avg
end

"""
   function spatial_mean(array)
   Computes spatial average of 3d array, IGNORING NANs
# Arguments#
- `array`: 3d float array in format (lat x lon x time)
# Output
- `avg`: 1d array in format (time), where every value is spatial average at that time
"""
function spatial_mean(array)
   avg = similar(array)[1,1,:]
   for i in range(1, stop = size(array)[3], step = 1)
      ts = array[:,:, i]
      ts = vcat(ts...)#flatten
      avg[i] = NaNMath.mean(ts)
   end
   return avg
end

"""
   function spatial_std(array)
   Computes spatial average of 3d array, IGNORING NANs
# Arguments#
- `array`: 3d float array in format (lat x lon x time)
# Output
- `avg`: 1d array in format (time), where every value is spatial standard deviation at that time
"""
function spatial_std(array)
   avg = similar(array)[1,1,:]
   for i in range(1, stop = size(array)[3], step = 1)
      ts = array[:,:, i]
      ts = vcat(ts...)#flatten
      avg[i] = NaNMath.std(ts)
   end
   return avg
end

function plot_eofs(solver, lon, lat, time, s)
   load_python()
   eofs = solver.eofs(neofs=2)
   pcs = solver.pcs(npcs = 2, pcscaling = 1)
   var = solver.varianceFraction()
   vmin = NaNMath.minimum(eofs)
   vmax = NaNMath.maximum(eofs)

   vmin = abs(vmin) > abs(vmax) ? vmin : -vmax
   vmax = abs(vmax) > abs(vmin) ? vmax : -vmin

   figure(figsize = (12, 5))
   ax = subplot(2,2,1, projection = ccrs.PlateCarree())
   gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=true,
                     linewidth=2, color="gray", alpha=0.5, linestyle="--")
   gl.right_labels = false
   gl.top_labels = false

   pc = contourf(lon, lat, transpose(eofs[1, :, :]), cmap = "coolwarm", transform = ccrs.PlateCarree(), vmin = vmin, vmax = vmax)
   frac = 0.013
   cb = colorbar(pc, fraction = frac)
   ax.coastlines()
   cb.set_label("SST Anomaly [σ]")
   title("EOF1, Exp. Var = " *string(round(var[1]*100, digits = 2)) * "%")
   ax = subplot(2,2,3)
   plot(time, pcs[:, 1])
   ax.set_ylim(-2.5, 2.5)
   xlabel("time [years]")
   ax.set_xlim(Dates.DateTime(1860,1,1), time[end])
   grid()
   title("PC1")

   ax = subplot(2,2,2, projection = ccrs.PlateCarree())
   gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=true,
                     linewidth=2, color="gray", alpha=0.5, linestyle="--")
   gl.right_labels = false
   gl.top_labels = false
   gl.left_labels = false
   pc = contourf(lon, lat, transpose(eofs[2, :, :]), cmap = "coolwarm", transform = ccrs.PlateCarree(), vmin = vmin, vmax = vmax )
   cb = colorbar(pc, fraction = frac)
   ax.coastlines()
   cb.set_label("SST Anomaly [σ]")
   title("EOF2, Exp. Var = " *string(round(var[2]*100, digits = 2)) * "%")
   ax = subplot(2,2,4)
   plot(time, pcs[:, 2])
   xlabel("time [years]")
   ax.set_ylim(-2.5, 2.5)
   ax.set_xlim(Dates.DateTime(1860,1,1), time[end])
   title("PC2")
   grid()
   savefig(s,bbox_inches="tight")
   close()
end

function plot_corr(xc, lat, lon, time, s)
   load_python()
   figure()
   ax = subplot(2,1,1, projection = ccrs.PlateCarree())
   ax.coastlines()
   halfway = (size(xc)[1] - 1)/2
   halfway = convert(Int32, halfway)
   vmin = NaNMath.minimum(xc)
   vmax = NaNMath.maximum(xc)

   vmin = abs(vmin) > abs(vmax) ? vmin : -vmax
   vmax = abs(vmax) > abs(vmin) ? vmax : -vmin

   pc = contourf(lon, lat, transpose(xc[halfway, :, :]), cmap = "coolwarm", transform = ccrs.PlateCarree(), vmin = vmin, vmax = vmax)
   cb = colorbar(pc, fraction = 0.013)
   cb.set_label("Correlation")
   gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=true,
                     linewidth=2, color="gray", alpha=0.5, linestyle="--")
   gl.right_labels = false
   gl.top_labels = false
   subplot(2,1,2)
   reshaped = zeros(size(xc)[2], size(xc)[3], size(xc)[1])

   permutedims!(reshaped, xc,[2,3,1])
   mean = spatial_mean(reshaped)
   lag = -halfway:halfway
   plot(lag, mean)
   xlabel("Lag [Months]")
   ylabel("Autocorrelation")
   savefig(s,bbox_inches="tight")
   close()
end

function plot_reg(reg, lat, lon, s)
   load_python()
   figure()
   ax = subplot(projection = ccrs.PlateCarree())
   ax.coastlines()
   vmin = NaNMath.minimum(reg)
   vmax = NaNMath.maximum(reg)

   vmin = vmin < 0 ? vmin : -vmin
   vmin = abs(vmin) > abs(vmax) ? vmin : -vmax
   vmax = abs(vmax) > abs(vmin) ? vmax : -vmin

   pc = contourf(lon, lat, transpose(reg), cmap = "coolwarm", vmin = vmin, vmax = vmax)
   cb = colorbar(pc, fraction = 0.013, format="%.0e")
   cb.set_label("Regression Coefficient")
   gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=true,
                     linewidth=2, color="gray", alpha=0.5, linestyle="--")
   gl.right_labels = false
   gl.top_labels = false

   savefig(s,bbox_inches="tight")
   close()
end

end
