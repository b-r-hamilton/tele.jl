module tele
using NCDatasets, Dates, PyCall, Revise, PyPlot, Statistics, NaNMath

export read_netcdf, contplot, anomaly, batch_open, character, filtseriesplot, fftplot

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
   function read_netcdf(path, varname, bbox, months)
   Read in NetCDF file and spatially subset into square bounding box and temporally subset into months of interest
# Arguments#
- `path`: string, path to netcdf file
- `varname`: string, name of variable you want to read
- `bbox`: geographic bounding box in format [latmin, latmax, lonmin, lonmax]
- `months`: array of months of interest (i.e. [7,8,9] means [july, aug, sept])
# Output
- `lat`: subsetted lat
- `lon`: subsetted lon
- `time`: subsetted time
- `data`: subsetted data
"""
function read_netcdf(path, varname, bbox)
   nc = NCDataset(path, "r")
   lat, lon, time = nc_dims(nc)
   #find lat lon box indices
   lat_flag = lat[begin] > lat[end] ? true : false
   a1 = searchsortedfirst(lat, bbox[1], rev = lat_flag)
   a2 = searchsortedfirst(lat, bbox[2], rev = lat_flag)
   o1 = searchsortedfirst(lon, bbox[3])
   o2 = searchsortedfirst(lon, bbox[4])
   println([a1,a2,o1,o2])

   subdata = nc[varname][minimum([o1,o2]):maximum([o1,o2]), minimum([a1, a2]):maximum([a1, a2]), :]
   na = nc[varname].attrib["missing_value"]
   replace!(subdata,na => 0)
   return lat[minimum([a1, a2]):maximum([a1, a2])], lon[minimum([o1,o2]):maximum([o1,o2])], time, subdata

end

function batch_open(paths, varname, bbox)
   times = []
   lat = []
   lon = []
   data = []
   for i in paths
      a, o, t, d = read_netcdf(i, varname, bbox)
      lat = a
      lon = o
      append!(times, t)
      if paths[1] == i
         data = d
      else
         data = cat(data, d, dims = 3)
      end
   return lat, lon, times, data
   end

end
#this is hacky!
function load_python()
   copy!(ccrs, pyimport("cartopy.crs"))
   copy!(cfeature, pyimport("cartopy.feature"))
   copy!(mpatches, pyimport("matplotlib.patches"))
   copy!(ticker, pyimport("matplotlib.ticker"))
end

function contplot(lat,lon, data,s)
   load_python()
   figure()
   ax = subplot(projection = ccrs.SouthPolarStereo())
   ax.gridlines()
   ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
   ax.coastlines()
   cf = contourf(lon,lat,transpose(data), cmap = "coolwarm", transform=ccrs.PlateCarree())
   colorbar(cf)
   savefig(s)
   close()
end

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
      marray[:, :, i] .= mean(data[:,:,t_ind], dims = 3)[:,:,1]
      sarray[:,:, i] .= std(data[:,:,t_ind], dims = 3)[:,:,1]
   end

   for i in range(1, stop = length(time), step = 1)
      m_index = Dates.month(time[i])
      anomaly[:,:,i] = (data[:,:,i] .- marray[:,:, m_index]) ./ sarray[:,:,m_index]
   end

   replace!(anomaly, NaN=>0)
   return anomaly

end

function character(series)
    min = NaNMath.minimum(series)
    max = NaNMath.maximum(series)
    mean = NaNMath.mean(series)
    median = NaNMath.median(series)
    std = NaNMath.std(series)
    dict = Dict("min"=>min, "max"=> max, "mean"=>mean, "median"=>median, "std"=>std)
    return dict
end

function filtseriesplot(series1, series1filt, series1time, series2, series2filt, series2time, t, s)
    figure(figsize = (10, 4))
    subplot(1,2,1)
    plot(series1time, series1)
    plot(series1time[2:end], series1filt[2:end])
    x_ind1 = range(1, stop = length(series1time), step = 1)
    series1filt_slope = normalequation(x_ind1, series1filt)
    series1filt_ne = series1filt_slope*x_ind1
    plot(series1time, series1filt_ne, color = "r")
    text(x = Date(1700,1,1), y = 3, s = "slope = " *string(round(series1filt_slope, digits =6)))
    legend(["CESM", "CESM smoothed", "CESM smoothed trend"])
    xlabel("time")
    ylabel("SST Anomaly [°C]")

    subplot(1,2,2)
    plot(series2time, series2)
    plot(series2time[2:end], series2filt[2:end])
    x_ind2 = range(1, stop = length(series2time[20:end]), step = 1)
    series2filt_slope = normalequation(x_ind2, series2filt[20:end])
    series2filt_ne = series2filt_slope*x_ind2
    plot(series2time[20:end], series2filt_ne, color = "r")
    text(x = Date(1880,1,1), y = 3, s = "slope = " *string(round(series2filt_slope, digits =6)))
    xlabel("time")
    legend(["ERSST", "ERSST smoothed", "ERSST smoothed trend"])
    suptitle(t)
    savefig(s)
    close()
end

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
   xt = range(ffte_freq[begin + 1], stop = ffte_freq[halfway], step = 2*10^(-10) * 100)
   ax.set_xticks(xt)
   ax.set_xlabel("Frequency [Hz]")
   ax.set_ylabel("Intensity")

   #ax.set_ylim([0, 750])
   ax.set_xlim(ffte_freq[begin+1], ffte_freq[halfway])
   ax.grid()
   title(t)
   savefig(s)
   close()

end

#linear regression with least square cost function
function normalequation(x,y)
   val = inv(transpose(x)*x)*transpose(x)*y
   return val
end

end
