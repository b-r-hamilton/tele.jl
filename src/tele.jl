module tele
using NCDatasets, Dates, PyCall, Revise, PyPlot, Statistics, NaNMath

export read_netcdf, contplot, anomaly, batch_open, character, filtseriesplot, fftplot, yearly_mean, spatial_mean, temp_mean

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
   replace!(subdata,na => NaN)
   replace!(subdata,missing=> NaN)
   subdata = Float64.(subdata)
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
   ax.set_extent([-180, 0, -90, -50], ccrs.PlateCarree())
   gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=true,
                     linewidth=2, color="gray", alpha=0.5, linestyle="--")
   gl.ylabels_right = false
   ax.coastlines()
   max = NaNMath.maximum(vcat(data...))
   min = NaNMath.minimum(vcat(data...))
   println(min, max)
   val = maximum([abs(max), abs(min)])

   cf = contourf(lon,lat,transpose(data), cmap = "coolwarm", levels = range(-val, stop = val, length = 10), transform=ccrs.PlateCarree())
   cb = colorbar(cf)
   cb.set_label("SST [°C]")
   tight_layout()
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
      marray[:, :, i] .= temp_mean(data[:,:,t_ind])
      sarray[:,:, i] .= temp_std(data[:,:,t_ind])
   end

   for i in range(1, stop = length(time), step = 1)
      m_index = Dates.month(time[i])
      anomaly[:,:,i] = (data[:,:,i] .- marray[:,:, m_index]) ./ sarray[:,:,m_index]
   end

   return anomaly

end

#assume a continuous series of monthly values
function yearly_mean(time, vals)
   initial_month = Dates.month(time[1])
   initial_year = Dates.year(time[1])
   stop_month = Dates.month(time[end])
   stop_year = Dates.year(time[end])
   years = range(initial_year, stop = stop_year, step = 1)
   ann_avg = zeros(length(years))
   ann_avg[1] = mean(vals[1:13-initial_month])
   ann_avg[end] = mean(vals[end - stop_month+1:end])

   for i in range(2, length(ann_avg)-1, step = 1)
      start_index = 12-initial_month + (i-1)*12
      stop_index = start_index + 12
      ann_avg[i] = mean(vals[start_index:stop_index])
   end

   return ann_avg, years
end

function character(series)
    min = NaNMath.minimum(series)
    max = NaNMath.maximum(series)
    mean = NaNMath.mean(series)
    median = NaNMath.median(series)
    std = NaNMath.std(series)
    dict = Dict(:min=>min, :max=> max, :mean=>mean, :median=>median, :std=>std)
    return dict
end

function filtseriesplot(series1, series1filt, series1time, series2, series2filt, series2time, t, s)
    figure(figsize = (10, 4))
    subplot(1,2,1)
    plot(series1time, series1)
    plot(series1time[2:end], series1filt[2:end])
    x_ind1 = range(1, stop = length(series1time), step = 1)
    s1_avg = mean(series1filt)
    series1filt_slope = normalequation(x_ind1, series1filt .- s1_avg)
    series1filt_ne = series1filt_slope*x_ind1 .+ s1_avg
    plot(series1time, series1filt_ne, color = "r")
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
    plot(series2time[50:end], series2filt[50:end])
    x_ind2 = range(1, stop = length(series2time[50:end]), step = 1)
    s2_avg = mean(series2filt[50:end])
    print(s2_avg)
    series2filt_slope = normalequation(x_ind2, series2filt[50:end] .- s2_avg)
    series2filt_ne = series2filt_slope*x_ind2 .+ s2_avg
    plot(series2time[50:end], series2filt_ne, color = "r")
    xpos = series2time[1000]
    ypos = s2_avg
    xlabel("time")
    legend(["ERSST", "ERSST smoothed", "ERSST smoothed trend: " * string(round(series2filt_slope, digits =6))])
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

#take a 3d array (lat, lon, time) and compute the spatial mean IGNORING NAN VALUES
#JULIA IS STUPID AND WONT DO THIS WITH "SKIPMISSING" FOR A 3D ARRAY
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

function spatial_mean(array)
   avg = similar(array)[1,1,:]
   for i in range(1, stop = size(array)[3], step = 1)
      ts = array[:,:, i]
      ts = vcat(ts...)#flatten
      avg[i] = NaNMath.mean(ts)
   end
   return avg
end

function spatial_std(array)
   avg = similar(array)[1,1,:]
   for i in range(1, stop = size(array)[3], step = 1)
      ts = array[:,:, i]
      ts = vcat(ts...)#flatten
      avg[i] = NaNMath.std(ts)
   end
   return avg
end
end
