module tele
using NCDatasets, Dates, PyCall, Revise, PyPlot, Statistics

export read_netcdf, contplot, anomaly, batch_open

const ccrs = PyNULL()
const cfeature = PyNULL()
const mpatches = PyNULL()

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
end

function contplot(lat,lon, data)
   load_python()
   figure()
   ax = subplot(projection = ccrs.SouthPolarStereo())
   ax.gridlines()
   ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
   ax.coastlines()
   cf = contourf(lon,lat,transpose(data), cmap = "coolwarm", transform=ccrs.PlateCarree())
   colorbar(cf)
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

end
