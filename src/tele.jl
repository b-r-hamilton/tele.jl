module tele
using NCDatasets, Dates, PyCall

export read_netcdf

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
   return lat, lon, time
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

   #find indices of months of interest
   bool_t_ind = findfirst.(isequal.(m), (months,)).== nothing
   t_ind = findall(x->x==0, bool_t_ind)

   nc = NCDataset(path, "r")
   subdata = nc[varname][o1:o2, a1:a2, t_ind]

   subdata = nc[varname][minimum([o1,o2]):maximum([o1,o2]), minimum([a1, a2]):maximum([a1, a2]), t_ind]
   na = nc[varname].attrib["missing_value"]
   replace!(subdata,na => 0)
   return lat[minimum([a1, a2]):maximum([a1, a2])], lon[minimum([o1,o2]):maximum([o1,o2])], time[t_ind], subdata

end

#this is hacky!
function load_python()
   copy!(ccrs, pyimport("cartopy.crs"))
   copy!(cfeature, pyimport("cartopy.feature"))
   copy!(mpatches, pyimport("matplotlib.patches"))

   return ccrs, cfeature, mpatches
end

end
