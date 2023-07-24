from core.algorithm import AlmanacSunrise, AlmanacSunset
from datetime import datetime, timedelta
from math import floor,ceil
import netCDF4 as nc
import numpy as np  

# inputs for this tool are duration and area on the map for which the sunrise_sunset should be calculated
# and for which grid size
start_date = datetime(2015, 1, 1)
end_date = datetime(2015, 12, 31)

min_lat = -60
max_lat = 61
min_lon = -180
max_lon = 181

grid_size = 5

lons = np.arange(min_lon, max_lon, grid_size)  
lats = np.arange(min_lat, max_lat, grid_size)  
print(lats)
print(lons)

# TODO: remove land cells from the mesh

current_date = start_date
t_days=(end_date - start_date).days + 1
increment_days=7
n_days=ceil(t_days/increment_days)

# create the dataset to create
ds_sunrise = np.zeros([n_days, len(lats), len(lons)])
ds_sunset = np.zeros([n_days, len(lats), len(lons)])
day = 0

while day < n_days:

    for j in range(len(lats)):
        for i in range(len(lons)):
            sunrise = AlmanacSunrise(lats[j], lons[i], current_date)
            sunset = AlmanacSunset(lats[j], lons[i], current_date)
            ds_sunrise[day, j, i] = sunrise
            ds_sunset[day, j, i] = sunset
    current_date += timedelta(increment_days)
    day += 1

print("sunrise, sunset calculated")

home_folder='/nethome/manra003/atlanteco_plankton_plastic_nutrients_dispersion/data/'

ds = nc.Dataset(home_folder + 'SunriseTime_5x5_7d_Atlantic_2015.nc', 'w', format='NETCDF4')
ds.description = "File to store sunrise time in UTC with grid size of 2 x 2 degree"
ds.history = "Created " + datetime.utcnow().strftime("%d/%m/%y")

ds.createDimension('lon', len(lons))
ds.createDimension('lat', len(lats))
ds.createDimension('time', n_days)

times = ds.createVariable('time', 'f8', ('time',))
lats1 = ds.createVariable('lat', 'f4', ('lat',))
lons1 = ds.createVariable('lon', 'f4', ('lon',))
sunrise_time = ds.createVariable('sunrise', 'f4', ('time', 'lat', 'lon',))

sunrise_time.units = 'UTC'
times.units = 'seconds since 2015-01-01'
times.calendar = 'standard'

dates = [start_date + n * timedelta(days=increment_days) for n in range(n_days)]
times[:] = nc.date2num(dates, units=times.units, calendar=times.calendar)
lats1[:] = lats
lons1[:] = lons
sunrise_time[::] = ds_sunrise

ds.close()

ds = nc.Dataset(home_folder + 'SunsetTime_5x5_7d_Atlantic_2015.nc', 'w', format='NETCDF4')
ds.description = "File to store sunset time in UTC with grid size of 2 x 2 degree"
ds.history = "Created " + datetime.utcnow().strftime("%d/%m/%y")

ds.createDimension('lon', len(lons))
ds.createDimension('lat', len(lats))
ds.createDimension('time', n_days)

times = ds.createVariable('time', 'f8', ('time',))
lats1 = ds.createVariable('lat', 'f4', ('lat',))
lons1 = ds.createVariable('lon', 'f4', ('lon',))
sunset_time = ds.createVariable('sunset', 'f4', ('time', 'lat', 'lon',))

sunset_time.units = 'UTC'
times.units = 'seconds since 2015-01-01'
times.calendar = 'standard'

dates = [start_date + n * timedelta(days=increment_days) for n in range(n_days)]
times[:] = nc.date2num(dates, units=times.units, calendar=times.calendar)
lats1[:] = lats
lons1[:] = lons
sunset_time[::] = ds_sunset

ds.close()

print("sunrise sunset saved in the NetCDF files")