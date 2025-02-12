from core.algorithm import AlmanacSunrise, AlmanacSunset
from datetime import datetime, timedelta
import netCDF4 as nc
import numpy as np  

# inputs for this tool are duration and area on the map for which the dawn_dusk should be calculated
# and for which grid size
start_date = datetime(2015, 1, 1)
end_date = datetime(2015, 12, 31)

min_lat = -60
max_lat = 61
min_lon = -180
max_lon = 181

grid_size = 2

lons = np.arange(min_lon, max_lon, grid_size)  
lats = np.arange(min_lat, max_lat, grid_size)  
print(lats)
print(lons)

sunrise_sunset_type= 'civil'

current_date = start_date
n_days = 365
# create the dataset to create
ds_dawn = np.zeros([n_days, len(lats), len(lons)])
ds_dusk = np.zeros([n_days, len(lats), len(lons)])
day = 0

while day < n_days:

    for j in range(len(lats)):
        for i in range(len(lons)):
            dawn = AlmanacSunrise(lats[j], lons[i], current_date, sunrise_sunset_type)
            dusk = AlmanacSunset(lats[j], lons[i], current_date, sunrise_sunset_type)
            ds_dawn[day, j, i] = dawn
            ds_dusk[day, j, i] = dusk
    current_date += timedelta(days=1)
    day += 1

print("dawn, dusk calculated")

# store the output of dawn and dusk in a netcdf format to enable reading by parcels kernel.

home_folder='/nethome/manra003/atlanteco_plankton_plastic_nutrients_dispersion/data/'

ds = nc.Dataset(home_folder + 'DawnTime_2x2_1d_2015.nc', 'w', format='NETCDF4')
ds.description = "File to store dawn time in UTC with grid size of 2 x 2 degree"
ds.history = "Created " + datetime.utcnow().strftime("%d/%m/%y")

ds.createDimension('lon', len(lons))
ds.createDimension('lat', len(lats))
ds.createDimension('time', n_days)

times = ds.createVariable('time', 'f8', ('time',))
lats1 = ds.createVariable('lat', 'f4', ('lat',))
lons1 = ds.createVariable('lon', 'f4', ('lon',))
dawn_time = ds.createVariable('dawn', 'f4', ('time', 'lat', 'lon'))

dawn_time.units = 'UTC'
times.units = 'seconds since 2015-01-01'
times.calendar = 'standard'

dates = [start_date + n * timedelta(days=1) for n in range(n_days)]
times[:] = nc.date2num(dates, units=times.units, calendar=times.calendar)
lats1[:] = lats
lons1[:] = lons
dawn_time[::] = ds_dawn

ds.close()

ds = nc.Dataset(home_folder + 'DuskTime_2x2_1d_2015.nc', 'w', format='NETCDF4')
ds.description = "File to store dusk time in UTC with grid size of 2 x 2 degree"
ds.history = "Created " + datetime.utcnow().strftime("%d/%m/%y")

ds.createDimension('lon', len(lons))
ds.createDimension('lat', len(lats))
ds.createDimension('time', n_days)

times = ds.createVariable('time', 'f8', ('time',))
lats1 = ds.createVariable('lat', 'f4', ('lat',))
lons1 = ds.createVariable('lon', 'f4', ('lon',))
dusk_time = ds.createVariable('dusk', 'f4', ('time', 'lat', 'lon'))

dusk_time.units = 'UTC'
times.units = 'seconds since 2015-01-01'
times.calendar = 'standard'

dates = [start_date + n * timedelta(days=1) for n in range(n_days)]
times[:] = nc.date2num(dates, units=times.units, calendar=times.calendar)
lats1[:] = lats
lons1[:] = lons
dusk_time[::] = ds_dusk

ds.close()

print("dawn dusk saved in the NetCDF files")