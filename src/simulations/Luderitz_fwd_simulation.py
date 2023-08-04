import netCDF4 as nc
from parcels import Field, FieldSet, ParticleSet, JITParticle, AdvectionRK4, AdvectionRK4_3D
from parcels.tools.statuscodes import ErrorCode
from datetime import timedelta, datetime, date
import numpy as np
import sys
import parcels
import kernels.utilities as util
import pandas as pd
from kernels.plankton import ZooplanktonDrift
import dask
import xarray as xr

dask.config.set({"array.slicing.split_large_chunks": False})
print(parcels.__version__)

# region arguments
args = sys.argv
assert len(args) == 5
start_year = np.int32(args[1])
start_mon = np.int32(args[2])
mon_name = date(1900, start_mon, 1).strftime('%b')
release_depth = np.float32(args[3])
rk_mode = args[4]

start_day = 1
simulation_dt = 100

if rk_mode != '2D':
    min_ind, max_ind = 0, 50
else:
    if release_depth <= 1:
        min_ind, max_ind = 0, 2
    elif release_depth == 100.0:
        min_ind, max_ind = 29, 31
    else:
        raise ValueError('Depth indices have not been setup.')
# endregion

data_path = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/'
project_data_path = '/nethome/manra003/atlanteco_plankton_plastic_nutrients_dispersion/data/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

simulation_start = datetime(start_year, start_mon, start_day, 0, 0, 0) 
#since we want the release to start at the same depth = 1 m, one day before the start is also added
days= [simulation_start-timedelta(days=1)]+[simulation_start+timedelta(days=i) for i in range(simulation_dt+1)]

ufiles = [data_path + 'ROMEO.01_1d_uo_{0}{1}{2}_grid_U.nc'.format(d.strftime("%Y"),d.strftime("%m"),d.strftime("%d")) for d in days]
vfiles = [data_path + 'ROMEO.01_1d_vo_{0}{1}{2}_grid_V.nc'.format(d.strftime("%Y"),d.strftime("%m"),d.strftime("%d")) for d in days]
wfiles = [data_path + 'ROMEO.01_1d_wo_{0}{1}{2}.nc'.format(d.strftime("%Y"),d.strftime("%m"),d.strftime("%d")) for d in days]

def define_2Dvariables():
    filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': ufiles[0], 'data': ufiles},
                'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': ufiles[0], 'data': vfiles}}

    variables = {'U': 'uo',
                'V': 'vo'}

    dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthu', 'time': 'time_counter'}
    return filenames, variables, dimensions


def define_3Dvariables():
    filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
                'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles},
                'W': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles}}

    variables = {'U': 'uo',
                'V': 'vo',
                'W': 'wo'}

    dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}
    return filenames, variables, dimensions

if rk_mode == '2D':
    print('load 2D fields')
    filenames, variables, dimensions = define_2Dvariables()
else:
    print('load 3D fields')
    filenames, variables, dimensions = define_3Dvariables()

#'lon': range(605,1903)= 60W,21E  'lat': range(0,1877) = 78S- 0Eq, (1177,1877)= 40S to Eq 'depth': range(min_ind, max_ind), 
fieldset = FieldSet.from_nemo(filenames, variables, dimensions, indices={'depth': range(min_ind, max_ind), 'lon': range(605,1903), 'lat': range(1177,1877)}, chunksize=False) 

if rk_mode == 'DVM':

    fieldset.add_constant('migration_speed', 0.035) # in m/s
    fieldset.add_constant('min_depth', 0.035) # in m 
    fieldset.add_constant('max_depth', 150) # in m

    # total number of seconds in that "day" from which the data is available example in GLOB16: 12:00:00.
    time_zero_totalseconds = (datetime.strptime(str(fieldset.U.grid.time_origin)[:19],'%Y-%m-%dT%H:%M:%S') - simulation_start).total_seconds()
    
    simulation_start.hour * 60 * 60 + simulation_start.minute * 60 + simulation_start.second
    fieldset.add_constant('start_time', time_zero_totalseconds)

    years= np.unique([d.year for d in days])

    # load a non leap-year values for sunrise and sunset
    sunrise_nc = nc.Dataset(project_data_path + 'SunriseTime_2x2_1d_2015.nc','r')
    sunrise_da = xr.DataArray(sunrise_nc['sunrise'][::], coords={"time": sunrise_nc['time'][:], "lat": sunrise_nc['lat'][:], "lon": sunrise_nc['lon'][:]})

    sunset_nc = nc.Dataset(project_data_path + 'SunsetTime_2x2_1d_2015.nc','r')
    sunset_da = xr.DataArray(sunset_nc['sunset'][::], coords={"time": sunset_nc['time'][:], "lat": sunset_nc['lat'][:], "lon": sunset_nc['lon'][:]})
    
    # for unique years in the data files, get the datetimes for all days and remove the leap year dates. 
    # It is okay to have missing values for a day in Feb, interpolation should work fine. It is just o have consistent values for all days
    dates = pd.date_range(date(years[0],1,1), date(years[0],12,31))
    sunrise_da['time']= dates[~dates.is_leap_year]
    sunset_da['time']= dates[~dates.is_leap_year]
    
    sunrise_nc.close()
    sunset_nc.close()
    sunrise_copy=sunrise_da.copy()
    sunset_copy=sunset_da.copy()

    for y in years[1:]:
        dates = pd.date_range(date(y, 1, 1), date(y, 12, 31))
        sunrise_copy['time']= dates[~dates.is_leap_year]
        sunrise_da = xr.combine_by_coords([sunrise_da, sunrise_copy], coords=['time'])

        sunset_copy['time']= dates[~dates.is_leap_year]
        sunset_da = xr.combine_by_coords([sunset_da, sunset_copy], coords=['time'])
    
    fieldset.add_field(Field.from_xarray(sunrise_da,
                                        'Sunrise', 
                                        dimensions={'lon':'lon','lat':'lat', 'time':'time'},
                                        mesh='spherical'))

    fieldset.add_field(Field.from_xarray(sunset_da,
                                        'Sunset', 
                                        dimensions={'lon':'lon','lat':'lat', 'time':'time'},
                                        mesh='spherical'))

u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

coords = pd.read_csv(project_data_path + 'Benguela_release_points_961x641_grid_015625.csv')

if release_depth == 0:
    depth_arg = None
else:
    depth_arg = [release_depth for i in range(len(coords['Longitude']))]

pset = ParticleSet.from_list(fieldset=fieldset, 
                             pclass=JITParticle,
                             lon=coords['Longitude'],
                             lat=coords['Latitude'],
                             depth=depth_arg,
                             time=simulation_start)
pset.populate_indices()                            
output_file = pset.ParticleFile(name="/nethome/manra003/analysis/dispersion/simulations/Fwd_{0}_Jul2023_BenguelaUpwR_117x117_{1}{2}_{3}z_{4}days.zarr".format(rk_mode, mon_name, start_year, int(release_depth), simulation_dt), 
                                outputdt=timedelta(days=1))

if rk_mode == '3D':
    kernels = pset.Kernel(AdvectionRK4_3D) + pset.Kernel(util.PreventThroughSurfaceError)
elif rk_mode == 'DVM':
    kernels = pset.Kernel(AdvectionRK4) + pset.Kernel(ZooplanktonDrift) 
else:
    kernels= pset.Kernel(AdvectionRK4)

pset.execute(kernels,                
             runtime=timedelta(days=simulation_dt),
             dt= timedelta(minutes=10),                       
             output_file=output_file,
             recovery={ErrorCode.ErrorOutOfBounds: util.delete_particle})

output_file.close()