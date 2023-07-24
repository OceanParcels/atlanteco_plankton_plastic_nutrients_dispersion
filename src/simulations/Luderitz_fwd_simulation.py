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
from parcels.tools.converters import TimeConverter
import dask

dask.config.set({"array.slicing.split_large_chunks": False})
print(parcels.__version__)

#region arguments
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
#endregion

data_path = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/'
project_data_path = '/nethome/manra003/atlanteco_plankton_plastic_nutrients_dispersion/data/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

simulation_start = datetime(start_year, start_mon, start_day, 0, 0, 0) #UTC +2
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

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, indices={'lon': range(605,1903), 'lat': range(1177,1877)}, chunksize=False) 

if rk_mode == 'DVM':
    sunrise_nc = nc.Dataset(project_data_path + 'SunriseTime_5x5_7d_Atlantic_2015.nc','r')
    sunset_nc = nc.Dataset(project_data_path + 'SunsetTime_5x5_7d_Atlantic_2015.nc','r')
        
    # total number of seconds in that "day" from which the data is available example in GLOB1: 12:00:00.
    time_zero_totalseconds = simulation_start.hour * 60 * 60 + simulation_start.minute * 60 + simulation_start.second
    fieldset.add_constant('start_time', time_zero_totalseconds)

    time_origin = TimeConverter(np.datetime64(nc.num2date(sunrise_nc['time'][0],'seconds since '+ str(simulation_start.year) + '-01-01',sunrise_nc['time'].calendar)))

    fieldset.add_field(Field("Sunrise",
                            data=sunrise_nc['sunrise'][::],
                            lon=sunrise_nc['lon'][:],
                            lat=sunrise_nc['lat'][:],
                            time=sunrise_nc['time'][:],
                            time_origin=time_origin,
                            transpose=False,
                            time_periodic=timedelta(days=365),
                            interp_method='linear'))

    fieldset.add_field(Field("Sunset",
                            data=sunset_nc['sunset'][::],
                            lon=sunset_nc['lon'][:],
                            lat=sunset_nc['lat'][:],
                            time=sunset_nc['time'][:],
                            time_origin=time_origin,
                            transpose=False,
                            time_periodic=timedelta(days=365),
                            interp_method='linear'))


u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

coords = pd.read_csv(project_data_path + 'Luderitz_cell_117x117.csv')

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
output_file = pset.ParticleFile(name="/nethome/manra003/analysis/dispersion/simulations/Fwd_{0}_Luderitz_117x117_{1}01-31_{2}_{3}z_{4}days.zarr".format(rk_mode, mon_name, start_year, int(release_depth), simulation_dt), 
                                outputdt=timedelta(days=1))

if rk_mode == '3D':
    kernels = pset.Kernel(AdvectionRK4_3D) + pset.Kernel(util.PreventThroughSurfaceError)
elif rk_mode == 'DVM':
    kernels = pset.Kernel(AdvectionRK4_3D) + pset.Kernel(ZooplanktonDrift) 
else:
    kernels= pset.Kernel(AdvectionRK4)

pset.execute(kernels,                
             runtime=timedelta(days=simulation_dt),
             dt= 300,                       
             output_file=output_file,
             recovery={ErrorCode.ErrorOutOfBounds: util.delete_particle})

output_file.close()

