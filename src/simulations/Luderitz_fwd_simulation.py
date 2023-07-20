import netCDF4 as nc
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4, AdvectionRK4_3D
from parcels.tools.statuscodes import ErrorCode
from datetime import timedelta, datetime, date
import numpy as np
import sys
import parcels
import kernels.utilities as util

print(parcels.__version__)

#region arguments
args = sys.argv
assert len(args) == 5
start_year = np.int32(args[1])
start_mon = np.int32(args[2])
mon_name = date(1900, start_mon, 1).strftime('%b')
start_day = 1
simulation_dt=100

release_depth = np.float32(args[3])
if release_depth < 1:
    min_ind, max_ind = 0, 1
elif release_depth == 100.0:
    min_ind, max_ind = 29, 30
else:
    raise ValueError('Depth indices have not been setup.')

rk_mode = args[4]
#endregion

data_path = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

simulation_start = datetime(start_year, start_mon, start_day, 12, 0, 0)
days=[simulation_start+timedelta(days=i) for i in range(simulation_dt+1)]

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

#'lon': range(605,1903)= 60W,21E  'lat': range(0,1877) = 78S- 0Eq, (1177,1877)= 40S to Eq

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, indices={'depth': [min_ind, max_ind], 'lon': range(605,1903), 'lat': range(0,1877)}, chunksize='auto') 

u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, chunksize='auto')

coords = np.load('/nethome/manra003/analysis/dispersion/Luderitz_cell_117x117.npz')

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

if rk_mode == '2D':
    kernels= pset.Kernel(AdvectionRK4)
else:
    kernels= pset.Kernel(AdvectionRK4_3D) + util.PreventThroughSurfaceError

pset.execute(kernels,                
             runtime=timedelta(days=simulation_dt),
             dt= 300,                       
             output_file=output_file,
             recovery={ErrorCode.ErrorOutOfBounds: util.delete_particle})

output_file.close()

