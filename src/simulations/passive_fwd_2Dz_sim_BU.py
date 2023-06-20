import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4
from parcels.tools.statuscodes import ErrorCode
from datetime import timedelta, datetime, date
import numpy as np
import sys

args = sys.argv
assert len(args) == 4
year = np.int32(args[1])
month_num = np.int32(args[2])
mon_name = date(1900, month_num, 1).strftime('%b')

r_depth = np.int32(args[3])
if r_depth == 0:
    min_ind, max_ind = 0, 1
elif r_depth == 100:
    min_ind, max_ind = 29, 30
else:
    raise ValueError('Depth indices have not been setup.')

data_path = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

simulation_start = datetime(year, month_num, 1, 12, 0, 0)
simulation_end = datetime(year, month_num, 31, 12, 0, 0)

ufiles =  sorted(glob(data_path + 'ROMEO.01_1d_uo_{0}{1}*_U.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))))

vfiles =  sorted(glob(data_path + 'ROMEO.01_1d_vo_{0}{1}*_V.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))))

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': ufiles[0], 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': ufiles[0], 'data': vfiles}}

variables = {'U': 'uo',
             'V': 'vo'}

dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthu', 'time': 'time_counter'}

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, indices={'depth': [min_ind, max_ind]}, chunksize=False)


u_file = nc.Dataset(ufiles[0])
ticks = u_file['time_counter'][:][0]
modeldata_start = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_start >= modeldata_start

u_file = nc.Dataset(ufiles[len(ufiles)-1])
ticks = u_file['time_counter'][:][0]
modeldata_end = datetime(1900, 1, 1) + timedelta(seconds=ticks)

assert simulation_end <= modeldata_end

fieldset = FieldSet.from_nemo(filenames, variables, dimensions, chunksize='auto')

coords = np.load('/nethome/manra003/analysis/dispersion/Benguela_0625_401x257_release_points.npz')

def delete_particle(particle, fieldset, time):
    particle.delete()

if r_depth == 0:
    depth_arg = None
else:
    depth_arg = [r_depth for i in range(len(coords['Longitude']))]

pset = ParticleSet.from_list(fieldset=fieldset, 
                             pclass=JITParticle,
                             lon=coords['Longitude'],
                             lat=coords['Latitude'],
                             depth=depth_arg,
                             time=simulation_start)
                            
output_file = pset.ParticleFile(name="/nethome/manra003/analysis/dispersion/simulations/Benguela_0625_401x257_{0}01-31_{1}_{2}z.zarr".format(mon_name, year, r_depth), 
                                outputdt=timedelta(days=1))

pset.execute(AdvectionRK4,                
             runtime=timedelta(days=30),
             dt=300,                       
             output_file=output_file,
             recovery={ErrorCode.ErrorOutOfBounds: delete_particle})

output_file.close()

