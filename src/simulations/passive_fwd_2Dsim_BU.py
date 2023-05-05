import netCDF4 as nc
from glob import glob
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4
from parcels.tools.statuscodes import ErrorCode
from datetime import timedelta, datetime, date
import numpy as np
import sys

args = sys.argv
assert len(args) == 3
year = np.int32(args[1])
month_num = np.int32(args[2])
mon_name = date(1900, month_num, 1).strftime('%b')

data_path = '/storage/shared/oceanparcels/input_data/NEMO16_CMCC/'
mesh_mask = data_path + 'GLOB16L98_mesh_mask_atlantic.nc'

simulation_start = datetime(year, month_num, 1, 12, 0, 0)
simulation_end = datetime(year, month_num, 31, 12, 0, 0)

ufiles =  sorted(glob(data_path + 'ROMEO.01_1d_uo_{0}{1}*_U.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))))

vfiles =  sorted(glob(data_path + 'ROMEO.01_1d_vo_{0}{1}*_V.nc'.\
                      format(simulation_start.strftime("%Y"), simulation_start.strftime("%m"))))

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': vfiles}}

variables = {'U': 'uo',
             'V': 'vo'}

dimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}

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

pset = ParticleSet.from_list(fieldset=fieldset, 
                             pclass=JITParticle,
                             lon=coords['Longitude'],
                             lat=coords['Latitude'],
                             time=simulation_start)
                            
output_file = pset.ParticleFile(name="/nethome/manra003/analysis/dispersion/simulations/Benguela_0625_401x257_{0}01-31_{1}.zarr".format(mon_name,year), 
                                outputdt=timedelta(days=1))

pset.execute(AdvectionRK4,                
             runtime=timedelta(days=30),
             dt=300,                       
             output_file=output_file,
             recovery={ErrorCode.ErrorOutOfBounds: delete_particle})

output_file.close()

